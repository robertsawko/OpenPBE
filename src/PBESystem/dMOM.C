/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include <cmath>
#include <sstream>
#include <string>

#include <cassert>

#include "addToRunTimeSelectionTable.H"
#include "dMOM.H"

#include "Integrator.H"
#include "Utility.H"
#include "fvScalarMatrix.H"
#include "fvcAverage.H"
#include "fvcFlux.H"
#include "fvcSnGrad.H"
#include "fvm.H"
#include "mathematicalConstants.H"
#include "twoPhaseSystem.C"
#include <boost/math/distributions/normal.hpp> // for normal_distribution

namespace Foam {
namespace PBEMethods {

defineTypeNameAndDebug(dMOM, 0);
addToRunTimeSelectionTable(PBEMethod, dMOM, dictionary);
// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
using constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dMOM::dMOM(const dictionary &pbeProperties, const phaseModel &phase)
    : PBEMethod(pbeProperties, phase),
      dMOMDict_(pbeProperties.subDict("dMOMCoeffs")), dispersedPhase_(phase),
      mesh_(dispersedPhase_.U().mesh()), S2_(IOobject("S2",
                                                      mesh_.time().timeName(),
                                                      mesh_,
                                                      IOobject::MUST_READ,
                                                      IOobject::AUTO_WRITE),
                                             mesh_),
      breakupSource_(IOobject("breakupSource",
                              mesh_.time().timeName(),
                              mesh_,
                              IOobject::NO_READ,
                              IOobject::AUTO_WRITE),
                     mesh_,
                     dimensionedScalar("0", S2_.dimensions() / dimTime, 0.0)),
      coalescenceSource_(
          IOobject("coalescenceSource",
                   mesh_.time().timeName(),
                   mesh_,
                   IOobject::NO_READ,
                   IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar("0", S2_.dimensions() / dimTime, 0.0)),
      sigma2_hat_(IOobject("sigma2_hat",
                           mesh_.time().timeName(),
                           mesh_,
                           IOobject::NO_READ,
                           IOobject::AUTO_WRITE),
                  mesh_,
                  dimensionedScalar("sigma2_hat", dimless, 0.0)),
      number_density_(IOobject("n",
                               mesh_.time().timeName(),
                               mesh_,
                               IOobject::NO_READ,
                               IOobject::AUTO_WRITE),
                      mesh_,
                      dimensionedScalar("n", dimless, 0.0)),
      Nf_(2.0),      // Both papers take 2. TODO: Could make user slectible
      C_alpha_(0.5), // TODO: This is not specified in any of the papers!
      We_cr_(0.31),  // Following paper [2], which cites Yo and Morel (2001)
      We0_(0.8 * We_cr_), // Paper [1] equation 51
      k_br_(0.2),         // Paper [2], comment below equation [30]
      Ca_cr_(0.534), // Following Bruijn model for aqueous potassium carbonate
                     // solution in kersone appropriate for S&A case.
      k_cl_1_(1.0),  // Paper [1], comment after eq 30
      k_cl_2_(12.7), // Paper [1], given by eq 50
      consistency_(  // TODO: Dimensional consistency is something not discussed
          "1m",      // in either of the papers. log(d) as such cannot be a
          dimLength, // dimensional quantity as log is transcendental and d is
          1.0),      // dimensional. There must be some sort of scaling.
      interfacialTension_( // Must be separate for 2phase Euler!
          readScalar(dMOMDict_.lookup("interfacialTension"))),
      minDiameter_("min(d)",
                   dimLength,
                   dMOMDict_.lookupOrDefault<scalar>("minDiameter", 1e-5)),
      k_collv_(sqrt(8 * pi / 3)),  // Paper [1] equation 36
      k_colli_(sqrt(2 * pi / 15)), // Paper [1] equation 46
      A_H_(5e-21), // Hamaker constant; given in paper [2] after in the comment
                   // to eq 46
      Fcl_(dMOMDict_.lookupOrDefault<scalar>("calibration", 1.0)) {

    printAvgMaxMin(mesh_, d());
}

dMOM::~dMOM() {}

using boost::math::normal;
using boost::math::pdf;

void dMOM::correct() {
    // Creating aliases to shorten notation
    // const auto S0 = 6.0 / pi * phase_ * pow(d(), -3.0);
    const volScalarField S3(
        "S3,", max(phase_ * 6.0 / pi, dimensionedScalar("s", dimless, 1e-3)));
    const auto d32 = d();

    const volScalarField &epsilon =
        dispersedPhase_.U().mesh().lookupObject<volScalarField>(
            "epsilon." + dispersedPhase_.name());

    Info << "updating size moments" << endl;
    // Inversion procedure
    // Moment inversion is given analytically (we allow minimum 1% variation)
    sigma2_hat_ = max(2.0 / 5.0 * log(S3 / S2_ * pow(d32, -1.0)), 0.001);
    Info << sigma2_hat_.dimensions() << " " << endl;
    // Inversion procedure
    const volScalarField S0("S0", exp(log(S2_ * d32) - sigma2_hat_));

    printAvgMaxMin(mesh_, S0);
    printAvgMaxMin(mesh_, S2_);
    printAvgMaxMin(mesh_, S3);
    Info << "Lognormal parameters:" << endl;
    // printAvgMaxMin(mesh_, S0);
    printAvgMaxMin(mesh_, sigma2_hat_);

    const scalar &sigma = interfacialTension_;

    forAll(dispersedPhase_, celli) {
        const scalar &muc = phase_.otherPhase().mu()()[celli];
        // **** BREAKUP CHARACTERISTICS CALCULATIONS ****
        // Turbulent length scale
        auto L_k = pow(
            pow(phase_.otherPhase().nu()()[celli], 4) / epsilon[celli], 0.25);

        auto shear_rate =
            sqrt(phase_.otherPhase().rho()[celli] * epsilon[celli] / muc);
        // Critical diameter for viscous breakup
        auto d_cr_viscous = 2 * sigma * Ca_cr_ / (muc * shear_rate);

        // Critical diameter for inertial breakup
        scalar d_cr_inertia =
            (1 + C_alpha_) *
            pow(2 * sigma * We_cr_ / phase_.rho()[celli], 3.0 / 5.0) *
            pow(epsilon[celli], -2.0 / 5.0);

        // Molecular viscosity ratio
        auto lambda = phase_.mu()()[celli] / muc;
        // Constants takens from Hill's/Brujin
        auto ftau =
            exp(1.43 + 0.267 * log(lambda) - 0.023 * pow(log(lambda), 2));
        auto tau_viscous_constant = muc / sigma * ftau;
        // Physical constants for inertial breakup
        auto tau_inertia_constant =
            2.0 * pi * k_br_ * sqrt((3.0 * phase_.rho()[celli] +
                                     2.0 * phase_.otherPhase().rho()[celli]) /
                                    (192.0 * sigma));

        // Modified distributions parameters
        auto &var = sigma2_hat_[celli];

        // **** COALESCENCE CHARACTERISTICS CALCULATIONS ****
        // Paper[1] equation (5)

        // TODO: My guess is that here a switching MUST occur between inertial
        // and viscous breakup.
        // Paper[1] equation (3) for gamma = 2
        scalar d_eq = k_cl_1_ * d32[celli];
        scalar u_relv = shear_rate * d_eq;
        // scalar u_reli = pow(epsilon[celli] * d_eq, 1.0 / 3.0);
        // Critical film thickness eq (46) in paper [2]
        auto h_cr = pow(A_H_ * d_eq / (24.0 * pi * sigma), 1.0 / 3.0);
        auto F_i = 3.0 * pi / 2.0 * phase_.otherPhase().mu()()[celli] *
                   shear_rate * pow(d_eq, 2);
        // Finally, ladies and gentlemen, I present to you:
        // Drainage model 3!
        auto td = pi * phase_.mu()()[celli] * sqrt(F_i) / (2 * h_cr) *
                  pow(d_eq / (4 * pi * sigma), 3.0 / 2.0);
        // Equation (40)
        scalar P_coalv = exp(-td * shear_rate);
        // Info << "P_coalv=" << P_coalv << endl;
        // auto h0 = 8.3 * h_cr; // Paper [1] equation 53
        // Info << "h0=" << h0 << " rho=" << phase_.otherPhase().rho()[celli]
        //<< " sigma=" << interfacialTension_ << " epsilon="
        //<< epsilon[celli] << " We0=" << We0_ << " d_eq=" << d_eq
        //<< " mu2=" << pow(phase_.mu()()[celli], 2.0) << endl;
        // auto Phi_max = 2.0 * pow(h0, 2.0) * phase_.otherPhase().rho()[celli]
        // *
        // interfacialTension_ /
        //(We0_ * pow(phase_.mu()()[celli], 2.0) * d_eq);
        // Info << "Phi_max=" << Phi_max << endl;
        // auto We = phase_.otherPhase().rho()[celli] *
        // pow(epsilon[celli], 2.0 / 3.0) * pow(d_eq, 5.0 / 3.0) /
        // interfacialTension_;
        // Info << "We=" << We << " We0=" << We0_ << " pow="
        //<< pow(k_cl_2_ * (We - We0_) / Phi_max, 2.0) << endl;
        // scalar P_coali =
        // Phi_max / pi *
        // pow(1.0 - pow(k_cl_2_ * (We - We0_) / Phi_max, 2.0), 0.5);
        // Info << "P_coali=" << P_coali << endl;

        // Equationsw with gamma = 2.0
        auto mu_viscous = var;
        auto mu_inertial = 0.5 * var;
        auto t_k = log(L_k / d32[celli]);
        auto t_crv = log(d_cr_viscous / d32[celli]);
        auto t_cri = log(d_cr_inertia / d32[celli]);

        normal inertial(mu_inertial, sqrt(var));
        normal viscous(mu_viscous, sqrt(var));

        // exp expression are correction factor that results from the
        // t = log(d/{\bar d}) substitution and recasting the integrals as erfc.
        // TODO: Needs more doc!
        breakupSource_[celli] =
            // DPD constant contribution; the remaining d^gamma goes into
            // correction factor
            (pow(Nf_, 1.0 / 3.0) - 1.0) *
                (
                    // Viscous breakup
                    1 / tau_viscous_constant *
                    (cdf(viscous, max(t_crv, t_k)) - cdf(viscous, t_crv))) *
                exp(var / 2.0) / d32[celli] +
            // Inertial breakup contribution
            1 / tau_inertia_constant *
                // TODO: does t_cr need limiting?
                (1 - cdf(inertial, max(t_cri, t_k))) * exp(var / 8.0) /
                d32[celli];

        // Note many terms simplify and S2 remains
        coalescenceSource_[celli] = Fcl_ * (pow(2.0, 2.0 / 3.0) - 2.0) *
                                    (k_collv_ * u_relv * P_coalv // viscous
                                     // + k_colli_ * u_reli * P_coali // inertia
                                     ) *
                                    pow(S2_[celli], 2.0) * // contribution from
                                    pow(k_cl_1_, -2.0);    // pow(d_eq, -2)
    }
    printAvgMaxMin(mesh_, breakupSource_);
    printAvgMaxMin(mesh_, coalescenceSource_);
    // TODO: is there a better way to assure that source terms on boundaries
    // are equal to 0?
    // maybe using internalField() instead of the whole field somewhere?
    // how source terms are handled in combustion/chemistry/turbulence codes in
    // OF

    // Paper [1] comment below equation 8: only S2 is solved
    fvScalarMatrix mEqn(fvm::ddt(S2_) + fvm::div(dispersedPhase_.phi(), S2_)
                        //+ fvm::div(limitedFlux_,moments_[i])
                        == breakupSource_ + coalescenceSource_);
    mEqn.relax();
    mEqn.solve();
    S2_ ==
        max(S2_,
            dimensionedScalar("0", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.001));

    printAvgMaxMin(mesh_, d());
}

const volScalarField dMOM::d() const {
    // NOTE: There is some ciruclar logic in all Lo's models regarding the
    // Sauter mean diameter calculation. All of his equations are valid if and
    // only if the drops are exact spheres in which case, the Sauter mean should
    // be exactly equal to mean diameter. A subtle difference comes only from
    // the averaging process. mean(d**3) / mean(d**2) is not the same as
    // mean(d). I do not think it is possible to estimate the difference using
    // Jensen's type arguments.
    //
    // The expression below is used for thresholding the value so that small or
    // zero phase fraction values do not compromise the stability of two fluid
    // model (zero diameter limit).
    return max(6.0 / pi * dispersedPhase_ / S2_, // Equation (6) from paper [1]
               minDiameter_);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

} // end namespace PBEMethods
} // end namespace Foam

// ************************************************************************* //
