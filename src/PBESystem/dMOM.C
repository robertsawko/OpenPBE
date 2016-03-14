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

#include <string>
#include <sstream>
#include <cmath>

#include <cassert>

#include "dMOM.H"
#include "addToRunTimeSelectionTable.H"

#include "fvScalarMatrix.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcAverage.H"
#include "fvm.H"
#include "Integrator.H"
#include "mathematicalConstants.H"
#include "Utility.H"
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
      mesh_(dispersedPhase_.U().mesh()), Sgammas_(NO_OF_MOMENTS),
      breakupSource_(NO_OF_MOMENTS), coalescenceSource_(NO_OF_MOMENTS),
      log_d_bar_(IOobject("logd_bar",
                          mesh_.time().timeName(),
                          mesh_,
                          IOobject::NO_READ,
                          IOobject::AUTO_WRITE),
                 mesh_,
                 dimensionedScalar("log_dbar", dimless, 0.0)),
      sigma_hat_(IOobject("sigma_hat",
                          mesh_.time().timeName(),
                          mesh_,
                          IOobject::NO_READ,
                          IOobject::AUTO_WRITE),
                 mesh_,
                 dimensionedScalar("sigma_hat", dimless, 0.0)),
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
      k_br_(0.2),    // Following paper [2], comment below equation [30]
      Ca_cr_(
          0.534),   // This follows Bruijn model for aqueous potassium carbonate
                    // solution in kersone appropriate for S&A case.
      k_cl_1_(1.0), // Paper [1], comment after eq 30
      consistency_( // TODO: Dimensional consistency is something not discussed
          "1m",     // in either of the papers. log(d) as such cannot be a
          dimLength, // quantity as log is transcendental and d dimensional.
          1.0),      // Therefore, there must be some sort of consistency scale.
      interfacialTension_( // Must be separate for 2phase Euler!
          readScalar(dMOMDict_.lookup("interfacialTension"))),
      minDiameter_("min(d)",
                   dimLength,
                   dMOMDict_.lookupOrDefault<scalar>("minDiameter", 1e-9)),
      k_coll_(sqrt(8 * pi / 3)),
      // Hamaker constant; given in paper [2] after in the comment to eq 46
      A_H_(5e-21)
{
    forAll(Sgammas_, i) {
        Sgammas_.set(i,
                     new volScalarField(IOobject("S" + std::to_string(i),
                                                 mesh_.time().timeName(),
                                                 mesh_,
                                                 IOobject::MUST_READ,
                                                 IOobject::AUTO_WRITE),
                                        mesh_));
        breakupSource_.set(
            i,
            new volScalarField(
                IOobject("breakup" + std::to_string(i) + "Source",
                         mesh_.time().timeName(),
                         mesh_,
                         IOobject::NO_READ,
                         IOobject::AUTO_WRITE),
                mesh_,
                dimensionedScalar(
                    "0", Sgammas_[i].dimensions() / dimTime, 0.0)));
        coalescenceSource_.set(
            i,
            new volScalarField(
                IOobject("coalescence" + std::to_string(i) + "Source",
                         mesh_.time().timeName(),
                         mesh_,
                         IOobject::NO_READ,
                         IOobject::AUTO_WRITE),
                mesh_,
                dimensionedScalar(
                    "0", Sgammas_[i].dimensions() / dimTime, 0.0)));
    }
    printAvgMaxMin(mesh_, d());
}

dMOM::~dMOM() {}

using boost::math::normal;
using boost::math::pdf;

void dMOM::correct() {
    // Creating aliases to shorten notation
    const volScalarField &S0 = Sgammas_[0];
    const volScalarField &S1 = Sgammas_[1];
    const volScalarField &S2 = Sgammas_[2];
    // Calculate moments (note equation (2) in paper [1])
    const volScalarField m1 = S1 / S0;
    const volScalarField m2 = S2 / S0;
    // Non-logarithmized variance
    // TODO: Investigate why variance become negative
    const volScalarField v = mag(m2 - pow(m1, 2));

    Info << "updating size moments" << endl;
    // Moment inversion is given analytically (see )
    log_d_bar_ = log(m1 / consistency_ / sqrt(1.0 + v / pow(m1, 2)));
    sigma_hat_ = sqrt(log(1.0 + v / pow(m1, 2)));

    Info << "Lognormal parameters:" << endl;
    printAvgMaxMin(mesh_, S0);
    printAvgMaxMin(mesh_, log_d_bar_);
    printAvgMaxMin(mesh_, sigma_hat_);

    std::vector<volScalarField> mSources_;
    const volScalarField &epsilon =
        dispersedPhase_.U().mesh().lookupObject<volScalarField>(
            "epsilon." + dispersedPhase_.name());
    const scalar &sigma = interfacialTension_;

    forAll(dispersedPhase_, celli) {
        // **** BREAKUP CHARACTERISTICS CALCULATIONS ****
        // Turbulent length scale
        auto L_k = pow(
            pow(phase_.otherPhase().nu()()[celli], 4) / epsilon[celli], 0.25);

        auto shear_rate =
            sqrt(phase_.otherPhase().rho()[celli] * epsilon[celli] /
                 phase_.otherPhase().mu()()[celli]);
        // Critical diameter for viscous breakup
        auto d_cr_viscous =
            2 * sigma * Ca_cr_ / phase_.otherPhase().mu()()[celli] / shear_rate;

        // Critical diameter for inertial breakup
        scalar d_cr_inertia =
            (1 + C_alpha_) *
            pow(2 * sigma * We_cr_ / phase_.rho()[celli], 3.0 / 5.0) *
            pow(epsilon[celli], -2.0 / 5.0);

        // Molecular viscosity ratio
        auto lambda = phase_.mu()()[celli] / phase_.otherPhase().mu()()[celli];
        auto ftau =
            exp(1.43 + 0.267 * log(lambda) - 0.023 * pow(log(lambda), 2));
        auto tau_viscous_constant =
            (phase_.mu()()[celli] * shear_rate) / sigma * ftau;
        // Physical constants for inertial breakup
        auto tau_inertia_constant =
            2 * pi * k_br_ * sqrt((3 * phase_.rho()[celli] +
                                   2 * phase_.otherPhase().rho()[celli]) /
                                  (192 * sigma));

        // Modified distributions parameters
        auto &std = sigma_hat_[celli];

        // **** COALESCENCE CHARACTERISTICS CALCULATIONS ****
        // Paper[1] equation (5)
        auto S3 = phase_[celli] * 6.0 / pi;
        scalar d_eq[NO_OF_MOMENTS], u_rel[NO_OF_MOMENTS];
        // Critical film thickness
        scalar h_cr[NO_OF_MOMENTS];
        // Interaction force
        scalar P_coal[NO_OF_MOMENTS];

        // TODO: My guess is that here a switching MUST occur between inertial
        // and viscous breakup.
        for (int gamma=0; gamma < NO_OF_MOMENTS; ++gamma){
            auto Sgamma = Sgammas_[gamma][celli];
            // Paper[1] equation (3)
            auto d3gamma = pow(S3 / Sgamma, 1 / (3.0 - gamma)); 
            d_eq[gamma] = k_cl_1_ * d3gamma;
            u_rel[gamma] = shear_rate * d_eq[gamma];
            // Critical film thickness eq (46) in paper [2]
            auto h_cr = pow(A_H_ * d_eq[gamma] / (25.0 * pi * sigma), 1.0 / 3.0);
            auto F_i =
                3.0 * pi / 2.0 * phase_.otherPhase().mu()()[celli] *
                shear_rate * pow(d_eq[gamma], 2);
            // Finally, ladies and gentlemen, I present to you:
            // Drainage model 3!
            auto td =
                pi * phase_.mu()()[celli] * sqrt(F_i) / (2 * h_cr) *
                pow(d_eq[gamma] / (4 * pi * sigma), 3.0 / 2.0);
            // Equation (40)
            P_coal[gamma] = exp(-td * shear_rate);
        }

        for (label gamma = 0; gamma < NO_OF_MOMENTS; ++gamma) {
            auto mu_viscous = log_d_bar_[celli] + (gamma - 1) * std;
            auto mu_inertial = log_d_bar_[celli] + (gamma - 1.5) * std;
            auto x_crit = log(d_cr_inertia);

            normal inertial(mu_inertial, std);
            normal viscous(mu_viscous, std);

            // exp expression are correction factor that results from the
            // x = log(d) substitution and recasting the integrals as erfc.
            // TODO: Needs more doc!
            breakupSource_[gamma][celli] =
                // DPD constant contribution; the remaining d^gamma goes into
                // correction factor
                (pow(Nf_, (3.0 - gamma) / 3.0) - 1.0) *
                (
                    // Viscous breakup
                    1 / tau_viscous_constant *
                        (cdf(viscous, log(L_k)) -
                         cdf(viscous, log(d_cr_viscous))) *
                        exp(pow(gamma * std, 2) + 2.0 * gamma * mu_viscous -
                            2.0 * gamma * pow(std, 2) - 2.0 * mu_viscous +
                            pow(sigma, 2)) +
                    // Inertial breakup contribution
                    1 / tau_inertia_constant *
                        // TODO: does x_crit need limiting?
                        (1 - cdf(inertial, max(x_crit, log(L_k)))) *
                        exp(pow(gamma * std, 2) + 2.0 * gamma * mu_inertial -
                            3.0 * gamma * pow(std, 2) - 3.0 * mu_inertial +
                            9.0 / 4.0 * pow(sigma, 2)));

            coalescenceSource_[gamma][celli] =
                (pow(2.0, gamma / 3.0) - 2.0) *
                pow(6.0 * phase_[celli] / pi, 2.0) *
                k_coll_ * u_rel[gamma] *
                P_coal[gamma] *
                pow(d_eq[gamma], gamma - 4);
        }
    }
    for (auto bs : breakupSource_)
        printAvgMaxMin(mesh_, bs);
    for (auto cs : coalescenceSource_)
        printAvgMaxMin(mesh_, cs);
    // TODO: is there a better way to assure that source terms on boundaries
    // are equal to 0?
    // maybe using internalField() instead of the whole field somewhere?
    // how source terms are handled in combustion/chemistry/turbulence codes in
    // OF

    forAll(Sgammas_, i) {
        fvScalarMatrix mEqn(
            fvm::ddt(Sgammas_[i]) +
            fvm::div(dispersedPhase_.phi(), Sgammas_[i])
            //+ fvm::div(limitedFlux_,moments_[i])
            ==
            breakupSource_[i] + coalescenceSource_[i]);
        mEqn.relax();
        mEqn.solve();
    }

    for (auto &moment : Sgammas_) {
        moment =
            max(moment,
                dimensionedScalar(moment.name(), moment.dimensions(), SMALL));
        // TODO: print a warning message
        printAvgMaxMin(mesh_, moment);
    }

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
    return max(
        6.0 / pi * dispersedPhase_ / Sgammas_[2],            // Equatio (6) from paper [1]
        Sgammas_[1] / Sgammas_[0]); // This could be user defined
                                    // mean diameter, instead.
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

} // end namespace PBEMethods
} // end namespace Foam

// ************************************************************************* //
