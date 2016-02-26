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

namespace Foam
{
namespace PBEMethods
{

defineTypeNameAndDebug(dMOM, 0);
addToRunTimeSelectionTable(PBEMethod, dMOM, dictionary);
// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
using constant::mathematical::pi;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dMOM::dMOM
(
    const dictionary& pbeProperties,
    const phaseModel& phase
)
:
    PBEMethod(pbeProperties, phase),
    dMOMDict_(pbeProperties.subDict("dMOMCoeffs")),
    dispersedPhase_(phase),
    mesh_(dispersedPhase_.U().mesh()),
    moments_(),
    log_d_bar_
    (
        IOobject
        (
            "log(d bar)",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("log(d bar)", dimless, 0.0)
    ),
    sigma_hat_
    (
        IOobject
        (
            "sigma hat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("sigma_hat", dimLength, 0.0)
    ),
    number_density_
    (
        IOobject
        (
            "n",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("n", dimless, 0.0)
    ),
    Nf_(2.0),       // Both papers take 2. TODO: Could make user slectible
    C_alpha_(0.5),  // TODO: This is not specified in any of the papers!
    We_cr_(0.31),   // Following paper [2], which cites Yo and Morel (2001)
    k_br_(0.2),     // Following paper [2], comment below equation [30]
    Ca_cr_(0.534)   // This follows Bruijn model for aqueous potassium carbonate
                    // solution in kersone appropriate for S&A case.
{
    for (std::size_t i = 0; i < 3; ++i){
        moments_.emplace_back
            (
                IOobject
                (
                    "m" + std::to_string(i),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            );
    }
}

dMOM::~dMOM()
{
}

void dMOM::correct()
{
    // Creating aliases to shorten notation
    const volScalarField& m1 = moments_[1];
    const volScalarField& m2 = moments_[2];
    // Non-logarithmized variance
    const volScalarField v = m2 - pow(m1, 2);

    Info<< "updating size moments" << endl;
    number_density_ = moments_[0];
    // Moment inversion is given analytically (see )
    log_d_bar_ = log(m1 / sqrt(1 + v / pow(m1, 2)));
    sigma_hat_ = sqrt(log(1.0 + v / pow(m1, 2)));

    Info<< "Lognormal parameters:" << endl;
    printAvgMaxMin(mesh_, number_density_);
    printAvgMaxMin(mesh_, log_d_bar_);
    printAvgMaxMin(mesh_, sigma_hat_);

    std::vector<volScalarField> mSources_;

    for (std::size_t i = 0; i<moments_.size(); ++i) {
        mSources_.emplace_back
                (
                    IOobject
                    (
                        "m" + std::to_string(i) + "Source",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                        ),
                    momentSourceTerm(i)
                    );
    }

    Info<< "moment sources:" << endl;
    //TODO: is there a better way to assure that source terms on boundaries
    //are equal to 0?
    //maybe using internalField() instead of the whole field somewhere?
    //how source terms are handled in combustion/chemistry/turbulence codes in
    //OF
    for (auto& mSource : mSources_)
    {
        mSource.boundaryField() = 0;

        //fix for nan values
        forAll(mesh_.C(), celli)
        {
            auto& mSource_i = mSource[celli];
            if ( std::isnan(mSource_i) )
            {
                mSource_i = 0.0;
            }
        }

        printAvgMaxMin(mesh_, mSource);
    }

    for (std::size_t i = 0; i<moments_.size(); ++i)
    {
        fvScalarMatrix mEqn
        (
            fvm::ddt(moments_[i])
            + fvm::div(dispersedPhase_.phi(),moments_[i])
            //+ fvm::div(limitedFlux_,moments_[i])
            ==
            mSources_[i]
        );
        mEqn.relax();
        mEqn.solve();
    }

    for (auto& moment : moments_)
    {
        moment = max(
            moment,
            dimensionedScalar(moment.name(), moment.dimensions(), SMALL)
        );
        //TODO: print a warning message
        printAvgMaxMin(mesh_, moment);
    }

    printAvgMaxMin(mesh_, d());
}

const volScalarField dMOM::d() const
{
    // Diameter is assumed to be Sauter mean by following equation (6) in [1]
    return 6.0 / pi * dispersedPhase_ / moments_[2];
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> dMOM::momentSourceTerm(label momenti)
{
    volScalarField bS = breakupSourceTerm(momenti);
    volScalarField cS = coalescenceSourceTerm(momenti);
    Info<< "breakup (moment " << momenti << ") "
        << bS.weightedAverage(mesh_.V()).value() << endl;
    Info<< "coalescence (moment " << momenti << ") "
        << cS.weightedAverage(mesh_.V()).value() << endl;
    return cS + bS;// breakupSourceTerm(momenti);
}

tmp<volScalarField> dMOM::coalescenceSourceTerm(label momenti)
{
    volScalarField toReturn
    (
        IOobject
        (
            "Sbc",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar
        (
            "Sbc",
            pow(dimVolume, momenti) /dimTime,
            0
        )
    );

    forAll(dispersedPhase_, celli)
    {
        /*
        GammaDistribution gamma(gamma_k_[celli], gamma_theta_[celli]);

        auto m0 = moments_[0][celli];

        auto deathIntegrand = [&](double xi){

            auto deathInnerIntegrand = [&](double xi_prime){

                return coalescence_->S(xi, xi_prime, celli) *
                       m0 * pdf(gamma, xi_prime);
            };

            return pow(xi, momenti) * m0 * pdf(gamma, xi) *
                   integrate(deathInnerIntegrand, 0.);
        };

        auto birthIntegrand = [&](double xi){

            auto birthInnerIntegrand = [&](double xi_prime){

                return coalescence_->S(xi, xi_prime, celli) *
                       pow(xi_prime + xi, momenti) *
                       m0 * pdf(gamma, xi_prime);
            };

            return m0 * pdf(gamma, xi) *
                   integrate(birthInnerIntegrand, 0.);
        };

        toReturn[celli] = integrate(birthIntegrand, 0.) / 2. -
                          integrate(deathIntegrand, 0.);
                          */
    }

    return tmp<volScalarField>( new volScalarField(toReturn));
}

using boost::math::normal;
using boost::math::pdf;

tmp<volScalarField> dMOM::breakupSourceTerm(label momenti)
{
    volScalarField toReturn
    (
        IOobject
        (
            "Sbr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar
        (
            "Sbr", 
            pow(dimVolume, momenti) /dimTime,
            0
        ) 
    );

    const volScalarField& epsilon =
        dispersedPhase_.U().mesh().lookupObject<volScalarField>(
            "epsilon."+ dispersedPhase_.name());
    const scalar sigma = phase_.fluid().sigma().value();
    const label& gamma = momenti;

    forAll(dispersedPhase_, celli)
    {
        // Turbulent length scale
        auto L_k = pow(
            pow(phase_.otherPhase().nu()()[celli], 4) / epsilon[celli],
            0.25);

        auto shear_rate = sqrt(
            phase_.otherPhase().rho()[celli] * epsilon[celli] /
            phase_.otherPhase().mu()()[celli]);
        // Critical diameter for viscous breakup
        auto d_cr_viscous =
            2 * sigma * Ca_cr_ / phase_.otherPhase().mu()()[celli] / shear_rate;

        // Critical diameter for inertial breakup
        scalar d_cr_inertia =
            (1 + C_alpha_) *
            pow(2 * sigma * We_cr_ / phase_.rho()[celli], 3.0 / 5.0) *
            pow(epsilon[celli], -2.0 / 5.0);
        
        // Physical constants for inertial breakup
        auto tau_inertia_constant = 2 * pi * k_br_ * sqrt(
            (3 * phase_.rho()[celli] + 2 * phase_.otherPhase().rho()[celli]) /
            (192 * sigma));


        // Modified distribution parameters
        auto& std = sigma_hat_[celli];
        auto mu = log_d_bar_[celli] + (gamma - 1.5) * std;
        auto x_crit = log(d_cr_inertia);

        normal distribution(mu, std);

        toReturn[celli] = 
            // DPD constant contribution; the remaining d^gamma goes into
            // correction factor
            (pow(Nf_, (3 - gamma) / 3) - 1.0) *
            (
                // Viscous breakup
                // ??? +
                // Inertial breakup contribution
                1 / tau_inertia_constant * 
                (1 - cdf(distribution, max(x_crit, log(L_k)))) * //TODO: does
                //  x_crit need limiting? This is the correction factor that
                // results from the x = log(d) substitution and the follow up
                // term rearrangment. TODO: Needs more doc!
                exp( 
                    pow(gamma * std, 2) +
                    2.0 * gamma * mu - 
                    3.0 * gamma * pow(std, 2) - 
                    3.0 * mu +
                    9.0 / 4.0 * pow(sigma, 2))
            );
    }

    return tmp<volScalarField>( new volScalarField(toReturn));
}


}//end namespace PBEMethods
}//end namespace Foam

// ************************************************************************* //
