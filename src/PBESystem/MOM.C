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

#include <boost/math/distributions/gamma.hpp>
#include <cassert>

#include "MOM.H"
#include "addToRunTimeSelectionTable.H"

#include "fvScalarMatrix.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcAverage.H"
#include "fvm.H"
#include "Integrator.H"
#include "mathematicalConstants.H"
#include "Utility.H"

namespace Foam
{
namespace PBEMethods
{


defineTypeNameAndDebug(MOM, 0);
addToRunTimeSelectionTable(PBEMethod, MOM, dictionary);
// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
using constant::mathematical::pi;

MOM::MOM
(
    const dictionary& pbeProperties,
    const phaseModel& phase
)
:
    PBEMethod(pbeProperties, phase),
    MOMDict_(pbeProperties.subDict("MOMCoeffs")),
    dispersedPhase_(phase),
    mesh_(dispersedPhase_.U().mesh()),
    moments_(),
    d_
    (
        IOobject
        (
            "diameter",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("diameter", dimLength, 0.0)
    ) ,
    gamma_k_
    (
        IOobject
        (
            "gamma_k",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("gamma_k", dimless, 0.0)
    ),
    gamma_theta_
    (
        IOobject
        (
            "gamma_theta",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("gamma_theta", dimVolume, 0.0)
    ),
    gamma_c0_
    (
        IOobject
        (
            "gamma_c0",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("c0", dimless, 0.0)
    )
{
    for (std::size_t i = 0; i<3; ++i){
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

    d_ = pow(6.0 / pi * moments_[1] / moments_[0], 1.0 / 3.0);
}

MOM::~MOM()
{

}
void MOM::correct()
{
    Info<< "updating size moments" << endl;
    gamma_c0_ = moments_[0];
    gamma_k_ = pow(moments_[1] , 2)
        / (moments_[2] * moments_[0] - pow(moments_[1],2)
           + dimensionedScalar("small", dimensionSet(0, 6, 0, 0, 0), SMALL));

    gamma_theta_ =
        (moments_[0] * moments_[2] - pow(moments_[1], 2))
        /
        max(
            moments_[0] * moments_[1],
            dimensionedScalar("small", dimensionSet(0,3,0,0,0), SMALL)
        );
    
    Info<< "gamma parameters:" << endl;
    printAvgMaxMin(mesh_, gamma_c0_);
    printAvgMaxMin(mesh_, gamma_k_);
    printAvgMaxMin(mesh_, gamma_theta_);

    std::vector<volScalarField> mSources_;

    for (std::size_t i = 0; i<moments_.size(); ++i){
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

    d_ = pow(6.0 / pi * moments_[1] / moments_[0], 1.0/3.0);

    printAvgMaxMin(mesh_, d_);
}

const volScalarField MOM::d() const
{
    return d_;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> MOM::momentSourceTerm(label momenti)
{
    volScalarField bS = breakupSourceTerm(momenti);
    volScalarField cS = coalescenceSourceTerm(momenti);
    Info<< "breakup (moment " << momenti << ") "
        << bS.weightedAverage(mesh_.V()).value() << endl;
    Info<< "coalescence (moment " << momenti << ") "
        << cS.weightedAverage(mesh_.V()).value() << endl;
    return cS + bS;// breakupSourceTerm(momenti);
}

using GammaDistribution = boost::math::gamma_distribution<>;
using boost::math::pdf;

tmp<volScalarField> MOM::coalescenceSourceTerm(label momenti)
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
        GammaDistribution gamma(gamma_k_[celli], gamma_theta_[celli]);

        auto m0 = moments_[0][celli];

        auto deathIntegrand = [&](double xi){

            auto deathInnerIntegrand = [&](double xi_prime){

                return coalescence_->S(xi, xi_prime).value() *
                       m0 * pdf(gamma, xi_prime);
            };

            return pow(xi, momenti) * m0 * pdf(gamma, xi) *
                   integrate(deathInnerIntegrand, 0.);
        };

        auto birthIntegrand = [&](double xi){

            auto birthInnerIntegrand = [&](double xi_prime){

                return coalescence_->S(xi, xi_prime).value() *
                       pow(xi_prime + xi, momenti) *
                       m0 * pdf(gamma, xi_prime);
            };

            return m0 * pdf(gamma, xi) *
                   integrate(birthInnerIntegrand, 0.);
        };

        toReturn[celli] = integrate(birthIntegrand, 0.) / 2. -
                          integrate(deathIntegrand, 0.);
    }

    return tmp<volScalarField>( new volScalarField(toReturn));
}

tmp<volScalarField> MOM::breakupSourceTerm(label momenti)
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

    forAll(dispersedPhase_, celli)
    {
        GammaDistribution gamma(gamma_k_[celli], gamma_theta_[celli]);

        auto m0 = moments_[0][celli];

        auto breakupSourceIntegrand = [&](double xi){
            scalar breakupDeath = breakup_->S(xi).value() * m0 * pdf(gamma, xi);

            auto breakupBirthIntegrand = [&](double xi_prime){
                return daughterParticleDistribution_->beta(xi, xi_prime).value()
                    * breakup_->S(xi_prime).value() * m0 * pdf(gamma, xi_prime);
            };

            double breakupBirth = integrate(breakupBirthIntegrand, xi);

            return pow(xi, momenti) * (breakupBirth - breakupDeath);
        };
        
        toReturn[celli] = integrate(breakupSourceIntegrand, 0.);
    }

    return tmp<volScalarField>( new volScalarField(toReturn));
}


}//end namespace PBEMethods
}//end namespace Foam

// ************************************************************************* //
