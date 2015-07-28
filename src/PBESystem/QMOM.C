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

#include "QMOM.H"
#include "addToRunTimeSelectionTable.H"

#include "fvScalarMatrix.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcAverage.H"
#include "fvm.H"
#include "Integrator.H"
#include "mathematicalConstants.H"
#include "Utility.H"
#include "MomentInversion.H"

namespace Foam
{
namespace PBEMethods
{


defineTypeNameAndDebug(QMOM, 0);
addToRunTimeSelectionTable(PBEMethod, QMOM, dictionary);
// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
using constant::mathematical::pi;

QMOM::QMOM
(
        const dictionary& pbeProperties,
        const phaseModel& phase
        )
    :
      PBEMethod(pbeProperties, phase),
      QMOMDict_(pbeProperties.subDict("QMOMCoeffs")),
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
    )
{
    for (std::size_t i = 0; i<3; ++i){
        moments_.emplace_back(
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

QMOM::~QMOM()
{

}
void QMOM::correct()
{
    std::vector<volScalarField> mSources_;

    for (std::size_t i = 0; i < moments_.size(); ++i){
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
    //clamp d_ here?


    printAvgMaxMin(mesh_, d_);
}

const volScalarField QMOM::d() const
{
    return d_;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> QMOM::momentSourceTerm(label momenti)
{
    volScalarField bS = breakupSourceTerm(momenti);
    volScalarField cS = coalescenceSourceTerm(momenti);
    Info<< "breakup (moment " << momenti << ") "
        << bS.weightedAverage(mesh_.V()).value() << endl;
    Info<< "coalescence (moment " << momenti << ") "
        << cS.weightedAverage(mesh_.V()).value() << endl;
    return cS + bS;// breakupSourceTerm(momenti);
}


tmp<volScalarField> QMOM::coalescenceSourceTerm(label momenti)
{
    volScalarField toReturn
    (
        IOobject
        (
            "Sc",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar
        (
            "Sc",
            pow(dimVolume, momenti) / dimTime,
            0
        )
    );

    forAll(dispersedPhase_, celli)
    {
        Eigen::VectorXd momentVector(moments_.size());
        for (int i=0; i<momentVector.size(); ++i)
            momentVector[i] = moments_[i][celli];

        auto quadrature = wheeler_inversion(momentVector);
        int N = quadrature.abcissas.size();

        double death{0.}, birth{0.};

        for (int i=0; i<N; ++i){
            double xi_i = quadrature.abcissas[i];
            double innerDeathTerm{0.}, innerBirthTerm{0.};

            for (int j=0; j<N; ++j){
                double xi_j = quadrature.abcissas[j];

                innerDeathTerm += quadrature.weights[j] *
                             coalescence_->S(xi_i, xi_j).value();

                innerBirthTerm += quadrature.weights[j] *
                        pow(xi_i + xi_j, momenti) *
                        coalescence_->S(xi_i, xi_j).value();
            }

            death += pow(xi_i, momenti) *
                     quadrature.weights[i] *
                     innerDeathTerm;

            birth += quadrature.weights[i] *
                     innerBirthTerm;
        }
        toReturn[celli] = birth/2. - death;
    }

    return tmp<volScalarField>( new volScalarField(toReturn));

}

tmp<volScalarField> QMOM::breakupSourceTerm(label momenti)
{

    //value of the integral
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
        Eigen::VectorXd momentVector(moments_.size());
        for (int i=0; i<momentVector.size(); ++i)
            momentVector[i] = moments_[i][celli];

        auto quadrature = wheeler_inversion(momentVector);
        int N = quadrature.abcissas.size();

        double death = 0.;

        for (int i=0; i<N; ++i){
            double xi_i = quadrature.abcissas[i];

            death += quadrature.weights[i] *
                    pow(xi_i, momenti) *
                    breakup_->S(xi_i).value();
        }

        auto birthIntegrand = [&](double xi){
            double result = 0.;

            for (int i=0;i<N; ++i){
                double xi_i = quadrature.abcissas[i];

                result += quadrature.weights[i] *
                        breakup_->S(xi_i).value() *
                        daughterParticleDistribution_->
                        beta(xi,xi_i).value();
            }

            return pow(xi, momenti)*result;
        };

        auto birth = integrate(birthIntegrand, 0.);

        toReturn[celli] = birth - death;

    }

    return tmp<volScalarField>( new volScalarField(toReturn));
}


}//end namespace PBEMethods
}//end namespace Foam

// ************************************************************************* //
