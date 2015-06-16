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

#include "MOM.H"
#include "addToRunTimeSelectionTable.H"

#include "PBESystems-internal.H"
#include "fvScalarMatrix.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcAverage.H"
#include "fvm.H"

namespace Foam
{
namespace PBEMethods
{


defineTypeNameAndDebug(MOM, 0);
addToRunTimeSelectionTable(PBEMethod, MOM, dictionary);
// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

MOM::MOM
(
    const dictionary& pbeProperties,
    const phaseModel& phase
)
:
    PBEMethod(pbeProperties, phase),
    MOMDict_(pbeProperties.subDict("MOMCoeffs")),
    dispersedPhase_(phase),
    moments_{{
    volScalarField
    (
        IOobject
        (
            "m0",
            dispersedPhase_.U().mesh().time().timeName(),
            dispersedPhase_.U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        dispersedPhase_.U().mesh()
    ),
    volScalarField
    (
        IOobject
        (
            "m1",
            dispersedPhase_.U().mesh().time().timeName(),
            dispersedPhase_.U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        dispersedPhase_.U().mesh()
    ),
    volScalarField
    (
        IOobject
        (
            "m2",
            dispersedPhase_.U().mesh().time().timeName(),
            dispersedPhase_.U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        dispersedPhase_.U().mesh()
    ),
    volScalarField
    (
        IOobject
        (
            "m3",
            dispersedPhase_.U().mesh().time().timeName(),
            dispersedPhase_.U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        dispersedPhase_.U().mesh()
    )
    }},
    scaleD_
    (
        dimensionedScalar
        (
            "scaleDiameter",
            dimensionSet(0,0,0,0,0),
            MOMDict_.lookup("scaleDiameter")
        )
    ),
    d32_
    (
        IOobject
        (
            "diameter",
            dispersedPhase_.U().mesh().time().timeName(),
            dispersedPhase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        moments_[3] * scaleD_ * scaleM3_  / moments_[2]
    ),
    gamma_alpha_
    (
        IOobject
        (
            "gamma_alpha",
            dispersedPhase_.U().mesh().time().timeName(),
            dispersedPhase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dispersedPhase_.U().mesh(),
        dimensionedScalar("gamma_alpha", dimLength, 0.0)
    ),
    gamma_beta_
    (
        IOobject
        (
            "gamma_beta",
            dispersedPhase_.U().mesh().time().timeName(),
            dispersedPhase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dispersedPhase_.U().mesh(),
        dimensionedScalar("gamma_beta", dimless, 0.0)
    ),
    gamma_c0_
    (
        IOobject
        (
            "gamma_c0",
            dispersedPhase_.U().mesh().time().timeName(),
            dispersedPhase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dispersedPhase_.U().mesh(),
        dimensionedScalar("c0", dimless, 0.0)
    ),
    scaleM3_(MOMDict_.lookupOrDefault<scalar>("scaleM3",1.0)),
    Nf_(MOMDict_.lookupOrDefault<scalar>("daughterDropletsNr",2.0)),
    maxD_
    (
        dimensionedScalar
        (
            "maxDiameter",
            dimLength,
            MOMDict_.lookup("maxDiameter")
        )
    ),
    minD_
    (
        dimensionedScalar
        (
            "minDiameter",
            dimLength, 
            MOMDict_.lookup("minDiameter")
        )
    ),
    minGammaAlpha_
    (
        "minGammaAlpha",
        dimLength,
        MOMDict_.lookupOrDefault<scalar>("minGammaAlpha", 1e-05)
    ),
    maxGammaBeta_(MOMDict_.lookupOrDefault<scalar>("maxGammaBeta", 10.0)),
    maxDiameterMultiplicator_
    (
        MOMDict_.lookupOrDefault<scalar>("maxDiameterMultiplicator", 10.0)
    ),
    integrationSteps_(MOMDict_.lookupOrDefault<scalar>("integrationSteps", 10)),
    bList_(integrationSteps_),
    cList_(pow(integrationSteps_, 2)),
    gammaList_(integrationSteps_)
{
}

MOM::~MOM()
{
}
void MOM::correct()
{
    Info<< "updating size moments" << endl;
    gamma_c0_ = moments_[0];
    gamma_alpha_ = ( moments_[0] * moments_[2] - pow(moments_[1], 2) )
        / (moments_[0] * moments_[1]
           + dimensionedScalar("small", dimensionSet(0,1,0,0,0),SMALL ) );
    gamma_beta_ = pow(moments_[1] , 2)
        / (moments_[2] * moments_[0] - pow(moments_[1],2)
           + dimensionedScalar("small", dimensionSet(0,2,0,0,0),SMALL ) );

    Info<< "gamma parameters:" << endl;
    printAvgMaxMin(gamma_c0_);
    printAvgMaxMin(gamma_alpha_);
    printAvgMaxMin(gamma_beta_);

    //TODO: get rid of 4
    std::array<volScalarField, 4> mSources_{{
        volScalarField
        (
            IOobject
            (
                "m0Source",
                dispersedPhase_.U().mesh().time().timeName(),
                dispersedPhase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            momentSourceTerm(0) / scaleD_.value()
        ),
        volScalarField
        (
            IOobject
            (
                "m1Source",
                dispersedPhase_.U().mesh().time().timeName(),
                dispersedPhase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            momentSourceTerm(1) / pow(scaleD_.value(), 2)
        ),
        volScalarField
        (
            IOobject
            (
                "m2Source",
                dispersedPhase_.U().mesh().time().timeName(),
                dispersedPhase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            momentSourceTerm(2) / pow(scaleD_.value(), 3)
        ),
        volScalarField
        (
            IOobject
            (
                "m3Source",
                dispersedPhase_.U().mesh().time().timeName(),
                dispersedPhase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dispersedPhase_.U().mesh(),
            dimensionedScalar("m3Source", pow(dimLength, 3) / dimTime, 0.0)
        )
    }};

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
        forAll(dispersedPhase_.U().mesh().C(), celli)
        {
            auto& mSource_i = mSource[celli];
            if ( std::isnan(mSource_i) )
            {
                mSource_i = 0.0;
            }
        }

        printAvgMaxMin(mSource);
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
        printAvgMaxMin(moment);
    }

    d32_ = moments_[3] * scaleD_ * scaleM3_ / moments_[2];
    d32_ = min(d32_, maxD_);
    d32_ = max(d32_,  minD_);

    printAvgMaxMin(d32_);
}

const volScalarField MOM::d() const
{
    return d32_;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> MOM::momentSourceTerm(label momenti)
{
    volScalarField bS = breakupSourceTerm(momenti);
    volScalarField cS = coalescenceSourceTerm(momenti) / scaleD_.value();
    Info << "breakup (moment "
        << momenti << ") " << bS.weightedAverage(dispersedPhase_.U().mesh().V()).value() << endl;
    Info << "coalescence (moment "
        << momenti << ") " << cS.weightedAverage(dispersedPhase_.U().mesh().V()).value() << endl;
    return cS + bS;// breakupSourceTerm(momenti);
}


tmp<volScalarField> MOM::coalescenceSourceTerm(label momenti)
{
    dimensionedScalar xiDim = dimensionedScalar("xiDim", dimLength, 1.0);

    scalar d32Mean = max( d32_.weightedAverage(dispersedPhase_.U().mesh().V()).value() , SMALL);
    scalar maxArg = maxDiameterMultiplicator_ * d32Mean ;
    label maxIter = integrationSteps_;
    dimensionedScalar dx = xiDim * maxArg / integrationSteps_;

    dimensionedScalar arg_i1 = xiDim;
    dimensionedScalar arg_j1 = xiDim;
    label iterator = 0;
    label iterator2 = 0;

    //value of the integral
    volScalarField sum
        (
            IOobject
            (
                "Sbr",
                dispersedPhase_.U().mesh().time().timeName(),
                dispersedPhase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            dispersedPhase_.U().mesh(),
            dimensionedScalar
            (
                "Sbr", 
                pow(dimLength, momenti - 1) / dimTime,
                0
            ) 
        );
    volScalarField sum_i0(sum);


    volScalarField sum2
        (
            IOobject
            (
                "Sbr",
                dispersedPhase_.U().mesh().time().timeName(),
                dispersedPhase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            dispersedPhase_.U().mesh(),
            dimensionedScalar
            (
                "Sbr", 
                pow(dimLength, momenti) / dimTime,
                0
            ) 
        );

    return tmp<volScalarField>( new volScalarField(sum2));

}

void MOM::printAvgMaxMin(const volScalarField &v) const
{
    Info<< v.name() << ": avg, max,min "
        << v.weightedAverage(dispersedPhase_.U().mesh().V()).value()
        << ", " << max(v).value()
        << ", " << min(v).value() << endl;
}

tmp<volScalarField> MOM::breakupSourceTerm(label momenti)
{
    //this function integrates the momentum kernel according to an assumed shape
    //of the distribution (gamma distribution with method of moments estimator
    //describing the shape with three independent parameters related to first
    //three moments (m0,m1,m2)
    
    //numerical integration is performed so integration limits must be assumed
    //integral is calculated on interval 0 < d < maxDiameterMultiplicator * (mean d32)
    //TODO: is there a beter way to identify integral limits?
    //integration step is calculated according to prescribed number of
    //integration steps

    scalar d32Mean = max(d32_.weightedAverage(dispersedPhase_.U().mesh().V()).value(), SMALL);
    scalar maxArg = maxDiameterMultiplicator_ * d32Mean ;
    label maxIter = integrationSteps_;
    scalar dx = maxArg / integrationSteps_;

    scalar arg_i1 = 0.0;
    label iterator = 0;

    //value of the integral
    volScalarField sum
        (
            IOobject
            (
                "Sbr",
                dispersedPhase_.U().mesh().time().timeName(),
                dispersedPhase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            dispersedPhase_.U().mesh(),
            dimensionedScalar
            (
                "Sbr", 
                pow(dimLength, momenti) /dimTime,
                0
            ) 
        );

    return tmp<volScalarField>( new volScalarField(sum));
}

tmp<volScalarField> MOM::gamma(const volScalarField& x)
{
    volScalarField result
        (
            IOobject
            (
                "gamma",
                dispersedPhase_.U().mesh().time().timeName(),
                dispersedPhase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dispersedPhase_.U().mesh(),
            dimensionedScalar("gamma", dimless, 0.0)
        );
    using constant::mathematical::pi;

    forAll(result, celli)
    {
        if (x[celli] < 0.5)
        {
            result[celli] = 
                pi / (sin(pi * x[celli]) * internal::gamma(1 - x[celli]));
        }
        else
        {
            result[celli] = internal::gamma(x[celli]);
        }
    }
    return tmp<volScalarField>( new volScalarField(result));
}


tmp<volScalarField> MOM::gammaDistribution(const dimensionedScalar xi)
{
	gamma_c0_ = max(SMALL, gamma_c0_);
	gamma_alpha_ = max(minGammaAlpha_, gamma_alpha_);
	gamma_beta_ = max(SMALL, gamma_beta_);
	gamma_beta_ = min(maxGammaBeta_, gamma_beta_);

	return tmp<volScalarField>
        (
            gamma_c0_ * pow(xi, gamma_beta_ - 1) 
            * pow(constant::mathematical::e, - xi / gamma_alpha_ )
            / 
            (
                pow(gamma_alpha_, gamma_beta_) 
                * (gamma(gamma_beta_) + SMALL) 
            )
        );
}


}//end namespace PBEMethods
}//end namespace Foam

// ************************************************************************* //
