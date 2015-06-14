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
    phase_(phase),
    moments_{{
    volScalarField
    (
        IOobject
        (
            "m0",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.U().mesh()
    ),
    volScalarField
    (
        IOobject
        (
            "m1",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.U().mesh()
    ),
    volScalarField
    (
        IOobject
        (
            "m2",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.U().mesh()
    ),
    volScalarField
    (
        IOobject
        (
            "m3",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.U().mesh()
    )
    }},
    mSources_{{
    volScalarField
    (
        IOobject
        (
            "m0Source",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.U().mesh(),
        dimensionedScalar("m0Source", dimless / dimTime, 0.0)
    ),
    volScalarField
    (
        IOobject
        (
            "m1Source",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.U().mesh(),
        dimensionedScalar("m1Source", dimLength / dimTime, 0.0)
    ),
    volScalarField
    (
        IOobject
        (
            "m2Source",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.U().mesh(),
        dimensionedScalar("m2Source", pow(dimLength, 2) / dimTime, 0.0)
    ),
    volScalarField
    (
        IOobject
        (
            "m3Source",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase_.U().mesh(),
        dimensionedScalar("m3Source", dimless / dimTime, 0.0)
    )
    }},
    d32_
    (
        //argument set to true applies this diameter in other places in the
        //solver (i.e. in drag models) leaving it empty or false prevent of
        //usign it in different places ('debug' mode)
        phase_.d()
        /*IOobject
        (
            "diameter",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.U().mesh(),
        dimensionedScalar("diameter", dimless, 0.0)*/
    ),
    expectedD_
    (
        IOobject
        (
            "expectedDiameter",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.U().mesh(),
        dimensionedScalar("expectedDiameter", dimLength, 0.0)
    ),
    gamma_alpha_
    (
        IOobject
        (
            "gamma_alpha",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.U().mesh(),
        dimensionedScalar("gamma_alpha", dimLength, 0.0)
    ),
    gamma_beta_
    (
        IOobject
        (
            "gamma_beta",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.U().mesh(),
        dimensionedScalar("gamma_beta", dimless, 0.0)
    ),
    gamma_c0_
    (
        IOobject
        (
            "gamma_c0",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.U().mesh(),
        dimensionedScalar("c0", dimless, 0.0)
    ),
    scaleD_
    (
        dimensionedScalar
        (
            "scaleDiameter",
            dimensionSet(0,3,0,0,0),
            MOMDict_.lookup("scaleDiameter")
        )
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
    gammaList_(integrationSteps_)/*,
    limitedFlux_(phase_.phi()),
    nAlphaSubCycles_
    (
        readLabel(phase_.U().mesh().solutionDict().subDict("PIMPLE").lookup("nAlphaSubCycles"))
    )*/
{
}

void MOM::correct()
{
    updateMoments();
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
        << momenti << ") " << bS.weightedAverage(phase_.U().mesh().V()).value() << endl;
    Info << "coalescence (moment "
        << momenti << ") " << cS.weightedAverage(phase_.U().mesh().V()).value() << endl;
    return cS + bS;// breakupSourceTerm(momenti);
}


tmp<volScalarField> MOM::coalescenceSourceTerm(label momenti)
{
    dimensionedScalar xiDim = dimensionedScalar("xiDim", dimLength, 1.0);

    scalar d32Mean = max( d32_.weightedAverage(phase_.U().mesh().V()).value() , SMALL);
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
                phase_.U().mesh().time().timeName(),
                phase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase_.U().mesh(),
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
                phase_.U().mesh().time().timeName(),
                phase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase_.U().mesh(),
            dimensionedScalar
            (
                "Sbr", 
                pow(dimLength, momenti) / dimTime,
                0
            ) 
        );

    //value of the function being integrated (for previous argument)
    volScalarField f_i0
        (
            IOobject
            (
                "f_i0",
                phase_.U().mesh().time().timeName(),
                phase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase_.U().mesh(),
            dimensionedScalar
            (
                "f_i0", 
                pow(dimLength, momenti - 2) / dimTime,
                0
            ) 
        );
    //value of the function being integrated (for current argument)
    volScalarField f_i1
        (
            IOobject
            (
                "f_i1",
                phase_.U().mesh().time().timeName(),
                phase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase_.U().mesh(),
            dimensionedScalar
            (
                "f_i1", 
                pow(dimLength, momenti - 2) / dimTime,
                0
            ) 
        );
    
    //in this method argument of integral is a single scalar value
    //this field extrapolates scalar argument to volScalarField
    //so it can be passed to 'S' method of breakup model
    volScalarField argField
        (
            IOobject
            (
                "arg",
                phase_.U().mesh().time().timeName(),
                phase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase_.U().mesh(),
            dimensionedScalar
            (
                "arg", 
                dimLength,
                0
            ) 
        );
    volScalarField argField2(argField);


    if (momenti == 0){
        //TODO: see next one, few lines below
        Info << "number of integration steps for coalescence module: " 
            << maxIter <<  " x " << maxIter << endl;
    }

    scalar avgCoalescenceRate = 0.0;

    for (label iter = 1; iter <= maxIter; iter++)
    {
        arg_i1 = iter * dx;
        argField = arg_i1;
        sum = dimensionedScalar("zero",pow(dimLength, momenti - 1) / dimTime , 0.0);

        for (label iter2 = 1; iter2 <= maxIter; iter2++)
        {
            arg_j1 = iter2 * dx;
            argField2 = arg_j1;

            if (momenti == 0){
                //TODO: can it be done better?
                //can we buffer breakup and gamma values insted of calculating them
                //few time without messing with the class structure?
                //
                //it could have been list of list but one long list is slightly
                //more convinient to use in this case
                cList_.set
                    (
                        (iter - 1)* maxIter + (iter2 - 1),
                        coalescence_-> S(argField, argField2)
                    );
                avgCoalescenceRate += 
                    cList_[(iter-1)*maxIter + (iter2 - 1)].weightedAverage
                    (
                        phase_.U().mesh().V()
                    ).value()
                    / pow(maxIter,2);
                if(iter2 == maxIter && iter==maxIter){

                    Info << "Average coalescence rate: " << avgCoalescenceRate << endl;

                }
                //in sourceTerm method we have:
                //return coalescenceSourceTerm(momenti) + breakupSourceTerm(momenti);
                //so breakup is calculated pravious to coalescence and there
                //is no need to update gammaList here
                /*gammaList_.set
                (
                    iter2 - 1,
                    gammaDistribution(arg_j1  / scaleD_.value()) 
                );*/

            }

            f_i1 =
                (
                    pow
                    (
                        pow(arg_i1, 3) + pow(arg_j1, 3),
                        momenti / 3.0
                    )
                    - pow(arg_i1, momenti)
                    - pow(arg_j1, momenti)
                )
                * cList_[(iter-1)*maxIter + (iter2 - 1)]//coalescence_ -> S(argField, argField2)  
                * gammaList_[iter-1] //gammaDistribution(arg_i1)
                * gammaList_[iter2-1] * xiDim;//gammaDistribution(arg_j1);
            sum += (f_i1 + f_i0) * dx * 0.5;
            f_i0 = f_i1;
        }

        sum2 = ( sum_i0 + sum ) * dx * 0.5;
        sum_i0 = sum;
    }

    return tmp<volScalarField>( new volScalarField(sum2));

}

void MOM::printAvgMaxMin(const volScalarField &v) const
{
    Info<< v.name() << ": avg, max,min "
        << v.weightedAverage(phase_.U().mesh().V()).value()
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

    scalar d32Mean = max(d32_.weightedAverage(phase_.U().mesh().V()).value(), SMALL);
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
                phase_.U().mesh().time().timeName(),
                phase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase_.U().mesh(),
            dimensionedScalar
            (
                "Sbr", 
                pow(dimLength, momenti) /dimTime,
                0
            ) 
        );

    //value of the function being integrated (for previous argument)
    volScalarField f_i0
    (
        IOobject
        (
            "f_i0",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        phase_.U().mesh(),
        dimensionedScalar
        (
            "f_i0", 
            pow(dimLength, momenti - 1) /dimTime,
            0
        ) 
    );
    //value of the function being integrated (for current argument)
    volScalarField f_i1
    (
        IOobject
        (
            "f_i1",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        phase_.U().mesh(),
        dimensionedScalar
        (
            "f_i1", 
            pow(dimLength, momenti - 1) /dimTime,
            0
        ) 
    );
    
    //in this method argument of integral is a single scalar value
    //this field extrapolates scalar argument to volScalarField
    //so it can be passed to 'S' method of breakup model
    volScalarField argField
    (
        IOobject
        (
            "arg",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        phase_.U().mesh(),
        dimensionedScalar
        (
            "arg", 
            dimLength,
            0
        ) 
    );

    dimensionedScalar xiDim = dimensionedScalar("xiDim", dimLength, 1.0);

    scalar avgBreakupRate = 0.0;

    scalar maxGamma = 0;
    if (momenti==0){
        //TODO: see next one, few lines below
        Info << "number of integration steps for breakup module: " 
            << maxIter << endl;
    }
    for (label iter = 1; iter <= maxIter; iter++){
        arg_i1 = iter * dx;
        argField = arg_i1 * xiDim;

        if (momenti == 0){
            //TODO: can it be done better?
            //can we buffer breakup and gamma values insted of calculating them
            //few time without messing with the class structure?
            bList_.set
            (
                iter - 1,
                breakup_-> S(argField)
            );

            avgBreakupRate += 
                bList_[iter-1].weightedAverage(phase_.U().mesh().V()).value()
                / maxIter;

            gammaList_.set
            (
                iter - 1,
                gammaDistribution(arg_i1 / scaleD_.value()) 
            );

            maxGamma = max(maxGamma, gammaList_[iter-1].weightedAverage(phase_.U().mesh().V()).value() );

            if (iter == maxIter)
            {
               Info << "Gamma distribution values for maximum considered diameter: " 
                    << gammaList_[iter -1].weightedAverage(phase_.U().mesh().V()).value()
                    << ", max: " << max(gammaList_[iter-1].internalField())
                    << ", min: " << min(gammaList_[iter-1]).value() << endl;
                Info << "Maximum value for the distribution is: " << maxGamma << endl;

                Info << "Average breakup rate: " << avgBreakupRate << endl;
            }
            /*Info << "breakup kernels:" 
                << bList[iter -1].weightedAverage(phase_.U().mesh().V()).value()
                << ", max: " << max(bList[iter-1].internalField())
                << ", min: " << min(bList[iter-1]).value() << endl;*/

        }

        f_i1 =
            pow(arg_i1*xiDim, momenti)
            * ( pow(Nf_, ((3.0 - momenti) / 3.0)) - 1)
            * bList_[iter-1] //breakup_-> S(argField)
            * gammaList_[iter-1]; //gammaDistribution(arg_i1);
        /*Info << "multiplier: " << 
            pow(arg_i1*xiDim, momenti)
            << endl;
        Info << "multiplier2: " << 
            ( pow(Nf_, ((3.0 - momenti) / 3.0)) - 1)
            << endl;*/
        sum += (f_i1 + f_i0) * dx * 0.5 * xiDim;
            /*Info << "f0: " 
                << f_i0.weightedAverage(phase_.U().mesh().V()).value()
                << ", max: " << max(f_i0.internalField())
                << ", min: " << min(f_i0.internalField()) << endl;
            Info << "f1: " 
                << f_i1.weightedAverage(phase_.U().mesh().V()).value()
                << ", max: " << max(f_i1.internalField())
                << ", min: " << min(f_i1.internalField()) << endl;
            Info << "sum: " 
                << sum.weightedAverage(phase_.U().mesh().V()).value()
                << ", max: " << max(sum.internalField())
                << ", min: " << min(sum.internalField()) << endl;*/
        //Info << max(sum.internalField()) << endl;
        f_i0 = f_i1;
    }


    return tmp<volScalarField>( new volScalarField(sum));

}

tmp<volScalarField> MOM::gamma(const volScalarField& x)
{
    volScalarField result
        (
            IOobject
            (
                "gamma",
                phase_.U().mesh().time().timeName(),
                phase_.U().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phase_.U().mesh(),
            dimensionedScalar("gamma", dimless, 0.0)
        );
    scalar pi = constant::mathematical::pi;

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

void MOM::updateMoments()
{
    Info<< "updating size moments" << endl;
    gamma_c0_ = moments_[0];
    gamma_alpha_ = ( moments_[0] * moments_[2] - pow(moments_[1], 2) )
        / (moments_[0] * moments_[1]
           + dimensionedScalar("small", dimensionSet(0,1,0,0,0),SMALL ) );
    gamma_beta_ = pow(moments_[1] , 2)
        / (moments_[2] * moments_[0] - pow(moments_[1],2)
           + dimensionedScalar("small", dimensionSet(0,2,0,0,0),SMALL ) );

    d32_ = moments_[3] * scaleD_ * scaleM3_  / moments_[2];
    d32_ = min(d32_,  maxD_);
    d32_ = max(d32_,  minD_);

    printAvgMaxMin(d32_);
    Info<< "gamma parameters:" << endl;
    printAvgMaxMin(gamma_c0_);
    printAvgMaxMin(gamma_alpha_);
    printAvgMaxMin(gamma_beta_);

    mSources_[0] = momentSourceTerm(0) / scaleD_.value();
    mSources_[1] = momentSourceTerm(1) / pow(scaleD_.value(), 2);
    mSources_[2] = momentSourceTerm(2) / pow(scaleD_.value(), 3);
    //TODO: Why is the last source term not scaled here?

    //TODO: is there a better way to assure that source terms on boundaries 
    //are equal to 0?
    //maybe using internalField() instead of the whole field somewhere?
    //TODO: loops can likely be merged
    for (auto& mSource : mSources_)
        mSource.boundaryField() = 0;

    //fix for nan values
    forAll(phase_.U().mesh().C(), celli)
    {
        for (auto& mSource : mSources_)
        {
            auto& mSource_i = mSource[celli];
            if ( std::isnan(mSource_i) )
            {
                mSource_i = 0.0;
            }
        }
    }

    Info<< "moment sources:" << endl;
    for (auto& mSource : mSources_)
        printAvgMaxMin(mSource);

    for (std::size_t i = 0; i<moments_.size(); ++i)
    {
        fvScalarMatrix mEqn
            (
                fvm::ddt(moments_[i])
                + fvm::div(phase_.phi(),moments_[i])
                //+ fvm::div(limitedFlux_,moments_[i])
                ==
                mSources_[i]
            );
        mEqn.relax();
        mEqn.solve();
    }

    for (auto& moment : moments_)
    {
        moment = max(moment, dimensionedScalar(moment.name(),
                                               moment.dimensions(), SMALL));
        printAvgMaxMin(moment);
    }

    d32_ = moments_[3] * scaleD_ * scaleM3_ / moments_[2];
    //d32_ = dispersedPhase() * scaleD_ *6.0 / (moments_[2] * constant::mathematical::pi);
    d32_ = min(d32_, maxD_);
    d32_ = max(d32_,  minD_);

    expectedD_ = gamma_alpha_ * gamma_beta_ * scaleD_.value();

    printAvgMaxMin(d32_);
    printAvgMaxMin(expectedD_);
}


}//end namespace PBEMethods
}//end namespace Foam

// ************************************************************************* //
