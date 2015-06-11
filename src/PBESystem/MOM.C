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
    m0_
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
    m1_
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
    m2_
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
    m3_
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
    ),
    m0Source_
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
    m1Source_
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
    m2Source_
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
    m3Source_
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
    ),
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



//void MOM::solveAlphas()
//{
//    //Solve continuos phase
//    PtrList<surfaceScalarField> phiAlphaCorrs(phases_.size());
//    int phasei = 0;

//    forAllIter(PtrDictionary<phaseModel>, phases_, iter)
//    {
//        phaseModel& phase1 = iter();
//        volScalarField& alpha1 = phase1;

//        phase1.phiAlpha() =
//            dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0);

//        phiAlphaCorrs.set
//        (
//            phasei,
//            new surfaceScalarField
//            (
//                fvc::flux
//                (
//                    phi_,
//                    phase1,
//                    "div(phi," + alpha1.name() + ')'
//                )
//            )
//        );

//        surfaceScalarField& phiAlphaCorr = phiAlphaCorrs[phasei];

//        forAllIter(PtrDictionary<phaseModel>, phases_, iter2)
//        {
//            phaseModel& phase2 = iter2();
//            volScalarField& alpha2 = phase2;

//            if (&phase2 == &phase1) continue;

//            surfaceScalarField phir(phase1.phi() - phase2.phi());

//            //RS: Do we need the cAlpha at all in PBE code?
//            scalarCoeffSymmTable::const_iterator cAlpha
//                (
//                    cAlphas_.find(interfacePair(phase1, phase2))
//                );

//            if (cAlpha != cAlphas_.end())
//            {
//                surfaceScalarField phic
//                    (
//                        (mag(phi_) + mag(phase1.phi() - phase2.phi()))
//                        / phase_.U().mesh().magSf()
//                    );

//                phir += min(cAlpha() * phic, max(phic)) * nHatf(phase1, phase2);
//            }

//            word phirScheme
//            (
//                "div(phir," + alpha2.name() + ',' + alpha1.name() + ')'
//            );

//            phiAlphaCorr += fvc::flux
//            (
//                -fvc::flux(-phir, phase2, phirScheme),
//                phase1,
//                phirScheme
//            );
//        }

//        MULES::limit
//        (
//            geometricOneField(),
//            phase1,
//            phi_,
//            phiAlphaCorr,
//            zeroField(),
//            zeroField(),
//            1,
//            0,
//            3,
//            true
//        );

//        phasei++;
//    }

//    MULES::limitSum(phiAlphaCorrs);

//    volScalarField sumAlpha
//    (
//        IOobject
//        (
//            "sumAlpha",
//            phase_.U().mesh().time().timeName(),
//            phase_.U().mesh()
//        ),
//        phase_.U().mesh(),
//        dimensionedScalar("sumAlpha", dimless, 0)
//    );

//    phasei = 0;

//    forAllIter(PtrDictionary<phaseModel>, phases_, iter)
//    {
//        phaseModel& phase1 = iter();

//        surfaceScalarField& phiAlpha = phiAlphaCorrs[phasei];
//        phiAlpha += upwind<scalar>(phase_.U().mesh(), phi_).flux(phase1);

//        MULES::explicitSolve
//            (
//                geometricOneField(),
//                phase1,
//                phiAlpha,
//                zeroField(),
//                zeroField()
//            );
//        /*

//           this provided correct prediction od diameter ( d32 = m3 / m2) but failed to
//           preserve boundedness of alpha
//           fvScalarMatrix alphaEqn
//            (
//                fvm::ddt(phase1)
//                + fvm::div(phiAlpha,phase1)
//            );

//        alphaEqn.relax();
//        alphaEqn.solve();*/

//        if (&phase1 == &dispersedPhase()){
//            //reuse this field in solution method for moments
//            limitedFlux_ = phiAlpha;
//        }

//        phase1.phiAlpha() += phiAlpha;

//        Info<< phase1.name() << " volume fraction, min, max = "
//            << phase1.weightedAverage(phase_.U().mesh().V()).value()
//            << ' ' << min(phase1).value()
//            << ' ' << max(phase1).value()
//            << endl;

//        sumAlpha += phase1;

//        phasei++;
//    }

//    Info<< "Phase-sum volume fraction, min, max = "
//        << sumAlpha.weightedAverage(phase_.U().mesh().V()).value()
//        << ' ' << min(sumAlpha).value()
//        << ' ' << max(sumAlpha).value()
//        << endl;

//    calcAlphas();
//}


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
    Info << "updating size moments" << endl;
    gamma_c0_ = m0_;
    gamma_alpha_ = ( m0_ * m2_ - pow(m1_, 2) ) 
        / (m0_ * m1_
           + dimensionedScalar("small", dimensionSet(0,1,0,0,0),SMALL ) );
    gamma_beta_ = pow(m1_ , 2) 
        / (m2_ * m0_ - pow(m1_,2)
           + dimensionedScalar("small", dimensionSet(0,2,0,0,0),SMALL ) );

    d32_ = m3_ * scaleD_ * scaleM3_  / m2_;
    d32_ = min(d32_,  maxD_);
    d32_ = max(d32_,  minD_);
    Info << "d32: avg, max,min "
        << d32_.weightedAverage(phase_.U().mesh().V()).value()
        << ", " << max(d32_).value()
        << ", " << min(d32_).value() << endl;

    Info << "gamma parameters:" << endl;
    Info << "C0: avg, max,min "
        << gamma_c0_.weightedAverage(phase_.U().mesh().V()).value()
        << ", " << max(gamma_c0_).value()
        << ", " << min(gamma_c0_).value() << endl;
    Info << "alpha: avg, max,min " 
        << gamma_alpha_.weightedAverage(phase_.U().mesh().V()).value()
        << ", " << max(gamma_alpha_).value()
        << ", " << min(gamma_alpha_).value() << endl;
    Info << "beta: avg, max,min " 
        << gamma_beta_.weightedAverage(phase_.U().mesh().V()).value()
        << ", " << max(gamma_beta_).value()
        << ", " << min(gamma_beta_).value() << endl;

    m0Source_ = momentSourceTerm(0) / scaleD_.value();
    m1Source_ = momentSourceTerm(1) / pow(scaleD_.value(), 2);
    m2Source_ = momentSourceTerm(2) / pow(scaleD_.value(), 3);

    //TODO: is there a better way to assure that source terms on boundaries 
    //are equal to 0?
    //maybe using internalField() instead of the whole field somewhere?
    m0Source_.boundaryField() = 0;
    m1Source_.boundaryField() = 0;
    m2Source_.boundaryField() = 0;
    m3Source_.boundaryField() = 0;

    //fix for nan values
    forAll(phase_.U().mesh().C(), celli)
    {
        if ( m0Source_[celli] != m0Source_[celli] )
        {
            m0Source_[celli] = 0.0;
        }
        if ( m1Source_[celli] != m1Source_[celli] )
        {
            m1Source_[celli] = 0.0;
        }
        if ( m2Source_[celli] != m2Source_[celli] )
        {
            m2Source_[celli] = 0.0;
        }
    }

    Info << "moment sources:" << endl;
    Info << "m0: avg, max,min " 
        << m0Source_.weightedAverage(phase_.U().mesh().V()).value()
        << ", " << max(m0Source_).value()
        << ", " << min(m0Source_).value() << endl;
    Info << "m1: avg, max,min " 
        << m1Source_.weightedAverage(phase_.U().mesh().V()).value()
        << ", " << max(m1Source_).value()
        << ", " << min(m1Source_).value() << endl;
    Info << "m2: avg, max,min " 
        << m2Source_.weightedAverage(phase_.U().mesh().V()).value()
        << ", " << max(m2Source_).value()
        << ", " << min(m2Source_).value() << endl;

    fvScalarMatrix m0Eqn
        (
            fvm::ddt(m0_)
            + fvm::div(phase_.phi(),m0_)
            //+ fvm::div(limitedFlux_,m0_)
            ==
            m0Source_
        ); 
    m0Eqn.relax();
    m0Eqn.solve();

    fvScalarMatrix m1Eqn
        (
            fvm::ddt(m1_)
            + fvm::div(phase_.phi(),m1_)
            //+ fvm::div(limitedFlux_,m1_)
            ==
            m1Source_
        ); 
    m1Eqn.relax();
    m1Eqn.solve();

    fvScalarMatrix m2Eqn
        (
            fvm::ddt(m2_)
            + fvm::div(phase_.phi(),m2_)
            //+ fvm::div(limitedFlux_,m2_)
            ==
            m2Source_
        ); 
    m2Eqn.relax();
    m2Eqn.solve();

    fvScalarMatrix m3Eqn
        (
            fvm::ddt(m3_)
            + fvm::div(phase_.phi(),m3_)
            //+ fvm::div(limitedFlux_,m3_)
            ==
            m3Source_
        ); 

    m3Eqn.relax();
    m3Eqn.solve();

    m0_ = max(m0_, dimensionedScalar("m0", m0_.dimensions(), SMALL) );
    m1_ = max(m1_, dimensionedScalar("m1", m1_.dimensions(), SMALL) );
    m2_ = max(m2_, dimensionedScalar("m2", m2_.dimensions(), SMALL) );
    m3_ = max(m3_, dimensionedScalar("m3", m3_.dimensions(), SMALL) );
    
    Info << "m0: avg, max,min " 
        << m0_.weightedAverage(m0_.mesh().V()).value() 
        << ", " << max(m0_).value()
        << ", " << min(m0_).value() << endl;

    Info << "m1: avg, max,min " 
        << m1_.weightedAverage(m1_.mesh().V()).value() 
        << ", " << max(m1_).value()
        << ", " << min(m1_).value() << endl;

    Info << "m2: avg, max,min " 
        << m2_.weightedAverage(m2_.mesh().V()).value() 
        << ", " << max(m2_).value()
        << ", " << min(m2_).value() << endl;

    Info << "m3: avg, max,min " 
        << m3_.weightedAverage(m3_.mesh().V()).value() 
        << ", " << max(m3_).value()
        << ", " << min(m3_).value() << endl;


    d32_ = m3_ * scaleD_ * scaleM3_ / m2_;
    //d32_ = dispersedPhase() * scaleD_ *6.0 / (m2_ * constant::mathematical::pi);
    d32_ = min(d32_, maxD_);
    d32_ = max(d32_,  minD_);

    expectedD_ = gamma_alpha_ * gamma_beta_ * scaleD_.value();

    Info << "diameter: avg, max,min " 
        << d32_.weightedAverage(phase_.U().mesh().V()).value()
        << ", " << max(d32_).value()
        << ", " << min(d32_).value() << endl;
    Info << "expected diameter: avg, max,min " 
        << expectedD_.weightedAverage(phase_.U().mesh().V()).value()
        << ", " << max(expectedD_).value()
        << ", " << min(expectedD_).value() << endl;
};

void MOM::solve()
{
//    forAllIter(PtrDictionary<phaseModel>, phases_, iter)
//    {
//        iter().correct();
//    }

//    const Time& runTime = phase_.U().mesh().time();

//    if (nAlphaSubCycles_ > 1)
//    {
//        dimensionedScalar totalDeltaT = runTime.deltaT();

//        /*volScalarField m0_0 = m0_.oldTime();
//        volScalarField m1_0 = m1_.oldTime();
//        volScalarField m2_0 = m2_.oldTime();
//        volScalarField m3_0 = m3_.oldTime();*/

//        PtrList<volScalarField> alpha0s(phases_.size());
//        PtrList<surfaceScalarField> phiSums(phases_.size());

//        int phasei = 0;
//        forAllIter(PtrDictionary<phaseModel>, phases_, iter)
//        {
//            phaseModel& phase = iter();
//            volScalarField& alpha = phase;

//            alpha0s.set
//            (
//                phasei,
//                new volScalarField(alpha.oldTime())
//            );

//            phiSums.set
//            (
//                phasei,
//                new surfaceScalarField
//                (
//                    IOobject
//                    (
//                        "phiSum" + alpha.name(),
//                        runTime.timeName(),
//                        phase_.U().mesh()
//                    ),
//                    phase_.U().mesh(),
//                    dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
//                )
//            );
//            phasei++;
//        }

//        for
//        (
//            subCycleTime alphaSubCycle
//            (
//                const_cast<Time&>(runTime),
//                nAlphaSubCycles_
//            );
//            !(++alphaSubCycle).end();
//        )
//        {
//            //Is this in the right order? phiSums and solveAlphas
//            solveAlphas();

//            int phasei = 0;
//            forAllIter(PtrDictionary<phaseModel>, phases_, iter)
//            {
//                phiSums[phasei] +=
//                    (runTime.deltaT() / totalDeltaT) * iter().phi();
//                phasei++;
//            }
//        }

//        phasei = 0;
//        forAllIter(PtrDictionary<phaseModel>, phases_, iter)
//        {
//            phaseModel& phase = iter();
//            volScalarField& alpha = phase;

//            //Is this supposed to correct velocity field? I don't think this is
//            //doing anything.
//            phase.phi() = phiSums[phasei];

//            // Correct the time index of the field
//            // to correspond to the global time
//            alpha.timeIndex() = runTime.timeIndex();

//            /*m0_.timeIndex() = runTime.timeIndex();
//            m1_.timeIndex() = runTime.timeIndex();
//            m2_.timeIndex() = runTime.timeIndex();
//            m3_.timeIndex() = runTime.timeIndex();*/

//            // Reset the old-time field value
//            alpha.oldTime() = alpha0s[phasei];
//            alpha.oldTime().timeIndex() = runTime.timeIndex();

//            /*m0_.oldTime() = m0_0;
//            m1_.oldTime() = m1_0;
//            m2_.oldTime() = m2_0;
//            m3_.oldTime() = m3_0;
//            m0_.oldTime().timeIndex() = runTime.timeIndex();
//            m1_.oldTime().timeIndex() = runTime.timeIndex();
//            m2_.oldTime().timeIndex() = runTime.timeIndex();
//            m3_.oldTime().timeIndex() = runTime.timeIndex();*/

//            phasei++;
//        }

//    }
//    else
//    {
//        solveAlphas();
//        //updateMoments();
//    }

    updateMoments();
}



}//end namespace PBEMethods
}//end namespace Foam

// ************************************************************************* //
