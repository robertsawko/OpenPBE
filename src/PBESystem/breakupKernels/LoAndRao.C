/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "LoAndRao.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace breakupKernels
{
defineTypeNameAndDebug(LoAndRao, 0);
addToRunTimeSelectionTable
(
    breakupKernel,
    LoAndRao,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LoAndRao::LoAndRao
(
    const dictionary& breakupDict,
    const phaseModel& continuousPhase,
    const phaseModel& dispersedPhase
)
:
   breakupKernel(breakupDict,continuousPhase,dispersedPhase)//,
   //TODO this lookup does not work; there is no object of type mEuler::RASModel
   //in the registry
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
//
LoAndRao::~LoAndRao()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> LoAndRao::S(const volScalarField& xi) const
{
    const multiphaseTurbulence::turbulenceModel& turbulence
    (
       continuousPhase_.mesh().lookupObject
       <multiphaseTurbulence::turbulenceModel>("turbulenceModel")
    );
    //critical capillary number; TODO read or default from dictionary
    scalar omegaCrit = 0.5;

    volScalarField epsilon = turbulence.epsilon();
    volScalarField k = turbulence.k();
    volScalarField nut = turbulence.nut();

    //for max and min functions
    dimensionedScalar epsilonDim("epsilonDim", epsilon.dimensions(), 1.0);
    dimensionedScalar kDim("kDim", k.dimensions(), 1.0);
    dimensionedScalar nuDim("nuDim", dimensionSet(0, 2, -1, 0, 0, 0, 0), 1.0);

    const dimensionedScalar nuc = continuousPhase_.nu();
    const dimensionedScalar nud = dispersedPhase_.nu();
    const dimensionedScalar rhoc = continuousPhase_.rho();
    const dimensionedScalar rhod = dispersedPhase_.rho();
    //surface tension TODO: read from dictionary
    scalar sigma = 0.047;

    volScalarField kolmogorovL(
        pow(
            pow(nuc,3.0) / max(epsilon, SMALL*epsilonDim),
            0.25
        )
    );

    volScalarField shearRate
    (
        pow
        (
            max
            (
                epsilon, 
                SMALL * epsilonDim
            )
            / max
            (
                nuc, SMALL * nuDim
            )
            , 0.5
        )
    );

    volScalarField inertialCriticalD(
        "inertialCriticalD",
        //correction factor; unused at this time
        //(1 + 4.6 * max(PBE_.getAlpha(),dimensionedScalar("zero",dimless,0)) ) *
        pow(
            2.0 * sigma * 0.23 / rhoc
            ,0.6
        ) *
        pow(
            max(epsilon, SMALL*epsilonDim),
            -0.4
        )
    );

    //Info << "inertia critical diameter ,min, max = "
    //    << inertialCriticalD.weightedAverage(k.mesh().V()).value() << ", "
    //    << min (inertialCriticalD).value() << ", "
    //    << max (inertialCriticalD).value() << endl;


    volScalarField viscousCriticalD(
        2.0 * omegaCrit * sigma 
        / (max(nuc,SMALL*nuDim)
           * rhoc * shearRate)
    );

    //Info << "viscous critical diameter ,min, max = "
    //    << viscousCriticalD.weightedAverage(k.mesh().V()).value() << ", "
    //    << min (viscousCriticalD).value() << ", "
    //    << max (viscousCriticalD).value() << endl;

    volScalarField viscousBreakupRate(
        IOobject
        (
            "viscousBreakupRate",
            k.mesh().time().timeName(),
            k.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        k.mesh(),
        max(nuc, SMALL * nuDim) * rhoc / sigma
    );

    volScalarField breakupKernel
    (
        IOobject
        (
            "breakupKernel",
            k.mesh().time().timeName(),
            k.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        k.mesh(),
        dimensionedScalar
        (
            "breakupKernel", 
            xi.dimensions() / dimTime,
            0.0
        )
    );

    scalar d_p = 0;
    forAll(k.mesh().C(), celli){
        d_p = xi[celli];

        if(d_p < kolmogorovL[celli])
        {
            if(d_p > viscousCriticalD[celli])
            {
                breakupKernel[celli] += viscousBreakupRate[celli] * d_p;
            }
        }
        if(d_p > inertialCriticalD[celli]){
            breakupKernel[celli] += 2.0 * constant::mathematical::pi * 0.2
                * pow
                (
                    (3.0 * rhod.value() + 2.0 * rhoc.value()) * pow(d_p , 3) 
                    / (192.0 * sigma)
                    ,0.5
                );
        }
    }

    return tmp<volScalarField>
        (
            new volScalarField
            (
                /*IOobject
                  (
                  "S",
                  xi.mesh().time().timeName(),
                  xi.mesh(),
                  IOobject::NO_READ,
                  IOobject::NO_WRITE,
                  false
                  ),
                  xi.mesh(),
                //TODO Is the dimension correct?
                dimensionedScalar("S", xi.dimensions()/dimTime, 0) */
                breakupKernel
            )
        );
}

} //End namespace breakupKernels
} //End namespace Foam
// ************************************************************************* //
