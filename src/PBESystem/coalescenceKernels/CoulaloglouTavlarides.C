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

#include "CoulaloglouTavlarides.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace coalescenceKernels
{
defineTypeNameAndDebug(CoulaloglouTavlaridesC, 0);
addToRunTimeSelectionTable
(
    coalescenceKernel,
    CoulaloglouTavlaridesC,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

CoulaloglouTavlaridesC::CoulaloglouTavlaridesC
(
    const dictionary& coalescenceDict,
    const phaseModel& continuousPhase,
    const phaseModel& dispersedPhase
)
:
    coalescenceKernel(coalescenceDict,continuousPhase,dispersedPhase),
    c1_( coalescenceDict.lookupOrDefault<scalar>("c1",0.00487) ),
    c2_( coalescenceDict.lookupOrDefault<scalar>("c2",0.008) ),
    gamma_( coalescenceDict.lookupOrDefault<scalar>("gamma",0.0) ),
    sigma_( coalescenceDict.lookupOrDefault<scalar>("sigma",0.047) )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
//
CoulaloglouTavlaridesC::~CoulaloglouTavlaridesC()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> CoulaloglouTavlaridesC::S
(
    const volScalarField& xi1,
    const volScalarField& xi2
) const
{


    const multiphaseTurbulence::turbulenceModel& turbulence
    (
       continuousPhase_.mesh().lookupObject
       <multiphaseTurbulence::turbulenceModel>("turbulenceModel")
    );

    volScalarField epsilon = turbulence.epsilon();
    volScalarField k = turbulence.k();
    volScalarField nut = turbulence.nut();

    const dimensionedScalar nud = dispersedPhase_.nu();
    const dimensionedScalar nuc = continuousPhase_.nu();
    const dimensionedScalar rhod = dispersedPhase_.rho();
    const dimensionedScalar rhoc = continuousPhase_.rho();
    dimensionedScalar epsilonDim("epsilonDim", epsilon.dimensions(), 1.0);
    dimensionedScalar kDim("kDim", k.dimensions(), 1.0);
    dimensionedScalar nuDim("nuDim", nud.dimensions(), 1.0);
    dimensionedScalar xiDim = dimensionedScalar("unitLength", dimLength, 1.0);

    volScalarField epsilonSafe = max(epsilon, SMALL * epsilonDim);

    dimensionedScalar sigma("sigma", dimensionSet(1,0,-2,0,0,0,0), sigma_);

    volScalarField frequency = 
        c1_ * pow(epsilonSafe, 1.0 / 3.0) * pow(xi1 + xi2, 2) * 
        pow
        ( 
            pow(xi1, 2.0 / 3.0) + pow(xi2, 2.0 / 3.0),
            0.5
        ) 
        / (1 + gamma_);


    volScalarField rate = 
        exp
        (
            - c2_ * rhod * epsilonSafe * (rhod * nud) * pow(xiDim,-2)
            / ( pow(sigma,2) * pow( 1 + gamma_,3) )
            * pow
                (
                    xi1 * xi2 / (xi1 + xi2),
                    4.0
                )
        );


    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "S",
                xi1.mesh().time().timeName(),
                xi1.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            frequency * rate * pow(xiDim, -2)
        )
    );
}

} //End namespace coalescenceKernels
} //End namespace Foam
// ************************************************************************* //
