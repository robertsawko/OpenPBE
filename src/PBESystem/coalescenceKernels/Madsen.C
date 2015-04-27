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

#include "Madsen.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace coalescenceKernels
{
defineTypeNameAndDebug(Madsen, 0);
addToRunTimeSelectionTable
(
    coalescenceKernel,
    Madsen,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Madsen::Madsen
(
    const dictionary& coalescenceDict,
    const phaseModel& continuousPhase,
    const phaseModel& dispersedPhase
)
:
    coalescenceKernel(coalescenceDict,continuousPhase,dispersedPhase),
    sigma_( coalescenceDict.lookupOrDefault<scalar>("sigma",0.047) )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
//
Madsen::~Madsen()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> Madsen::S
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

    volScalarField dRel
        (
            0.5 * (xi1 + xi2)
        );
    //TODO: if we want to take polycelerity into account we have to have
    //relative velocity here:
    //volScalarField Urel = phase1.U() - phase2.U()

    volScalarField Urel = 
        sqrt
        (
            2.0 * 3.6 * pow
            (
                epsilonSafe * dRel,
                2.0 / 3.0
            )
        );

    volScalarField WeCrit = 
        rhod * pow(Urel,2) * min(xi1, xi2) 
        / sigma;

    volScalarField gamma = max(xi1, xi2) / ( min(xi1,xi2) + SMALL*xiDim );
    volScalarField fGamma = pow(gamma,3) - 2.4 * pow(gamma, 2) + 2.7 * gamma;

    volScalarField oneField = volScalarField 
    (
        IOobject
        (
            "oneField",
            k.mesh().time().timeName(),
            k.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        k.mesh(),
        dimensionedScalar("oneField", fGamma.dimensions() / WeCrit.dimensions(), 1.0)
    );

    volScalarField EBounce = 
        min
        (
            oneField, 
            pow
            (
                WeCrit / (4.8 * fGamma + SMALL),
                1.0 / 3.0
            )
        );

    volScalarField ECoalescence = 
        min
        (
            oneField, 
            4.8 * fGamma / ( WeCrit + SMALL)
        );

    //return min(EBounce, ECoalescence) * constant::mathematical::pi
        //* pow(dRel, 2) * Urel * pow(xiDim,-2);

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
            min(EBounce, ECoalescence) * constant::mathematical::pi
            * pow(dRel, 2) * Urel * pow(xiDim,-2)
        )
    );
}

} //End namespace coalescenceKernels
} //End namespace Foam
// ************************************************************************* //
