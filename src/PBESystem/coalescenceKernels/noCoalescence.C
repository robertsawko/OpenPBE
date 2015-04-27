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

#include "noCoalescence.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace coalescenceKernels
{
defineTypeNameAndDebug(noCoalescence, 0);
addToRunTimeSelectionTable
(
    coalescenceKernel,
    noCoalescence,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noCoalescence::noCoalescence
(
    const dictionary& coalescenceDict,
    const phaseModel& dispersedPhase
)
:
   coalescenceKernel(coalescenceDict, dispersedPhase)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
//
noCoalescence::~noCoalescence()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> noCoalescence::S
(
    const volScalarField& xi1,
    const volScalarField& xi2
) const
{
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
            xi1.mesh(),
            //TODO Is the dimension correct?
            dimensionedScalar("S", xi1.dimensions()/dimTime, 0) 
        )
    );
}

const dimensionedScalar noCoalescence::S
(
    const dimensionedScalar& xi1,
    const dimensionedScalar& xi2
) const
{
    return dimensionedScalar("S", xi1.dimensions()/dimTime, 0);
}
} //End namespace coalescenceKernels
} //End namespace Foam
// ************************************************************************* //
