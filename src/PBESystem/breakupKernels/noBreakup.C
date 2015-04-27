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

#include "noBreakup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace breakupKernels
{
defineTypeNameAndDebug(noBreakup, 0);
addToRunTimeSelectionTable
(
    breakupKernel,
    noBreakup,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noBreakup::noBreakup
(
    const dictionary& breakupDict,
    const phaseModel& dispersedPhase
)
:
   breakupKernel(breakupDict, dispersedPhase)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
//
noBreakup::~noBreakup()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> noBreakup::S(const volScalarField& xi) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
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
            dimensionedScalar("S", xi.dimensions()/dimTime, 0) 
        )
    );
}

dimensionedScalar noBreakup::S(const dimensionedScalar& xi) const
{
    return dimensionedScalar("S", xi.dimensions()/dimTime, 0);
}

} //End namespace breakupKernels
} //End namespace Foam
// ************************************************************************* //
