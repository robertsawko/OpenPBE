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

#include "dummy.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

namespace Foam
{
namespace PBEMethods
{

defineTypeNameAndDebug(dummy, 0);
addToRunTimeSelectionTable(PBEMethod, dummy, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dummy::dummy
(
    const dictionary& pbeProperties,
    const phaseModel& phase
)
:
    PBEMethod(pbeProperties, phase),
    d_
    (
        IOobject
        (
            IOobject::groupName("d", phase.name()),
            phase.U().time().timeName(),
            phase.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase.U().mesh(),
        dimensionedScalar(
            "d",
            dimLength,
            readScalar(
                pbeProperties.subDict("dummyCoeffs").lookup("constantDiameter")
            )
        )
    )
{
}

} // End namespace PBEMethods
} // End namespace Foam
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
