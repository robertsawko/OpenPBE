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

#include "uniformBinaryBreakup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace daughterParticleDistributions
{
defineTypeNameAndDebug(uniformBinaryBreakup, 0);
addToRunTimeSelectionTable
(
    daughterParticleDistribution,
    uniformBinaryBreakup,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uniformBinaryBreakup::uniformBinaryBreakup
(
    const dictionary& dpdDict
)
:
   daughterParticleDistribution(dpdDict)
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar uniformBinaryBreakup::beta(scalar xi1,
    scalar xi2
) const
{
    return (2.0 / xi2)  * pos(xi2 - xi1);
}

dimensionedScalar uniformBinaryBreakup::moment
(
        const dimensionedScalar &xi,
        label k
) const
{
    return 2*pow(xi,k)/(k+1);
}
} //End namespace daughterParticleDistributions
} //End namespace Foam
// ************************************************************************* //
