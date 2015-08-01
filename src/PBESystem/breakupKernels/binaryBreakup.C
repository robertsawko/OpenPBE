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

#include "binaryBreakup.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace breakupKernels
{
defineTypeNameAndDebug(binaryBreakup, 0);
addToRunTimeSelectionTable
(
    breakupKernel,
    binaryBreakup,
    dictionary
);

binaryBreakup::binaryBreakup
(
    const dictionary& breakupDict,
    const phaseModel& dispersedPhase
)
:
    breakupKernel(breakupDict, dispersedPhase)
{
}

dimensionedScalar binaryBreakupImpl::S(const dimensionedScalar &xi) const {
    return dimensionedScalar("S", dimless / dimTime, pow(xi.value(), 2));
}

} //End namespace breakupKernels
} //End namespace Foam
