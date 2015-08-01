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

#include "constantCoalescence.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace coalescenceKernels
{
defineTypeNameAndDebug(constantCoalescence, 0);
addToRunTimeSelectionTable
(
    coalescenceKernel,
    constantCoalescence,
    dictionary
);

constantCoalescence::constantCoalescence
(
    const dictionary& coalescenceDict,
    const phaseModel& dispersedPhase
)
:
    coalescenceKernel(coalescenceDict, dispersedPhase),
    impl_(readScalar(coalescenceDict.lookup("constant")))
{
}

constantCoalescenceImpl::constantCoalescenceImpl(scalar constant):
    constant_(constant)
{

}

const dimensionedScalar constantCoalescenceImpl::S
(
        const dimensionedScalar& xi1,
        const dimensionedScalar& xi2
) const
{
    return dimensionedScalar("S", dimless/dimTime, constant_);
}

} //End namespace coalescenceKernels
} //End namespace Foam
