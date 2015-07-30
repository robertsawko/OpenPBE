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

#include "coalescenceKernel.H"

namespace Foam
{
defineTypeNameAndDebug(coalescenceKernel, 0);
defineRunTimeSelectionTable(coalescenceKernel, dictionary)
}

Foam::coalescenceKernel::coalescenceKernel
(
    const dictionary& coalescenceDict,
    const phaseModel& dispersedPhase
) 
:
    coalescenceDict_(coalescenceDict),
    continuousPhase_(dispersedPhase.otherPhase()),
    dispersedPhase_(dispersedPhase)
{
}

Foam::autoPtr<Foam::coalescenceKernel> Foam::coalescenceKernel::New
(
    const dictionary& pbeDict,
    const phaseModel& dispersedPhase
)
{
    word coalescenceKernelType(pbeDict.lookup("coalescenceKernel"));
    Info << "Selecting coalescence model: " << coalescenceKernelType << endl;
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(coalescenceKernelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("coalescenceKernel::New")
            << "Unknown coalescenceKernelType type "
            << coalescenceKernelType << endl << endl
            << "Valid coalescenceKernel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return cstrIter()
    (
        pbeDict.subDict(coalescenceKernelType + "Coeffs"),
        dispersedPhase
    );
}
