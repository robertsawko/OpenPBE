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

#include "breakupKernel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
defineTypeNameAndDebug(breakupKernel, 0);
defineRunTimeSelectionTable(breakupKernel, dictionary);
}
// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::breakupKernel::breakupKernel
(
    const dictionary& breakupDict,
    const phaseModel& dispersedPhase
) 
:
    breakupDict_(breakupDict),
    continuousPhase_(dispersedPhase.otherPhase()),
    dispersedPhase_(dispersedPhase)
{
};

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
Foam::autoPtr<Foam::breakupKernel> Foam::breakupKernel::New
(
    const dictionary& pbeDict,
    const phaseModel& dispersedPhase
)
{
    word breakupKernelType(pbeDict.lookup("breakupKernel"));
    Info << "Selecting breakup model: " << breakupKernelType << endl;
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(breakupKernelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("breakupKernel::New")
            << "Unknown breakupKernelType type "
            << breakupKernelType << endl << endl
            << "Valid breakupKernel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return cstrIter()(
        pbeDict.subDict(breakupKernelType+"Coeffs"), dispersedPhase);
};


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::breakupKernel::~breakupKernel()
{
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
