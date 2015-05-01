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

#include "PBEMethod.H"

namespace Foam
{
    defineTypeNameAndDebug(PBEMethod, 0);
    defineRunTimeSelectionTable(PBEMethod, dictionary);
}
// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PBEMethod::PBEMethod(
    const dictionary& pbeProperties,
    const phaseModel& phase
) : 
    coalescence_(coalescenceKernel::New(pbeProperties, phase)),
    breakup_(breakupKernel::New(pbeProperties, phase)),
    daughterParticleDistribution_(
        daughterParticleDistribution::New(pbeProperties)
    ) 
{
    //TODO add some error checking. Do phases exist etc.
}

/*
namespace Foam
{
const phaseModel& PBEMethod::continuousPhase() const
{
    return *phases_[continuousPhaseName_];
}
const phaseModel& PBEMethod::dispersedPhase() const
{
    return *phases_[dispersedPhaseName_];
}
phaseModel& PBEMethod::continuousPhase()
{
    return *phases_[continuousPhaseName_];
}
phaseModel& PBEMethod::dispersedPhase()
{
    return *phases_[dispersedPhaseName_];
}

*/
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//bool Foam::PBEMethod::read()
//{
//}

Foam::autoPtr<Foam::PBEMethod> Foam::PBEMethod::New(
    const dictionary& pbeProperties,
    const phaseModel& phase
)
{
    word pbeMethod
    (
        pbeProperties.lookup("method")
    );

    Info << "Selecting population balance model: "
        << pbeMethod << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(pbeMethod);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("PBEMethod::New")
            << "Unknown pbeModelType type "
            << pbeMethod << endl << endl
            << "Valid PBEMethod types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(pbeProperties, phase);
}

// ************************************************************************* //
