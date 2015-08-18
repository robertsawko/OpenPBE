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

#include "daughterParticleDistribution.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
defineTypeNameAndDebug(daughterParticleDistribution, 0);
defineRunTimeSelectionTable(daughterParticleDistribution, dictionary)
}
// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::daughterParticleDistribution::daughterParticleDistribution
(
    const dictionary& dpdDict
) 
:
    dpdDict_(dpdDict)
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
Foam::autoPtr<Foam::daughterParticleDistribution>
Foam::daughterParticleDistribution::New
(
    const dictionary& pbeDict
)
{
    word daughterParticleDistributionType(
        pbeDict.lookup("daughterParticleDistribution")
    );
    Info<< "Selecting daughter particle model: " <<
        daughterParticleDistributionType << endl;
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(daughterParticleDistributionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("daughterParticleDistribution::New")
            << "Unknown daughterParticleDistributionType type "
            << daughterParticleDistributionType << endl << endl
            << "Valid daughterParticleDistribution types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return cstrIter()
    (
        pbeDict.subDict(daughterParticleDistributionType + "Coeffs")
    );
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
