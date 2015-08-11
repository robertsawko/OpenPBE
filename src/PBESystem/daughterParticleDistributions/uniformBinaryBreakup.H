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

Class
    Foam::uniformBinaryBreakup

Description

SourceFiles
    uniformBinaryBreakupI.H
    uniformBinaryBreakup.C
    uniformBinaryBreakupIO.C

\*---------------------------------------------------------------------------*/

#ifndef uniformBinaryBreakup_H
#define uniformBinaryBreakup_H

#include "daughterParticleDistribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace daughterParticleDistributions
{

/*---------------------------------------------------------------------------*\
                         Class uniformBinaryBreakup Declaration
\*---------------------------------------------------------------------------*/

class uniformBinaryBreakup
:
    public daughterParticleDistribution
{
 public:

  // Static data members

  //- Static data staticData
  TypeName("uniformBinaryBreakup")

  // Constructors

  uniformBinaryBreakup(
      const dictionary& pbeDict
  );

  //- Destructor
  virtual ~uniformBinaryBreakup() = default;

  // Member Functions

  dimensionedScalar beta(
      const dimensionedScalar& xi1,
      const dimensionedScalar& xi2) const override;

  dimensionedScalar moment(
            const dimensionedScalar& xi, label k) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
