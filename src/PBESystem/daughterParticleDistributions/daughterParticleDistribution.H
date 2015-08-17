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
    Foam::daughterParticleDistribution

Description

SourceFiles
    daughterParticleDistributionI.H
    daughterParticleDistribution.C
    daughterParticleDistributionIO.C

\*---------------------------------------------------------------------------*/

#ifndef daughterParticleDistribution_H
#define daughterParticleDistribution_H

#include "phaseModel.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
  Class daughterParticleDistribution Declaration
\*---------------------------------------------------------------------------*/

class daughterParticleDistribution
{
 protected:
  dictionary dpdDict_;

 public:
  daughterParticleDistribution(const daughterParticleDistribution&) = delete;
  daughterParticleDistribution& operator=(const daughterParticleDistribution&)
                                          = delete;

  TypeName("daughterParticleDistribution")

  declareRunTimeSelectionTable
  (
      autoPtr,
      daughterParticleDistribution,
      dictionary,
      (
          dictionary dpdDict
      ),
      (dpdDict)
  )

  daughterParticleDistribution(const dictionary& dpdDict);

  static autoPtr<daughterParticleDistribution> New(
      const dictionary& pbeDict
  );

  virtual ~daughterParticleDistribution() = default;

  virtual dimensionedScalar beta(
      const dimensionedScalar& xi1,const dimensionedScalar& xi2) const = 0;

  virtual dimensionedScalar moment(
          const dimensionedScalar& xi, label k) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
