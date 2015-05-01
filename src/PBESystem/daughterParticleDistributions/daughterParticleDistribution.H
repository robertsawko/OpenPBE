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

//#include "fvCFD.H"
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
 private:
  // Private Member Functions

  //- Disallow default bitwise copy construct
  daughterParticleDistribution(const daughterParticleDistribution&);

  //- Disallow default bitwise assignment
  void operator=(const daughterParticleDistribution&);

  // Protected data
 protected:
  dictionary dpdDict_;

 public:
  TypeName("daughterParticleDistribution");

  declareRunTimeSelectionTable
  (
      autoPtr,
      daughterParticleDistribution,
      dictionary,
      (
          dictionary dpdDict
      ),
      (dpdDict)
  );

  // Static data members

  //- Static data staticData


  // Constructors

  //- Construct null
  //daughterParticleDistribution();

  //- Construct from components
  daughterParticleDistribution(const dictionary& dpdDict);


  // Selectors
  static autoPtr<daughterParticleDistribution> New(
      const dictionary& pbeDict
  );

  //- Destructor
  virtual ~daughterParticleDistribution();

  // Member Functions

  //- Source function for the coalescence kernel used in PBE
  virtual tmp<volScalarField> beta(
      const volScalarField& xi1,const volScalarField& xi2) const = 0;

  virtual const dimensionedScalar beta(
      const dimensionedScalar& xi1,const dimensionedScalar& xi2) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
