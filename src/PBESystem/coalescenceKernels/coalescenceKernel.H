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
    Foam::coalescenceKernel

Description

SourceFiles
    coalescenceKernelI.H
    coalescenceKernel.C
    coalescenceKernelIO.C

\*---------------------------------------------------------------------------*/

#ifndef coalescenceKernel_H
#define coalescenceKernel_H

//#include "fvCFD.H"
#include "phaseModel.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
  Class coalescenceKernel Declaration
\*---------------------------------------------------------------------------*/

class coalescenceKernel
{
 private:
  // Private Member Functions

  //- Disallow default bitwise copy construct
  coalescenceKernel(const coalescenceKernel&);

  //- Disallow default bitwise assignment
  void operator=(const coalescenceKernel&);

  // Protected data
 protected:
  dictionary coalescenceDict_;
  const phaseModel& continuousPhase_;
  const phaseModel& dispersedPhase_;

 public:
  TypeName("coalescenceKernel");

  declareRunTimeSelectionTable
  (
      autoPtr,
      coalescenceKernel,
      dictionary,
      (
          dictionary coalescenceDict, 
          const phaseModel& dispersedPhase
      ),
      (coalescenceDict, dispersedPhase)
  );

  // Static data members

  //- Static data staticData


  // Constructors

  //- Construct null
  //coalescenceKernel();

  //- Construct from components
  coalescenceKernel
  (
      const dictionary& coalescenceDict,
      const phaseModel& dispersedPhase
  );


  // Selectors
  static autoPtr<coalescenceKernel> New
  (
      const dictionary& coalescenceDict,
      const phaseModel& dispersedPhase
  );

  //- Destructor
  virtual ~coalescenceKernel();

  // Member Functions

  //- Source function for the coalescence kernel used in PBE
  virtual tmp<volScalarField> S(
      const volScalarField& xi1,const volScalarField& xi2) const = 0;

  virtual const dimensionedScalar S(
      const dimensionedScalar& xi1,const dimensionedScalar& xi2) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
