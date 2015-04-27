/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::PBEMethod

Description
    Incompressible multi-phase mixture with built in solution for the
    phase fractions with interface compression for interface-capturing.

    Derived from transportModel so that it can be unsed in conjunction with
    the incompressible turbulence models.

    Surface tension and contact-angle is handled for the interface
    between each phase-pair.

SourceFiles
    PBEMethod.C

\*---------------------------------------------------------------------------*/

#ifndef PBEMethod_H
#define PBEMethod_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H"
#include "breakupKernels/breakupKernel.H"
#include "coalescenceKernels/coalescenceKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class PBEMethod Declaration
\*---------------------------------------------------------------------------*/

class PBEMethod
{

protected:

  // Protected data
  //- coalescence kernel
  autoPtr<coalescenceKernel> coalescence_;
  //- breakup kernel
  autoPtr<breakupKernel> breakup_;

public:

  //- Runtime type information
  TypeName("PBEMethod");


  // Declare runtime construction
  declareRunTimeSelectionTable
  (
     autoPtr,
     PBEMethod,
     dictionary,
     (
         const dictionary& pbeProperties,
         const phaseModel& phase
     ),
     (pbeProperties, phase)
  );

  // Constructors

  //- Construct from components
  PBEMethod(
      const dictionary& pbeProperties,
      const phaseModel& phase
  );

  static autoPtr<PBEMethod> New
  (
      const dictionary& pbeProperties,
      const phaseModel& phase
  );

  //- Destructor
  virtual ~PBEMethod()
  {}

  // Member Functions

  //- Dummy correct
  virtual void correct() = 0;

  ////- Return the Sauter-mean diameter
  virtual const volScalarField d() const = 0;

  ////- Read base transportProperties dictionary
  //virtual bool read(const dictionary& phaseProperties)
  //{
      //return true;
  //};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
