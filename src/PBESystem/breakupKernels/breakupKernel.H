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
    Foam::breakupKernel

Description

SourceFiles
    breakupKernelI.H
    breakupKernel.C
    breakupKernelIO.C

\*---------------------------------------------------------------------------*/

#ifndef breakupKernel_H
#define breakupKernel_H

//#include "fvCFD.H"
#include "phaseModel.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
  Class breakupKernel Declaration
\*---------------------------------------------------------------------------*/

class breakupKernel
{
 private:
  // Private Member Functions

  //- Disallow default bitwise copy construct
  breakupKernel(const breakupKernel&);

  //- Disallow default bitwise assignment
  void operator=(const breakupKernel&);

  // Protected data
 protected:
  dictionary breakupDict_;
  const phaseModel* continuousPhase_;
  const phaseModel* dispersedPhase_;

 public:
  TypeName("breakupKernel");

  declareRunTimeSelectionTable
  (
      autoPtr,
      breakupKernel,
      dictionary,
      (
          dictionary breakupDict, 
          const phaseModel& dispersedPhase
      ),
      (breakupDict, dispersedPhase)
  );

  // Static data members

  //- Static data staticData


  // Constructors

  //- Construct null
  breakupKernel();

  //- Construct from components
  breakupKernel
  (
      const dictionary& breakupDict,
      const phaseModel& dispersedPhase
  );


  // Selectors
  static autoPtr<breakupKernel> New
  (
      const dictionary& breakupDict,
      const phaseModel& dispersedPhase
  );

  //- Destructor
  virtual ~breakupKernel();

  // Member Functions

  //- Source function for the breakup kernel used in PBE
  virtual tmp<volScalarField> S(const volScalarField& xi) const = 0;

  //- Source function for the breakup kernel used in PBE
  virtual dimensionedScalar S(const dimensionedScalar& xi) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
