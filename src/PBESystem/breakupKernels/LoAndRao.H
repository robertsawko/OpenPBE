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
    Foam::LoAndRao

Description

SourceFiles
    LoAndRaoI.H
    LoAndRao.C
    LoAndRaoIO.C

\*---------------------------------------------------------------------------*/

#ifndef LoAndRao_H
#define LoAndRao_H

#include "breakupKernel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace breakupKernels
{

/*---------------------------------------------------------------------------*\
                         Class LoAndRao Declaration
\*---------------------------------------------------------------------------*/

class LoAndRao
:
    public breakupKernel
{
 private:
  // Private data
  //const multiphaseTurbulence::turbulenceModel& turbulence_;
//  const volScalarField& turbulence;

  //- Description of data_

  // Private Member Functions

  //- Disallow default bitwise copy construct
  LoAndRao(const LoAndRao&);

  //- Disallow default bitwise assignment
  void operator=(const LoAndRao&);

 public:

  // Static data members

  //- Static data staticData
  TypeName("LoAndRao");

  // Constructors

  LoAndRao
  (
      const dictionary& breakupDict,
      const phaseModel& continuousPhase,
      const phaseModel& dispersedPhase
  );

  //- Destructor
  virtual ~LoAndRao();

  // Member Functions

  tmp<volScalarField> S(const volScalarField& xi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakupModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
