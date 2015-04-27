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
    Foam::Madsen

Description

SourceFiles
    MadsenI.H
    Madsen.C
    MadsenIO.C

\*---------------------------------------------------------------------------*/

#ifndef Madsen_H
#define Madsen_H

#include "coalescenceKernel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceKernels
{

/*---------------------------------------------------------------------------*\
                         Class Madsen Declaration
\*---------------------------------------------------------------------------*/

class Madsen
:
    public coalescenceKernel
{
 private:
  // Private data
  const scalar sigma_;

  //- Description of data_

  // Private Member Functions

  //- Disallow default bitwise copy construct
  Madsen(const Madsen&);

  //- Disallow default bitwise assignment
  void operator=(const Madsen&);

 public:

  // Static data members

  //- Static data staticData
  TypeName("Madsen");

  // Constructors

  Madsen
  (
      const dictionary& coalescenceDict,
      const phaseModel& continuousPhase,
      const phaseModel& dispersedPhase
  );

  //- Destructor
  virtual ~Madsen();

  // Member Functions

  tmp<volScalarField> S(
      const volScalarField& xi1,
      const volScalarField& xi2) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
