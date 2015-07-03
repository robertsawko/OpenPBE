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
    Foam::constantCoalescence

Description

SourceFiles
    constantCoalescenceI.H
    constantCoalescence.C
    constantCoalescenceIO.C

\*---------------------------------------------------------------------------*/

#ifndef constantCoalescence_H
#define constantCoalescence_H

#include "coalescenceKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceKernels
{

/*---------------------------------------------------------------------------*\
                         Class constantCoalescence Declaration
\*---------------------------------------------------------------------------*/

class constantCoalescence
:
    public coalescenceKernel
{
 private:
  scalar constant_;
  // Private data

  //- Description of data_

  // Private Member Functions

  //- Disallow default bitwise copy construct
  constantCoalescence(const constantCoalescence&);

  //- Disallow default bitwise assignment
  void operator=(const constantCoalescence&);

 public:

  // Static data members

  //- Static data staticData
  TypeName("constantCoalescence");

  // Constructors

  constantCoalescence
  (
      const dictionary& coalescenceDict,
      const phaseModel& dispersedPhase
  );

  //- Destructor
  virtual ~constantCoalescence();

  // Member Functions

  tmp<volScalarField> S(
      const volScalarField& xi1,
      const volScalarField& xi2) const;

  const dimensionedScalar S(
      const dimensionedScalar& xi1,
      const dimensionedScalar& xi2) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
