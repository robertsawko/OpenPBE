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
    Foam::cTimesD

Description

SourceFiles
    cTimesDI.H
    cTimesD.C
    cTimesDIO.C

\*---------------------------------------------------------------------------*/

#ifndef cTimesD_H
#define cTimesD_H

#include "breakupKernel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace breakupKernels
{

/*---------------------------------------------------------------------------*\
                         Class cTimesD Declaration
\*---------------------------------------------------------------------------*/

class cTimesD
:
    public breakupKernel
{
 private:
  // Private data
  const dimensionedScalar c_;
  //const multiphaseTurbulence::turbulenceModel& turbulence_;
//  const volScalarField& turbulence;

  //- Description of data_

  // Private Member Functions

  //- Disallow default bitwise copy construct
  cTimesD(const cTimesD&);

  //- Disallow default bitwise assignment
  void operator=(const cTimesD&);

 public:

  // Static data members

  //- Static data staticData
  TypeName("cTimesD");

  // Constructors

  cTimesD
  (
      const dictionary& breakupDict,
      const phaseModel& continuousPhase,
      const phaseModel& dispersedPhase
  );

  //- Destructor
  virtual ~cTimesD();

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
