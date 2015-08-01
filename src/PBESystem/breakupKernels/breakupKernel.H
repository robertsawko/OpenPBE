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

#include "phaseModel.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

namespace Foam
{

class breakupKernel
{
 protected:
  dictionary breakupDict_;
  const phaseModel& continuousPhase_;
  const phaseModel& dispersedPhase_;

 public:
  breakupKernel(const breakupKernel&) = delete;
  breakupKernel& operator=(const breakupKernel&) = delete;

  TypeName("breakupKernel")

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
  )

  breakupKernel
  (
      const dictionary& breakupDict,
      const phaseModel& dispersedPhase
  );

  static autoPtr<breakupKernel> New
  (
      const dictionary& breakupDict,
      const phaseModel& dispersedPhase
  );

  virtual ~breakupKernel() = default;

  virtual dimensionedScalar S(const dimensionedScalar &xi,
                              label celli) const = 0;
};


} // End namespace Foam

#endif
