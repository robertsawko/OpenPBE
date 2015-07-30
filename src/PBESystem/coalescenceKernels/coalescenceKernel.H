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

namespace Foam
{

class coalescenceKernel
{
 protected:
  dictionary coalescenceDict_;
  const phaseModel& continuousPhase_;
  const phaseModel& dispersedPhase_;

 public:
  coalescenceKernel(const coalescenceKernel&) = delete;
  coalescenceKernel& operator=(const coalescenceKernel&) = delete;

  TypeName("coalescenceKernel")

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
  )

  coalescenceKernel
  (
      const dictionary& coalescenceDict,
      const phaseModel& dispersedPhase
  );

  static autoPtr<coalescenceKernel> New
  (
      const dictionary& coalescenceDict,
      const phaseModel& dispersedPhase
  );

  virtual ~coalescenceKernel() = default;

  virtual const dimensionedScalar S(
      const dimensionedScalar& xi1,const dimensionedScalar& xi2) const = 0;
};


} // End namespace Foam


#endif
