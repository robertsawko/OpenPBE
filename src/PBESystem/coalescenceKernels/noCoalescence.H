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
    Foam::noCoalescence

Description

SourceFiles
    noCoalescenceI.H
    noCoalescence.C
    noCoalescenceIO.C

\*---------------------------------------------------------------------------*/

#ifndef noCoalescence_H
#define noCoalescence_H

#include "coalescenceKernel.H"

namespace Foam
{
namespace coalescenceKernels
{

class noCoalescenceImpl
{
public:
    scalar S() const
    {
        return 0.;
    }
};

class noCoalescence
:
    public coalescenceKernel
{
 public:
  TypeName("noCoalescence")

  noCoalescence
  (
      const dictionary& coalescenceDict,
      const phaseModel& dispersedPhase
  );

  virtual ~noCoalescence() = default;

  scalar S(
      const dimensionedScalar& xi1,
      const dimensionedScalar& xi2,
      label celli) const override
  {
      return impl_.S();
  }

 private:
  noCoalescenceImpl impl_;
};


} // End namespace coalescenceModels
} // End namespace Foam

#endif
