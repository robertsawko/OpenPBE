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

namespace Foam
{
namespace coalescenceKernels
{

class constantCoalescenceImpl
{
public:
    constantCoalescenceImpl(scalar constant);

    const dimensionedScalar S(const dimensionedScalar &xi1,
                              const dimensionedScalar &xi2) const;

private:
    scalar constant_;

};

class constantCoalescence
:
    public coalescenceKernel
{


 public:
  TypeName("constantCoalescence")

  constantCoalescence
  (
      const dictionary& coalescenceDict,
      const phaseModel& dispersedPhase
  );

  virtual ~constantCoalescence() = default;

  dimensionedScalar S(
      const dimensionedScalar& xi1,
      const dimensionedScalar& xi2,
      label celli) const override
  {
      return impl_.S(xi1, xi2);
  }

 private:
  constantCoalescenceImpl impl_;
};


} // End namespace coalescenceModels
} // End namespace Foam


#endif
