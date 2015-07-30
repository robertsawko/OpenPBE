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
    Foam::CoulaloglouTavlaridesC

Description

SourceFiles
    CoulaloglouTavlaridesCI.H
    CoulaloglouTavlaridesC.C
    CoulaloglouTavlaridesCIO.C

\*---------------------------------------------------------------------------*/

#ifndef CoulaloglouTavlaridesC_H
#define CoulaloglouTavlaridesC_H

#include "coalescenceKernel.H"

namespace Foam
{
namespace coalescenceKernels
{

class CoulaloglouTavlaridesCImpl
{
public:
 CoulaloglouTavlaridesCImpl(scalar c1, scalar c2, scalar gamma, scalar sigma);

 dimensionedScalar S
 (
         const dimensionedScalar& xi1,
         const dimensionedScalar& xi2,
         const dimensionedScalar& epsilon,
         const dimensionedScalar& rhod,
         const dimensionedScalar& nud
 ) const;

private:
 const scalar c1_;
 const scalar c2_;
 const scalar gamma_;
 const scalar sigma_;

};

class CoulaloglouTavlaridesC
:
    public coalescenceKernel
{
 public:  
  TypeName("CoulaloglouTavlaridesC")

  CoulaloglouTavlaridesC
  (
      const dictionary& coalescenceDict,
      const phaseModel& dispersedPhase
  );

  virtual ~CoulaloglouTavlaridesC() = default;

  dimensionedScalar S(
      const dimensionedScalar& xi1,
      const dimensionedScalar& xi2) const override;

private:
  CoulaloglouTavlaridesCImpl impl_;
};

} // End namespace coalescenceModels
} // End namespace Foam

#endif
