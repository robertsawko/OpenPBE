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
    Foam::CoulaloglouTavlarides

Description

SourceFiles
    CoulaloglouTavlaridesI.H
    CoulaloglouTavlarides.C
    CoulaloglouTavlaridesIO.C

\*---------------------------------------------------------------------------*/

#ifndef CoulaloglouTavlarides_H
#define CoulaloglouTavlarides_H

#include "breakupKernel.H"

namespace Foam
{
namespace breakupKernels
{

class CoulaloglouTavlaridesImp
{
public:
    CoulaloglouTavlaridesImp(scalar c1, scalar c2, scalar gamma, scalar sigma);

    scalar S(scalar xi,
                        scalar rho_d,
                        scalar epsilon) const;

private:
    const scalar c1_;
    const scalar c2_;
    const scalar gamma_;
    const scalar sigma_;
};

class CoulaloglouTavlarides
:
    public breakupKernel
{
 public:
  TypeName("CoulaloglouTavlarides")

  CoulaloglouTavlarides
  (
      const dictionary& breakupDict,
      const phaseModel& dispersedPhase
  );

  virtual ~CoulaloglouTavlarides() = default;

  scalar S(const dimensionedScalar& xi, label celli) const override;

private:
  CoulaloglouTavlaridesImp impl_;
  const phaseModel& phase_;
  const volScalarField& rhod_;
};


} // End namespace breakupModels
} // End namespace Foam

#endif
