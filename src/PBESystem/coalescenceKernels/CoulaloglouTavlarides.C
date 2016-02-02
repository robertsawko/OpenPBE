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

\*---------------------------------------------------------------------------*/

#include "CoulaloglouTavlarides.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {
namespace coalescenceKernels {
defineTypeNameAndDebug(CoulaloglouTavlaridesC, 0);
addToRunTimeSelectionTable(coalescenceKernel,
                           CoulaloglouTavlaridesC,
                           dictionary);

CoulaloglouTavlaridesC::CoulaloglouTavlaridesC(
    const dictionary &coalescenceDict, const phaseModel &dispersedPhase)
    : coalescenceKernel(coalescenceDict, dispersedPhase),
      impl_(coalescenceDict.lookupOrDefault<scalar>("c1", 0.00487),
            coalescenceDict.lookupOrDefault<scalar>("c2", 0.008),
            coalescenceDict.lookupOrDefault<scalar>("gamma", 0.0),
            coalescenceDict.lookupOrDefault<scalar>("sigma", 0.047)),
      phase_(dispersedPhase),
      rhod_(dispersedPhase.rho()),
      nud_(dispersedPhase.nu()),
      epsilonName_("epsilon." + phase_.name()) {}

void CoulaloglouTavlaridesC::update(){
    nud_ = dispersedPhase_.nu()();
}

scalar CoulaloglouTavlaridesC::S(const dimensionedScalar &xi1,
                                 const dimensionedScalar &xi2,
                                 label celli) const {

    static const volScalarField &epsilonField(
        phase_.U().mesh().lookupObject<volScalarField>(epsilonName_));

    return impl_.S(xi1.value(),
                   xi2.value(),
                   epsilonField[celli],
                   rhod_[celli],
                   nud_[celli],
                   phase_.U().mesh().V()[celli]);
}

CoulaloglouTavlaridesCImpl::CoulaloglouTavlaridesCImpl(scalar c1,
                                                       scalar c2,
                                                       scalar gamma,
                                                       scalar sigma)
    : c1_(c1), c2_(c2), gamma_(gamma), sigma_(sigma) {}

scalar CoulaloglouTavlaridesCImpl::S(scalar xi1,
                                     scalar xi2,
                                     scalar epsilon,
                                     scalar rhoc,
                                     scalar nuc,
                                     scalar Vcell) const {
    scalar xi1pow13 = pow(xi1, 1.0 / 3.0);
    scalar xi2pow13 = pow(xi2, 1.0 / 3.0);

    scalar frequency = c1_ * pow(epsilon, 1.0 / 3.0) *
                       pow(xi1pow13 + xi2pow13, 2) *
                       pow(pow(xi1, 2.0 / 9.0) + pow(xi2, 2.0 / 9.0), 0.5) /
                       (1 + gamma_) / Vcell;

    scalar rate = exp((-c2_ * rhoc * epsilon * (rhoc * nuc) /
                       (pow(sigma_, 2) * pow(1 + gamma_, 3)) *
                       pow(xi1pow13 * xi2pow13 / (xi1pow13 + xi2pow13), 4.0)));

    return frequency * rate;
}

} // End namespace coalescenceKernels
} // End namespace Foam
