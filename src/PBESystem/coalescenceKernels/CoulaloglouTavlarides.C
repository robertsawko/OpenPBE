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

namespace Foam
{
namespace coalescenceKernels
{
defineTypeNameAndDebug(CoulaloglouTavlaridesC, 0);
addToRunTimeSelectionTable
(
    coalescenceKernel,
    CoulaloglouTavlaridesC,
    dictionary
);

CoulaloglouTavlaridesC::CoulaloglouTavlaridesC
(const dictionary& coalescenceDict,
    const phaseModel& dispersedPhase
)
:
    coalescenceKernel(coalescenceDict, dispersedPhase),
    impl_( coalescenceDict.lookupOrDefault<scalar>("c1",0.00487),
           coalescenceDict.lookupOrDefault<scalar>("c2",0.008),
           coalescenceDict.lookupOrDefault<scalar>("gamma",0.0),
           coalescenceDict.lookupOrDefault<scalar>("sigma",0.047) ),
    phase_(dispersedPhase),
    rhod_(dispersedPhase.rho()),
    nud_(dispersedPhase.nu()),
    epsilonName_("epsilon." + phase_.name())
{
}

dimensionedScalar CoulaloglouTavlaridesC::S(const dimensionedScalar &xi1,
                                            const dimensionedScalar &xi2,
                                            label celli) const {

    const volScalarField& epsilonField (
        phase_.U().mesh().lookupObject<volScalarField>(epsilonName_));

    dimensionedScalar epsilon(
        "epsilon", epsilonField.dimensions(), epsilonField[celli]),
        rhod("rhod", rhod_.dimensions(), rhod_[celli]),
        nud("nud", nud_.dimensions(), nud_[celli]),
        Vcell("V", dimVolume, phase_.U().mesh().V()[celli]);

    return impl_.S(xi1, xi2, epsilon, rhod, nud, Vcell);
}

CoulaloglouTavlaridesCImpl::CoulaloglouTavlaridesCImpl(scalar c1,
                                                       scalar c2,
                                                       scalar gamma,
                                                       scalar sigma)
    :
      c1_(c1), c2_(c2), gamma_(gamma), sigma_(sigma)
{}

dimensionedScalar CoulaloglouTavlaridesCImpl::S
(
        const dimensionedScalar& xi1,
        const dimensionedScalar& xi2,
        const dimensionedScalar& epsilon,
        const dimensionedScalar& rhoc,
        const dimensionedScalar& nuc,
        const dimensionedScalar& Vcell
) const
{
    dimensionedScalar xi1pow13 = pow(xi1, 1.0 / 3.0);
    dimensionedScalar xi2pow13 = pow(xi2, 1.0 / 3.0);

    dimensionedScalar frequency =
        c1_ * pow(epsilon, 1.0 / 3.0) *
        pow(xi1pow13 + xi2pow13, 2) *
        pow
        (
            pow(xi1, 2.0 / 9.0) + pow(xi2, 2.0 / 9.0),
            0.5
        )
        / (1 + gamma_) / Vcell;


    dimensionedScalar rate =
        exp
        (
            (- c2_ * rhoc.value() * epsilon.value() * (rhoc * nuc).value()
            / ( pow(sigma_, 2) * pow(1 + gamma_, 3) )
            * pow
                (
                    xi1pow13 * xi2pow13 / (xi1pow13 + xi2pow13),
                    4.0
                ).value())
        );

    return frequency * rate;
}

} //End namespace coalescenceKernels
} //End namespace Foam
