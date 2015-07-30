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
           coalescenceDict.lookupOrDefault<scalar>("sigma",0.047) )
{
}

dimensionedScalar CoulaloglouTavlaridesC::S(const dimensionedScalar &xi1,
                                            const dimensionedScalar &xi2) const
{
    scalar epsilon(1), nud(1), rhod(1);
    return impl_.S(xi1, xi2, epsilon, rhod, nud);
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
        const dimensionedScalar& rhod,
        const dimensionedScalar& nud
) const
{
    dimensionedScalar frequency =
        c1_ * pow(epsilon, 1.0 / 3.0) * pow(xi1 + xi2, 2) *
        pow
        (
            pow(xi1, 2.0 / 3.0) + pow(xi2, 2.0 / 3.0),
            0.5
        )
        / (1 + gamma_);


    dimensionedScalar rate =
        exp
        (
            (- c2_ * rhod * epsilon * (rhod * nud)
            / ( pow(sigma_,2) * pow( 1 + gamma_,3) )
            * pow
                (
                    xi1 * xi2 / (xi1 + xi2),
                    4.0
                ).value())
        );


    return frequency * rate;
}

} //End namespace coalescenceKernels
} //End namespace Foam
