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
namespace breakupKernels
{
defineTypeNameAndDebug(CoulaloglouTavlarides, 0);
addToRunTimeSelectionTable
(
    breakupKernel,
    CoulaloglouTavlarides,
    dictionary
);

CoulaloglouTavlarides::CoulaloglouTavlarides(const dictionary &breakupDict,
                                             const phaseModel &dispersedPhase)
    :
      breakupKernel(breakupDict, dispersedPhase),
      impl_(
          breakupDict.lookupOrDefault<scalar>("c1",0.00487),
          breakupDict.lookupOrDefault<scalar>("c2",0.008),
          breakupDict.lookupOrDefault<scalar>("gamma",0.0),
          breakupDict.lookupOrDefault<scalar>("sigma",0.047)
          ),
      phase_(dispersedPhase),
      rhod_(dispersedPhase.rho()),
      epsilonName_("epsilon." + phase_.name())
{
}

dimensionedScalar CoulaloglouTavlarides::S(const dimensionedScalar &xi,
                                           label celli) const {

    const volScalarField& epsilonField = 
        phase_.U().mesh().lookupObject<volScalarField>(epsilonName_);

    dimensionedScalar epsilon(
        "epsilon", epsilonField.dimensions(), epsilonField[celli]);
    dimensionedScalar rho_d("rho_d", rhod_.dimensions(), rhod_[celli]);

    Info << "S(" << xi.value() << ", " << epsilon.value() << ", "
         << rho_d.value() << ")" << impl_.S(xi, rho_d, epsilon) << nl;
    return impl_.S(xi, rho_d, epsilon);
}

CoulaloglouTavlaridesImp::CoulaloglouTavlaridesImp(scalar c1, scalar c2,
                                                   scalar gamma, scalar sigma)
    :
      c1_(c1), c2_(c2), gamma_(gamma), sigma_(sigma)
{
//    Info << "CoulaloglouTavlaride model constants: \n{\n";
//    Info << "\tc1: " << c1_ << endl;
//    Info << "\tc2: " << c2_ << endl;
//    Info << "\tgamma: " << gamma_ << endl;
//    Info << "\tsigma: " << sigma_ << endl;
}

dimensionedScalar CoulaloglouTavlaridesImp::S(
        const dimensionedScalar &xi,
        const dimensionedScalar &rho_d,
        const dimensionedScalar &epsilon
        ) const
{
    return c1_ * pow(epsilon, 1.0/3.0)
            * exp
            (
                -c2_ * sigma_ * pow(1 + gamma_,2)
                /(
                    rho_d.value() * pow(epsilon.value(), 2.0/3.0)
                    * pow(xi.value(), 5.0/9.0)
                 )
            )
            /((1 + gamma_) * pow(xi, 1.0/3.0));
}

} //End namespace breakupKernels
} //End namespace Foam
