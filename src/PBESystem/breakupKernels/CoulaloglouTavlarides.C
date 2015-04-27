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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

CoulaloglouTavlarides::CoulaloglouTavlarides
(
    const dictionary& breakupDict,
    const phaseModel& continuousPhase,
    const phaseModel& dispersedPhase
)
:
   breakupKernel(breakupDict,continuousPhase,dispersedPhase),
    c1_( breakupDict.lookupOrDefault<scalar>("c1",0.00487) ),
    c2_( breakupDict.lookupOrDefault<scalar>("c2",0.008) ),
    gamma_( breakupDict.lookupOrDefault<scalar>("gamma",0.0) ),
    sigma_( breakupDict.lookupOrDefault<scalar>("sigma",0.047) )
   //TODO this lookup does not work; there is no object of type mEuler::RASModel
   //in the registry
{
    Info << "CoulaloglouTavlaride model constants: \n{\n"; 
    Info << "\tc1: " << c1_ << endl;
    Info << "\tc2: " << c2_ << endl;
    Info << "\tgamma: " << gamma_ << endl;
    Info << "\tsigma: " << sigma_ << endl;

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
//
CoulaloglouTavlarides::~CoulaloglouTavlarides()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> CoulaloglouTavlarides::S(const volScalarField& xi) const
{
    const multiphaseTurbulence::turbulenceModel& turbulence
    (
       continuousPhase_.mesh().lookupObject
       <multiphaseTurbulence::turbulenceModel>("turbulenceModel")
    );

    volScalarField epsilon = turbulence.epsilon();
    //volScalarField k = turbulence.k();
    dimensionedScalar epsilonDim("epsilonDim", epsilon.dimensions(), 1.0);
    dimensionedScalar xiDim("xiDim", xi.dimensions(), 1.0);

    dimensionedScalar sigma("sigma", dimensionSet( 1, 0, -2, 0, 0), sigma_);
    dimensionedScalar c1("c1", dimensionSet( 0, 1, 0, 0, 0), c1_);

    const dimensionedScalar rhod = dispersedPhase_.rho();
        
    //Info << breakupKernel.weightedAverage(k.mesh().V()) << endl;
    volScalarField xiSafe = max(xi, SMALL * xiDim);
    volScalarField epsilonSafe = max(epsilon, SMALL * epsilonDim);


    return tmp<volScalarField>
        (
            c1 * pow( epsilonSafe , 1.0 / 3.0 )
            * exp
            (
                - c2_ * sigma * pow( 1 + gamma_, 2)
                / (
                    rhod * pow( epsilonSafe, 2.0 / 3.0) 
                    * pow(xiSafe, 5.0 / 3.0)
                )
            )
            /( (1 + gamma_) * pow(xiSafe, 2.0 / 3.0) )
        );
}

} //End namespace breakupKernels
} //End namespace Foam
// ************************************************************************* //
