/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::MOM

Description
    Incompressible multi-phase mixture with built in solution for the
    phase fractions with interface compression for interface-capturing.

    Derived from transportModel so that it can be unsed in conjunction with
    the incompressible turbulence models.

    Surface tension and contact-angle is handled for the interface
    between each phase-pair.

SourceFiles
    MOM.C

\*---------------------------------------------------------------------------*/

#ifndef MOMPBE_H
#define MOMPBE_H

#include "PBEMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PBEMethods
{
/*---------------------------------------------------------------------------*\
                      Class MOM Declaration
\*---------------------------------------------------------------------------*/

class MOM : public PBEMethod
{
private:
 const dictionary MOMDict_;
 const phaseModel& phase_;
 //first three moments
 volScalarField m0_;
 volScalarField m1_;
 volScalarField m2_;
 //m3_ should not be here since it is the same as alpha multiplied by constant
 //factor (in theory)
 //but convecting alpha with MULES give different results than convecting m2_
 //which produce error in diameter d32_ = m3 / m2
 volScalarField m3_;

 //source terms are stored so they do not have to be recalculated in every
 //iteration of alpha sub-cycle
 volScalarField m0Source_;
 volScalarField m1Source_;
 volScalarField m2Source_;
 volScalarField m3Source_;

 //Sauter mean diameter
 volScalarField d32_;
 //expected diameter TODO: it's only for debugging purposes; to be removed later
 volScalarField expectedD_;

 //gamma distribution parameters
 volScalarField gamma_alpha_;
 volScalarField gamma_beta_;
 volScalarField gamma_c0_;

 //scaling factor for diameter
 const dimensionedScalar scaleD_;
 const scalar scaleM3_;
 //number of daughter droplets
 scalar Nf_;
 //maximum diameter
 dimensionedScalar maxD_;
 //minimum diameter
 dimensionedScalar minD_;
 //minimum alpha parameter
 dimensionedScalar minGammaAlpha_;
 //maximum beta parameter
 scalar maxGammaBeta_;
 //parameters of numerical integration
 scalar maxDiameterMultiplicator_;
 //how many step to perform in numerical integration of source terms
 label integrationSteps_;

 // values of breakup and coalescence kernels
 PtrList<volScalarField> bList_;
 PtrList<volScalarField> cList_;
 // values of gamma distribution
 PtrList<volScalarField> gammaList_;

 //used to copy flux calculated in solveAlphas() to updateMoments() method
 //surfaceScalarField limitedFlux_;

 //- Return a source term due to breakup
 tmp<volScalarField> breakupSourceTerm(label momenti);
 //- Return a source term due to coalescence
 tmp<volScalarField> coalescenceSourceTerm(label momenti);
 //- Return a source term due to coalescence
 tmp<volScalarField> momentSourceTerm(label momenti);

 //- Sub-cycle number
 //label nAlphaSubCycles_;

 // Private member functions
 
 //- Update values of moments and diameter
 void updateMoments();
 //overwrites virtual solve from multiphaseSystem
 void solve();

 // Private member functions
 
 //- Return gamma function values of given diameter field
 tmp<volScalarField> gamma(const volScalarField& x);
 
 //- Return gamma distribution values for a given diameter
 tmp<volScalarField> gammaDistribution(const dimensionedScalar xi);

protected:
 // Protected data

 // Protected member functions

 //! Solve the PBE and convect the classes
 /*!
  *  Method of classes will only add the source terms
  *  *\author Robert Sawko
  *\return void
  */
 void solveAlphas();
public:
 TypeName("methodOfMoments");

 // Constructors

 //- Construct from components
 MOM
 (
     const dictionary& pbeProperties,
     const phaseModel& phase
 );

 //- Dummy correct
 virtual void correct() override;

 ////- Return the Sauter-mean diameter
 virtual const volScalarField d() const override;



 //- Destructor
 virtual ~MOM()
 {}


 // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PBEMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
