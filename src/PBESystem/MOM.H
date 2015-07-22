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
    m3_ should not be here since it is the same as alpha multiplied by constant
    factor (in theory)
    but convecting alpha with MULES give different results than convecting m2_
    which produce error in diameter d32_ = m3 / m2

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
 const phaseModel& dispersedPhase_;
 const fvMesh& mesh_;

 std::vector<volScalarField> moments_;

 volScalarField d_;

 //gamma distribution parameters
 volScalarField gamma_k_;
 volScalarField gamma_theta_;
 volScalarField gamma_c0_;

 //- Return a source term due to breakup
 tmp<volScalarField> breakupSourceTerm(label momenti);
 //- Return a source term due to coalescence
 tmp<volScalarField> coalescenceSourceTerm(label momenti);
 //- Return a source term due to coalescence
 tmp<volScalarField> momentSourceTerm(label momenti);

 // Private member functions
 
 //- Update values of moments and diameter
 void updateMoments();
 //overwrites virtual solve from multiphaseSystem
 void solve();

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
 virtual ~MOM();

 // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PBEMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
