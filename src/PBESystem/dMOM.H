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
    Foam::dMOM

Description

   Method of moment solution with diamter as an internal variable. The
   implementation is based on:

   [1] Lo, S. and Rao, P. "Modelling of droplet breakup and coalescence in an
   oil-water pipeline." 6th International Conference of Multiphase Flow,
   2007

   [2] Lo, S. and Zhang, D. "Modelling of breakup and coalescence in vertical
   bubbly two-phase flows" 6th International Conference on CFD in Oil and Gas,
   Metalurgical and Process Inustries, 2008

   Assumptions are:
    * specific breakup and coalescence models (integrated analytically),
    * log-normal distribution of the diameter.


SourceFiles
    dMOM.C

\*---------------------------------------------------------------------------*/

#ifndef dMOMPBE_H
#define dMOMPBE_H

#include "PBEMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PBEMethods
{
/*---------------------------------------------------------------------------*\
                      Class dMOM Declaration
\*---------------------------------------------------------------------------*/

class dMOM : public PBEMethod
{
private:
 const dictionary dMOMDict_;
 const phaseModel& dispersedPhase_;
 const fvMesh& mesh_;

 std::vector<volScalarField> Sgammas_;
 std::vector<volScalarField> breakupSource_;

 //log-normal distribution parameters
 volScalarField log_d_bar_;
 volScalarField sigma_hat_;
 volScalarField number_density_;

 //Modelling parameters
 const scalar Nf_;
 const scalar C_alpha_;
 const scalar We_cr_;
 const scalar k_br_;
 const scalar Ca_cr_;
 const scalar k_cl_1_;
 const dimensionedScalar consistency_;
 const scalar interfacialTension_;
 const dimensionedScalar minDiameter_;

 //- Return a source term due to breakup
 tmp<volScalarField> breakupSourceTerm(label momenti);
 //- Return a source term due to coalescence
 tmp<volScalarField> coalescenceSourceTerm(label momenti);

 // Private member functions
 
 //- Update values of moments and diameter
 void updateMoments();
 //overwrites virtual solve from multiphaseSystem
 void solve();

public:
 TypeName("diameterMethodOfMoments");

 // Constructors

 //- Construct from components
 dMOM
 (
     const dictionary& pbeProperties,
     const phaseModel& phase
 );

 //- Dummy correct
 virtual void correct() override;

 ////- Return the Sauter-mean diameter
 virtual const volScalarField d() const override;



 //- Destructor
 virtual ~dMOM();

 // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PBEMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
