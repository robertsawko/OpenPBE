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
    Foam::MOC

Description
    Incompressible multi-phase mixture with built in solution for the
    phase fractions with interface compression for interface-capturing.

    Derived from transportModel so that it can be unsed in conjunction with
    the incompressible turbulence models.

    Surface tension and contact-angle is handled for the interface
    between each phase-pair.

SourceFiles
    MOC.C

\*---------------------------------------------------------------------------*/

#ifndef MOCPBE_H
#define MOCPBE_H

#include "PBEMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PBEMethods
{
/*---------------------------------------------------------------------------*\
                      Class MOC Declaration
\*---------------------------------------------------------------------------*/

class MOC
:
    public PBEMethod
{
private:
 //- class sizes
 PtrList<volScalarField> classNumberDensity_;

 //- class velocities
 PtrList<volVectorField> classVelocity_;

 //- uniform spacing of classes
 dimensionedScalar deltaXi_;

 //- diameters for each class 
 std::vector<dimensionedScalar> xi_;

 bool usingMULES_;

 //TODO: combine volScalarField with
 
 //- Return a source term due to breakup
 volScalarField breakupSourceTerm(label classi);
 //- Return a source term due to coalescence
 volScalarField coalescenceSourceTerm(label classi);
 //- Return a source term due to coalescence
 volScalarField classSourceTerm(label classi);

 void solveWithFVM(const PtrList<volScalarField>& sources);
 void solveWithMULES(const PtrList<volScalarField>& sources);

 std::vector<scalar> breakupCache_;

protected:
 // Protected data

 // Protected member functions
public:
 TypeName("methodOfClasses");

 // Constructors

 //- Construct from components
 MOC
 (
     const dictionary& pbeProperties,
     const phaseModel& phase
 );

 void correct();
 const volScalarField d() const;


 //- Destructor
 virtual ~MOC()
 {}


 // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PBESystems
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
