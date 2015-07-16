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
    Foam::QMOM

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
    QMOM.C

\*---------------------------------------------------------------------------*/

#ifndef QMOMPBE_H
#define QMOMPBE_H

#include "PBEMethod.H"
#include <array>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
using Eigen::VectorXd;
using Eigen::MatrixXd;

struct pair{
    VectorXd abcissas, weights;
};

pair wheeler_inversion(const VectorXd& moments)
{
    int N = moments.size() / 2;  // number of nodes of the quadrature approximation

    MatrixXd sigma = MatrixXd::Zero(N + 1, 2 * N);
    sigma.row(1) = moments.head(2 * N);

    VectorXd a = VectorXd::Zero(N);
    VectorXd b = VectorXd::Zero(N);
    a(0) = moments(1) / moments(0);
    b(0) = moments(0);  // This value is insignificant as it's not being used.

    for(int k = 1; k < N; ++k){
        int Ncol = 2 * (N - k);
        sigma.block(k + 1, k, 1, Ncol)
            = sigma.block(k, k + 1, 1, Ncol)
            - a(k - 1) * sigma.block(k, k, 1, Ncol)
            - b(k - 1) * sigma.block(k - 1, k, 1, Ncol);

        a(k) = -sigma(k, k) / sigma(k, k - 1)
            +
            sigma(k+1, k+1) / sigma(k+1, k);

        b(k) = sigma(k+1, k) / sigma(k, k-1);
    }
    
    VectorXd b_diag = - b.array().abs().sqrt(); //-Eigen::sqrt(Eigen::abs(b));
    //TODO: find a tridiagonal representation
    MatrixXd jacobi = a.asDiagonal();
    for(int i = 0; i < b_diag.size() - 1; ++i){
        jacobi(i, i + 1) = b_diag(i + 1);
        jacobi(i + 1, i) = b_diag(i + 1);
    }

    // eigval, eigvec = np.linalg.eig(jacobi)
    //VectorXd eivals = A.selfadjointView<Eigen::Lower>().eigenvalues();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenresult(jacobi);
    VectorXd abcissas = eigenresult.eigenvalues();
    VectorXd weights = moments(0) * eigenresult.eigenvectors().row(0).array().pow(2);

    return pair{abcissas, weights};
};

namespace Foam
{
namespace PBEMethods
{
/*---------------------------------------------------------------------------*\
                      Class QMOM Declaration
\*---------------------------------------------------------------------------*/

class QMOM : public PBEMethod
{
private:
 const dictionary QMOMDict_;
 const phaseModel& dispersedPhase_;
 const fvMesh& mesh_;

 std::vector<volScalarField> moments_;

 volScalarField d_;

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

 void printAvgMaxMin(const volScalarField &v) const;
public:
 TypeName("methodOfMoments");

 // Constructors

 //- Construct from components
 QMOM
 (
     const dictionary& pbeProperties,
     const phaseModel& phase
 );

 //- Dummy correct
 virtual void correct() override;

 ////- Return the Sauter-mean diameter
 virtual const volScalarField d() const override;



 //- Destructor
 virtual ~QMOM();

 // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PBEMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
