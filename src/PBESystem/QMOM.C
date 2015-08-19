/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include <string>
#include <sstream>
#include <cmath>

#include <boost/math/distributions/gamma.hpp>
#include <cassert>

#include "QMOM.H"
#include "addToRunTimeSelectionTable.H"

#include "fvScalarMatrix.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcAverage.H"
#include "fvm.H"
#include "Integrator.H"
#include "mathematicalConstants.H"
#include "MULES.H"
#include "Utility.H"
#include "MomentInversion.H"

namespace Foam {
namespace PBEMethods {

defineTypeNameAndDebug(QMOM, 0);
addToRunTimeSelectionTable(PBEMethod, QMOM, dictionary);
// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
using constant::mathematical::pi;

QMOM::QMOM(const dictionary &pbeProperties, const phaseModel &phase)
    : PBEMethod(pbeProperties, phase),
      QMOMDict_(pbeProperties.subDict("QMOMCoeffs")), dispersedPhase_(phase),
      mesh_(dispersedPhase_.U().mesh()),
      moments_(2 * readLabel(QMOMDict_.lookup("quadratureOrder"))),
      d_(IOobject("diameter",
                  mesh_.time().timeName(),
                  mesh_,
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE),
         mesh_,
         dimensionedScalar("diameter", dimLength, 0.0)),
      usingMULES_(QMOMDict_.lookupOrDefault<Switch>("usingMULES", false)) {
    forAll(moments_, momenti)
        moments_.set(momenti,
                     new volScalarField(IOobject("m" + std::to_string(momenti),
                                                 mesh_.time().timeName(),
                                                 mesh_,
                                                 IOobject::MUST_READ,
                                                 IOobject::AUTO_WRITE),
                                        mesh_));

    d_ = pow(6.0 / pi * moments_[1] / moments_[0], 1.0 / 3.0);
}

QMOM::~QMOM() {}
void QMOM::correct() {
    std::vector<volScalarField> mSources;

    for (label i = 0; i < moments_.size(); ++i) {
        mSources.emplace_back(IOobject("m" + std::to_string(i) + "Source",
                                       mesh_.time().timeName(),
                                       mesh_,
                                       IOobject::NO_READ,
                                       IOobject::AUTO_WRITE),
                              momentSourceTerm(i));
    }

    Info << "moment sources:" << endl;
    for (auto &mSource : mSources) {
        mSource.boundaryField() = 0;

        // fix for nan values
        forAll(mesh_.C(), celli) {
            auto &mSource_i = mSource[celli];
            if (std::isnan(mSource_i)) {
                mSource_i = 0.0;
            }
        }

        printAvgMaxMin(mesh_, mSource);
    }

    if (usingMULES_)
        solveWithMULES(mSources);
    else
        solveWithFVM(mSources);

    for (auto &moment : moments_) {
        moment =
            max(moment,
                dimensionedScalar(moment.name(), moment.dimensions(), SMALL));
        // TODO: print a warning message
        printAvgMaxMin(mesh_, moment);
    }

    d_ = pow(6.0 / pi * moments_[1] / moments_[0], 1.0 / 3.0);
    // clamp d_ here?

    printAvgMaxMin(mesh_, d_);
}

void QMOM::solveWithMULES(const std::vector<volScalarField> &S) {
    const fvMesh &mesh = moments_[0].mesh();
    const surfaceScalarField &phi = phase_.phi();
    PtrList<surfaceScalarField> phimkCorrs(moments_.size());

    forAll(moments_, k) {
        volScalarField &mk = moments_[k];

        phimkCorrs.set(k,
                       new surfaceScalarField(
                           "phi" + mk.name() + "Corr",
                           fvc::flux(phi, mk, "div(phi," + mk.name() + ')')));

        surfaceScalarField &phimkCorr = phimkCorrs[k];

        /*
        // Ensure that the flux at inflow BCs is preserved
        forAll(phiNkCorr.boundaryField(), patchi)
        {
            fvsPatchScalarField& phiNkCorrp =
                phiNkCorr.boundaryField()[patchi];

            if (!phiNkCorrp.coupled())
            {
                const scalarField& phi1p = phase1.phi().boundaryField()[patchi];
                const scalarField& alpha1p = alpha1.boundaryField()[patchi];

                forAll(phiAlphaCorrp, facei)
                {
                    if (phi1p[facei] < 0)
                    {
                        phiAlphaCorrp[facei] = alpha1p[facei]*phi1p[facei];
                    }
                }
            }
        }
        */

        // This limited assumed that moment cannot be greater than it's current
        // value + 20%
        scalar maxmk = 1.2 * max(mk).value();

        MULES::limit(1.0 / mesh.time().deltaT().value(),
                     geometricOneField(),
                     mk,
                     phi,
                     phimkCorr,
                     zeroField(), // implicit source?
                     S[k],
                     maxmk,
                     0,
                     3,
                     true);
    }

    MULES::limitSum(phimkCorrs);

    forAll(moments_, k) {
        volScalarField &mk = moments_[k];

        surfaceScalarField &phimk = phimkCorrs[k];
        phimk += upwind<scalar>(mesh, phi).flux(mk);

        MULES::explicitSolve(geometricOneField(),
                             mk,
                             phimk,
                             zeroField(), // implicit source?
                             S[k]);
    }
}

void QMOM::solveWithFVM(const std::vector<volScalarField> &S) {
    for (label i = 0; i < moments_.size(); ++i) {
        fvScalarMatrix mEqn(fvm::ddt(moments_[i]) +
                                fvm::div(dispersedPhase_.phi(), moments_[i]) ==
                            S[i]);
        mEqn.relax();
        mEqn.solve();
    }
}

const volScalarField QMOM::d() const { return d_; }

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

volScalarField QMOM::momentSourceTerm(label momenti) {
    volScalarField toReturn(
        IOobject("S",
                 mesh_.time().timeName(),
                 mesh_,
                 IOobject::NO_READ,
                 IOobject::NO_WRITE,
                 false),
        mesh_,
        dimensionedScalar("S", pow(dimVolume, momenti) / dimTime, 0));

    forAll(dispersedPhase_, celli) {
        Eigen::VectorXd momentVector(moments_.size());
        for (int i = 0; i < momentVector.size(); ++i)
            momentVector[i] = moments_[i][celli];

        auto quadrature = wheeler_inversion(momentVector);
        int N = quadrature.abcissas.size();

        double coalescenceSource{0.}, breakupSource{0.};

        for (int i = 0; i < N; ++i) {
            double xi_i = quadrature.abcissas[i];
            double innerCoalescenceDeathTerm{0.}, innerCoalescenceBirthTerm{0.};

            for (int j = 0; j < N; ++j) {
                double xi_j = quadrature.abcissas[j];

                innerCoalescenceDeathTerm +=
                    quadrature.weights[j] *
                    coalescence_->S(xi_i, xi_j, celli).value();

                innerCoalescenceBirthTerm +=
                    quadrature.weights[j] * pow(xi_i + xi_j, momenti) *
                    coalescence_->S(xi_i, xi_j, celli).value();
            }

            coalescenceSource +=
                quadrature.weights[i] *
                (innerCoalescenceBirthTerm / 2. -
                 pow(xi_i, momenti) * innerCoalescenceDeathTerm);

            breakupSource +=
                quadrature.weights[i] * breakup_->S(xi_i, celli) *
                (daughterParticleDistribution_->moment(xi_i, momenti).value() -
                 pow(xi_i, momenti));
        }

        toReturn[celli] = coalescenceSource + breakupSource;
    }

    return toReturn;
}

} // end namespace PBEMethods
}//end namespace Foam

// ************************************************************************* //
