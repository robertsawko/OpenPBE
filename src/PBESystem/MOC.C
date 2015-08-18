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
#include "fvCFD.H"

#include <string>

#include "MOC.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "MULES.H"
#include "Utility.H"


namespace Foam
{
namespace PBEMethods
{
defineTypeNameAndDebug(MOC, 0);
addToRunTimeSelectionTable(PBEMethod, MOC, dictionary);
// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

MOC::MOC(const dictionary &pbeProperties, const phaseModel &phase)
    : PBEMethod(pbeProperties, phase),
      MOCDict_(pbeProperties.subDict("MOCCoeffs")),
      numberOfClasses_(readLabel(MOCDict_.lookup("numberOfClasses"))),
      classNumberDensity_(numberOfClasses_), classVelocity_(numberOfClasses_),
      deltaXi_("deltaXi", dimVolume, readScalar(MOCDict_.lookup("xi1"))),
      xi_(numberOfClasses_),
      usingMULES_(MOCDict_.lookupOrDefault<Switch>("usingMULES", false)),
      breakupCache_(numberOfClasses_*phase.size())
{
    Info << "Creating " << numberOfClasses_ << " class";
    // Taking pedantry one step too far!
    if (numberOfClasses_ > 1)
        Info << "es";
    Info << endl;

    // int phasei = 0;
    forAll(classNumberDensity_, i) {
        classNumberDensity_.set(
            i,
            new volScalarField(
                   IOobject(
                       "n" + std::to_string(i),
                       phase.U().mesh().time().timeName(),
                       phase.U().mesh(),
                       IOobject::MUST_READ,
                       IOobject::AUTO_WRITE),
                   phase.U().mesh()));
        classVelocity_.set(i, phase.U());
        xi_.set(i, new dimensionedScalar(
            "xi" + std::to_string(i), deltaXi_ * (i + 1)));
        Info << i << " has volume " << xi_[i] << endl;
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

volScalarField MOC::classSourceTerm(label i)
{
    return coalescenceSourceTerm(i) + breakupSourceTerm(i);
}
volScalarField MOC::coalescenceSourceTerm(label i)
{
    volScalarField coalescenceField
    (
        IOobject
        (
            "Scoal",
            classNumberDensity_[i].mesh().time().timeName(),
            classNumberDensity_[i].mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        classNumberDensity_[i].mesh(),
        //TODO Is the dimension correct?
        dimensionedScalar
        (
            "Scoal", 
            classNumberDensity_[i].dimensions() / dimTime,
            0
        ) 
    );

    forAll(phase_, celli){
        // the upper limit stops from removing mass through coalescence that
        // results in drops that are outside of PBE domain
        for(int j = 0; j < classNumberDensity_.size() - i - 1; ++j){
            coalescenceField[celli] -=
                coalescence_().S(xi_[i], xi_[j], celli).value()
                * classNumberDensity_[j][celli];
        }
        coalescenceField[celli] *= classNumberDensity_[i][celli];
        //-1 in pairs account for zero-based numbering
        for(int j = 0; j < i; ++j){
            coalescenceField[celli] +=
                0.5 * coalescence_().S(xi_[i - j - 1], xi_[j], celli).value()
                * classNumberDensity_[i - j - 1][celli]
                * classNumberDensity_[j][celli];
        }
    }
    return coalescenceField;

}

volScalarField MOC::breakupSourceTerm(label i)
{
    volScalarField breakupField(
        IOobject
        (
            "Sbr",
            classNumberDensity_[i].mesh().time().timeName(),
            classNumberDensity_[i].mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        classNumberDensity_[i].mesh(),
        //TODO Is the dimension correct?
        dimensionedScalar
        (
            "Sbr", 
            classNumberDensity_[i].dimensions() / dimTime,
            0.0
        ) 
    );

    auto phaseSize = phase_.size();

    forAll(phase_, celli){
        breakupField[celli] -=
            breakupCache_[i*phaseSize+celli] * classNumberDensity_[i][celli];

        for (label j = i + 1; j < numberOfClasses_; ++j){
            // deltaXi comes from the application of mean value theorem on the
            // second integral see Kumar and Ramkrishna (1996) paper
            breakupField[celli] +=
                daughterParticleDistribution_().beta(xi_[i], xi_[j]).value()
                * breakupCache_[j*phaseSize+celli]
                * classNumberDensity_[j][celli] * deltaXi_.value();
        }
    }
    return breakupField;
}

void MOC::correct(){

    auto phaseSize = phase_.size();

    for (int i=0; i<numberOfClasses_; ++i){
        for (int j=0; j<phaseSize; ++j)
            breakupCache_[i*phaseSize+j] = breakup_->S(xi_[i],j).value();
    }

    PtrList<volScalarField> S(numberOfClasses_);
    forAll(classNumberDensity_, k){
        S.set(k, classSourceTerm(k));
    }

    if(usingMULES_)
        solveWithMULES(S);
    else
        solveWithFVM(S);
}

void MOC::solveWithFVM(const PtrList<volScalarField>& S){
    const fvMesh& mesh = classNumberDensity_[0].mesh();

    forAll(classNumberDensity_, k){
        surfaceScalarField phi = fvc::interpolate(classVelocity_[k]) & mesh.Sf();
        fvScalarMatrix nEqn
        (
            fvm::ddt(classNumberDensity_[k])
            + 
            fvm::div(phi, classNumberDensity_[k], "div(U,nk)")
            ==
            S[k]
        );
        nEqn.relax();
        nEqn.solve();
    }
}

void MOC::solveWithMULES(const PtrList<volScalarField>& S){
    const fvMesh& mesh = classNumberDensity_[0].mesh();
    const surfaceScalarField& phi = phase_.phi();
    PtrList<surfaceScalarField> phiNkCorrs(classNumberDensity_.size());

    forAll(classNumberDensity_, k)
    {
        volScalarField& Nk = classNumberDensity_[k];

        phiNkCorrs.set
        (
            k,
            new surfaceScalarField
            (
                "phi" + Nk.name() + "Corr",
                fvc::flux
                (
                    phi,
                    Nk,
                    "div(phi," + Nk.name() + ')'
                )
            )
        );

        surfaceScalarField& phiNkCorr = phiNkCorrs[k];

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
        scalar maxNk = max(phase_ * mesh.V() / xi_[k]).value();

        MULES::limit
        (
            1.0 / mesh.time().deltaT().value(),
            geometricOneField(),
            Nk,
            phi,
            phiNkCorr,
            zeroField(), //implicit source?
            S[k], 
            maxNk,
            0,
            3,
            true
        );
    }

    MULES::limitSum(phiNkCorrs);

    volScalarField volumeSum
    (
        IOobject
        (
            "sumVolumeNk",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("sumVolumeNk", xi_[0].dimensions(), 0)
    );

    forAll(classNumberDensity_, k)
    {
        volScalarField& Nk = classNumberDensity_[k];

        surfaceScalarField& phiNk = phiNkCorrs[k];
        phiNk += upwind<scalar>(mesh, phi).flux(Nk);

        MULES::explicitSolve
        (
            geometricOneField(),
            Nk,
            phiNk,
            zeroField(), //implicit source?
            S[k]
        );

        volumeSum += xi_[k] * Nk;
    }

    volScalarField alphaPBE = phase_;
    alphaPBE.internalField() = volumeSum.internalField() / mesh.V();
    volScalarField continuityError("continuityError", mag(alphaPBE - phase_));

    Info<< "Max continuity error between advection and PBE: "
        << max(continuityError).value() << endl;
        
}

inline dimensionedScalar volumeToDiameter(const dimensionedScalar& v)
{
    return pow(6.0 / constant::mathematical::pi * v, 1.0/3.0);
}
const volScalarField MOC::d() const
{
    volScalarField sum = volumeToDiameter(xi_[0]) * classNumberDensity_[0];
    volScalarField norm = classNumberDensity_[0];
    for(int i = 1; i < numberOfClasses_; ++i)
    {
        sum += volumeToDiameter(xi_[i]) * classNumberDensity_[i];
        norm += classNumberDensity_[i];
    }
    return sum / norm;
}
} //end of PBEMethods namespace
} //end of Foam namespace

// ************************************************************************* //
