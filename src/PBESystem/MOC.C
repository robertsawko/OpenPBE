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
#include <sstream>

#include "MOC.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "MULES.H"


namespace Foam
{
namespace PBEMethods
{
defineTypeNameAndDebug(MOC, 0);
addToRunTimeSelectionTable(PBEMethod, MOC, dictionary);
// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

MOC::MOC
(
    const dictionary& pbeProperties,
    const phaseModel& phase
)
:
    PBEMethod(pbeProperties, phase),
    MOCDict_(pbeProperties.subDict("MOCCoeffs")),
    numberOfClasses_(readLabel(MOCDict_.lookup("numberOfClasses"))),
    classNumberDensity_(numberOfClasses_),
    classVelocity_(numberOfClasses_),
    deltaXi_("deltaXi", dimVolume, readScalar(MOCDict_.lookup("xi1"))),
    xi_(numberOfClasses_),
    usingMULES_(MOCDict_.found("usingMULES"))
{
    Info << "Creating " << numberOfClasses_ << " class";
    //Taking pedantry one step too far!
    if(numberOfClasses_ > 1)
        Info << "es";
    Info << endl;

    //int phasei = 0;
    forAll(classNumberDensity_, i)
    {
        //TODO: Is it possible to do the same with OF string?
        std::stringstream className;
        std::stringstream xiName;
        className << "n" << i; 
        //TODO: There MUST be a way of doing it more easily
        xiName << "xi " << deltaXi_.value() * i; 
        classNumberDensity_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    className.str(),
                    phase.U().mesh().time().timeName(),
                    phase.U().mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                phase.U().mesh()
            )
        );
        classVelocity_.set(i, phase.U());
        xi_.set
        (
            i,
            new dimensionedScalar(xiName.str(), deltaXi_ * (i + 1))
        );
        Info<< className.str().c_str() << " has volume " << xi_[i]
            << endl;
    }
    //TODO: Add error checking classes must add to dispersed or to 1



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
            classNumberDensity_[i].dimensions()/dimTime,
            0
        ) 
    );

    forAll(xi_, classj)
    {
        coalescenceField -=
            coalescence_().S(xi_[i], xi_[classj])
            * classNumberDensity_[classj];
    }

    coalescenceField *= classNumberDensity_[i];

    //-1 in pairs account for zero-based numbering
    for(int j = 0; j < i; ++j)
    {
        coalescenceField +=
            0.5 * coalescence_().S(xi_[i - j - 1], xi_[j])
            * classNumberDensity_[i - j - 1] * classNumberDensity_[j];
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
    breakupField -=
        breakup_().S(xi_[i]) * classNumberDensity_[i];

    for (label j = i + 1; j < numberOfClasses_; ++j)
        // deltaXi comes from the application of mean value theorem on the
        // second integral see Kumar and Ramkrishna (1996) paper
        breakupField +=
            daughterParticleDistribution_().beta(xi_[i], xi_[j])
            * breakup_().S( xi_[j]) * classNumberDensity_[j] * deltaXi_;

    return breakupField;
}

void MOC::correct(){

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

    volScalarField sumNk
    (
        IOobject
        (
            "sumNk",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("sumNk", dimless, 0)
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

        Info<< Nk.name() << " volume average, min, max = "
            << Nk.weightedAverage(mesh.V()).value()
            << ' ' << min(Nk).value()
            << ' ' << max(Nk).value()
            << endl;

        sumNk += Nk;
    }
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
