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
    xi_(numberOfClasses_)
{
    Info << "Creating " << numberOfClasses_ << " class";
    //Taking pedantry one step too far!
    if(numberOfClasses_ > 1)
        Info << "es";
    Info << endl;

    scalar xi1 = readScalar(MOCDict_.lookup("xi1"));
    //int phasei = 0;
    forAll(classNumberDensity_, i)
    {
        //TODO: Is it possible to do the same with OF string?
        std::stringstream className;
        std::stringstream xiName;
        className << "n" << i; 
        //TODO: There MUST be a way of doing it more easily
        xiName << "xi " << xi1 * i; 
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
            new dimensionedScalar
            (
                xiName.str(),
                dimVolume,
                xi1 * i
            )
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

    //TODO: Fix the odd and even controversy
    for(int j = 0; j < i; ++j)
    {
        coalescenceField +=
            0.5 * coalescence_().S(xi_[i - j], xi_[j])
            * classNumberDensity_[i - j] * classNumberDensity_[j];
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
            classNumberDensity_[i].dimensions()/dimTime,
            0.0
        ) 
    );
    breakupField -=
        breakup_().S( xi_[i]) * classNumberDensity_[i];

    if (2 * i < numberOfClasses_)
        breakupField +=
            2 * breakup_().S( xi_[2 * i]) * classNumberDensity_[2 * i];

    return breakupField;
}

void MOC::correct()
{
    const fvMesh& m = classNumberDensity_[0].mesh();

    forAll(classNumberDensity_, k)
    {
        surfaceScalarField phi = fvc::interpolate(classVelocity_[k]) & m.Sf();
        volScalarField S = classSourceTerm(k);
        solve(
            fvm::ddt(classNumberDensity_[k])
            + 
            fvm::div(phi, classNumberDensity_[k], "div(U,nk)")
            ==
            S
        );
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
    for(int i = 1; i < numberOfClasses_; i++)
    {
        sum += volumeToDiameter(xi_[i]) * classNumberDensity_[i];
        norm += classNumberDensity_[i];
    }
    return sum / norm;
}
} //end of PBEMethods namespace
} //end of Foam namespace

// ************************************************************************* //
