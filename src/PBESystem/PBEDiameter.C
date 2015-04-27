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

#include "PBEDiameter.H"
#include "addToRunTimeSelectionTable.H"
//#include "alphaContactAngleFvPatchScalarField.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(PBEDiameter, 0);

    addToRunTimeSelectionTable
    (
        diameterModel,
        PBEDiameter,
        dictionary
    );
}
}
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::PBEDiameter::PBEDiameter
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterModel(diameterProperties, phase),
    method_(
        PBEMethod::New(
            diameterProperties,
            phase
        )
    ),
    d_
    (
        IOobject
        (
            IOobject::groupName("d", phase.name()),
            phase_.U().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        method_->d()
    )
    //continuousPhaseName_(pbeDict_.lookup("continuousPhase")),
    //dispersedPhaseName_(pbeDict_.lookup("dispersedPhase")),
    //breakup_(breakupKernel::New(pbeDict_, continuousPhase(), dispersedPhase())) 
{
    //TODO add some error checking. Do phases exist etc.
}

/*
namespace Foam
{
const phaseModel& PBEDiameter::continuousPhase() const
{
    return *phases_[continuousPhaseName_];
}
const phaseModel& PBEDiameter::dispersedPhase() const
{
    return *phases_[dispersedPhaseName_];
}
phaseModel& PBEDiameter::continuousPhase()
{
    return *phases_[continuousPhaseName_];
}
phaseModel& PBEDiameter::dispersedPhase()
{
    return *phases_[dispersedPhaseName_];
}

*/
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
