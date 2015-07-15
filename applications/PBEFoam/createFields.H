    Info<< "Creating twoPhaseSystem\n" << endl;

    twoPhaseSystem fluid(mesh, g);

    phaseModel& phase1 = fluid.phase1();
    phaseModel& phase2 = fluid.phase2();

    volScalarField& alpha1 = phase1;
    volScalarField& alpha2 = phase2;

    volVectorField& U1 = phase1.U();
    surfaceScalarField& phi1 = phase1.phi();
    surfaceScalarField& alphaPhi1 = phase1.alphaPhi();
    surfaceScalarField& alphaRhoPhi1 = phase1.alphaRhoPhi();

    volVectorField& U2 = phase2.U();
    surfaceScalarField& phi2 = phase2.phi();
    surfaceScalarField& alphaPhi2 = phase2.alphaPhi();
    surfaceScalarField& alphaRhoPhi2 = phase2.alphaRhoPhi();

    rhoThermo& thermo1 = phase1.thermo();
    rhoThermo& thermo2 = phase2.thermo();

    volScalarField& p = thermo1.p();

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.U()
    );
