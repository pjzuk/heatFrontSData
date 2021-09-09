IOdictionary stationaryParameters
(
    IOobject
    (
        "stationaryParameters",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

scalar r1
(
    stationaryParameters.lookupOrDefault("x1",0.)
);
scalar r2
(
    stationaryParameters.lookupOrDefault("x2",1.)
);
scalar t1
(
    stationaryParameters.lookupOrDefault("t1",300.)
);
scalar t2
(
    stationaryParameters.lookupOrDefault("t2",300.)
);

volVectorField Ubc
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
volVectorField U
(
    IOobject
    (
        "U_st",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Ubc
);

volScalarField Tbc
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
volScalarField T
(
    IOobject
    (
        "T_st",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Tbc
);

volScalarField pbc
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
volScalarField p
(
    IOobject
    (
        "p_st",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pbc
);

volScalarField rho
(
    IOobject
    (
        "rho_st",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimMass/dimLength/dimLength/dimLength
);


autoPtr<basicPsiThermo> pThermo
(
    basicPsiThermo::New(mesh)
);
basicPsiThermo& thermo = pThermo();

rho=thermo.rho();

#   include "compressibleCreatePhi.H"
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

