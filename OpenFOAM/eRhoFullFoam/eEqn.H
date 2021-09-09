{

    // goes with viscous heating
    volScalarField symmGradU2("symmGradU2", twoSymm(Foam::T(fvc::grad(U))) && Foam::T(fvc::grad(U)));
    Info << max(turbulence->muEff()*symmGradU2) << endl;

    Info << gMax(rho) << " " << gMin(rho) << endl;

    fvScalarMatrix eEqn
    (
        fvm::ddt(rho, e)
      + fvm::div(phi, e)
      - fvm::laplacian(turbulence->alphaEff(), e)
     ==
       - fvm::SuSp
         (
             p*fvc::div(phi/fvc::interpolate(rho))/e,
             e
         )
     + fvm::SuSp(turbulence->muEff()*symmGradU2/e, e)
     + fvm::SuSp(lambda/e, e)
    );

    eEqn.relax
    (
        mesh.solutionDict().equationRelaxationFactor(e.select(pimple.finalIter()))
    );

    eEqn.solve();

    thermo.correct();

    // Recalculate density
    rho = thermo.rho();
    rho.correctBoundaryConditions();

}
