{

    Info << gMax(rho) << " " << gMin(rho) << endl;

    fvScalarMatrix eEqn
    (
        fvm::ddt(rho, e)
      - fvm::laplacian(turbulence->alphaEff(), e)
     ==
        fvm::SuSp(lambda/e, e)
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
