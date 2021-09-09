// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff() // this approach takes bulk viscosity as equal to 0
);

UEqn().relax
(
    mesh.solutionDict().equationRelaxationFactor(U.select(pimple.finalIter()))
);

volScalarField rUA = 1.0/UEqn().A();

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn() == -fvc::grad(p),
        mesh.solutionDict().solver((U.select(pimple.finalIter())))
    );
}
else
{
    U = rUA*(UEqn().H() - fvc::grad(p));
    U.correctBoundaryConditions();
}
