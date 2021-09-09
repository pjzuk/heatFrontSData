rho = thermo.rho();

rho.relax();
Info<< "rho max/min : " << max(rho).value()
    << " " << min(rho).value() << endl;

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity

dimensionedScalar currentMass = fvc::domainIntegrate(rho);

if (closedVolume)
{
    p *= initialMass/currentMass;
}

