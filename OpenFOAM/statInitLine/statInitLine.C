/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

Application
    zukSetCylinderEq

Description
     sets stationary flow to concentric spinning cylinders with specifed boundary temperatures

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OSspecific.H"
#include "fixedValueFvPatchFields.H"
#include "Random.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#include "setRootCase.H"

#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"


dimensionedScalar baseHs = fvc::domainIntegrate(rho*thermo.e()*t1/T);
dimensionedScalar initialMass = fvc::domainIntegrate(rho);
#include "medium.H"
thermo.correct();
dimensionedScalar finalMass = fvc::domainIntegrate(rho);

p *= initialMass/finalMass;
rho *= initialMass/finalMass;

dimensionedScalar statHs = fvc::domainIntegrate(rho*thermo.Cv()*T);

Info << "mass distortion: "<< 1. - (fvc::domainIntegrate(rho)).value()/initialMass.value() << endl;
Info << "base Hs:         "<< baseHs.value() << endl;
Info << "stationary Hs:   "<< statHs.value() << endl;


U.correctBoundaryConditions();
T.correctBoundaryConditions();
p.correctBoundaryConditions();
rho.correctBoundaryConditions();

U.write();
T.write();
p.write();
rho.write();
    
Info<< "End\n" << endl;

return 0;
}



//*************************************************************************//

