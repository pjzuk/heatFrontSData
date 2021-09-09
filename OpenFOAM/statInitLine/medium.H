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
    medium.H

Description
    template for porous medium shape
\*---------------------------------------------------------------------------*/


// volVectorField containig coordinates of cell centers
const volVectorField& Ccells = U.mesh().C();
scalar r = r1;

forAll(Ccells, celli)
{
   r = Ccells[celli].x() - r1;
   T[celli] = t1+ r/(r2-r1)*(t2-t1);
}
T.correctBoundaryConditions();
rho = rho*Tbc/T;
