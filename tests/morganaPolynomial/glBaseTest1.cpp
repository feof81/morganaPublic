/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <cmath>
#include <iostream>

#include "glBase.h"
#include "polyDynamic.h"

// mpirun -np 2 ./bin/morgana



int main(int argc, char *argv[])
{
  glBase spectral;
  UInt nx = 2, ny = 0, nz = 0;
  UInt ix = 3, iy = 1, iz = 1;
  
  polyDynamicCard outPoly = spectral.getPolynomial(nx,ny,nz,ix,iy,iz);
  cout << "polynomial" << endl;
  cout << outPoly << endl << endl;
  
  polyDynamic polyVal;
  polyVal.setPolyDynamicCard(outPoly);
  
  point3d P(1.0, 0.0, 0.0);
  cout << "Val : " <<  polyVal.evaluate(P) << endl;
}