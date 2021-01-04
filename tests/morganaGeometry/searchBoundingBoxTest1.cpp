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
#include <fstream>

#include "simpleFormats.hpp"
#include "typesInterface.hpp"

int main(int argc, char *argv[])
{
  //Data---------------------
  Real maxX =  100.0e-3;
  Real minX = -100.0e-3;
  
  Real maxY =  100.0e-3;
  Real minY = -100.0e-3;
  
  Real maxZ = 3.3 + 50e-3;
  Real minZ = 3.3 - 50e-3;
  
  Real h = 5e-3;
  string fileName = "box.msz";
  
  //Derived data-------------
  Real lX = maxX - minX;
  Real lY = maxY - minY;
  Real lZ = maxZ - minZ;
  
  UInt nx = lX / h;
  UInt ny = lY / h;
  UInt nz = lZ / h;
  
  //Create coordinates-------
  point3d P;
  sVect<point3d> points;
  
  for(UInt ix=0; ix <= nx; ++ix)
  {
    for(UInt iy=0; iy <= ny; ++iy)
    {
      for(UInt iz=0; iz <= nz; ++iz)
      {
	P.setX(minX + lX * (Real(ix) / Real(nx)) );
	P.setY(minY + lY * (Real(iy) / Real(ny)) );
	P.setZ(minZ + lZ * (Real(iz) / Real(nz)) );
	
	points.push_back(P);
      }
    }
  }
  
  cout << "Num points: " << points.size() << endl;
  
  ofstream myfile;
  myfile.open (fileName.c_str());
  
  myfile << points.size() << endl;
  for(UInt i=1; i <= points.size(); ++i)
  {
    myfile << points(i).getX() << " " << points(i).getY() << " " << points(i).getZ() << " " << h << endl;
  }
  
  myfile << 0 << endl;
  myfile.close();
}
