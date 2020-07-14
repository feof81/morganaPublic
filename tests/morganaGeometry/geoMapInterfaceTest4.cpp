/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "time.h"
#include <cmath>
#include <iostream>

#include "typesInterface.hpp"
#include "geoMapInterface.hpp"


/*! PERFORMACE PERFORMACE */
int main(int argc, char *argv[])
{
  time_t start, end;
  
  geoMapInterface<linearTriangle> interface;
  
  point3d Y;
  point3d Y1(1.0,1.0,0.0);
  point3d Y2(1.0,2.0,0.0);
  point3d Y3(2.0,1.0,0.0);
  
  sVect<point3d> points(3);  
  points(1) = Y1;
  points(2) = Y2;
  points(2) = Y3;
  
  UInt n=1e8;
  
  
  //Position test
  point3d Z;
  
  cout << "Eval" << " "; time(&start); cout << endl;
  for(UInt i=1; i <= n; ++i)
  {
    Y = Y1 * (Real(i) / Real(n)) + Y2 * (1.0 - Real(i) / Real(n));
    Z = interface.getPosition(points,Y);
  }  
  time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;
  
  
  //Jacobian test
  tensor3d T;
  
  cout << "Jacobian" << " "; time(&start); cout << endl;
  for(UInt i=1; i <= n; ++i)
  {
    Y = Y1 * (Real(i) / Real(n)) + Y2 * (1.0 - Real(i) / Real(n));
    T = interface.getGradient(points,Y);
  } 
  time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;
  
  
  //Jacobian det test
  Real D;
  
  cout << "Det" << " "; time(&start); cout << endl;
  for(UInt i=1; i <= n; ++i)
  {
    Y = Y1 * (Real(i) / Real(n)) + Y2 * (1.0 - Real(i) / Real(n));
    D = interface.getGradientDet(points,Y);
  } 
  time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;
}
