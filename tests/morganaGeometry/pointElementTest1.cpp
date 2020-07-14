/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "pointElement.hpp"
#include "geoShapes.h"

using namespace std;
using namespace boost::mpi;

int main(int argc, char *argv[])
{  
  pointElement<linearTetra> E1, E2;
  
  E1.setGeoId(10);
  E1.setPoint(1, point3d(0.0, 0.0, 0.0));
  E1.setPoint(2, point3d(1.0, 0.0, 0.0));
  E1.setPoint(3, point3d(0.0, 1.0, 0.0));
  E1.setPoint(4, point3d(0.0, 0.0, 1.0));
  
  E2.setGeoId(10);
  E2.setPoint(1, point3d(0.0, 0.0, 0.0));
  E2.setPoint(2, point3d(1.0, 0.0, 0.0));
  E2.setPoint(3, point3d(0.0, 1.0, 0.0));
  E2.setPoint(4, point3d(0.0, 0.0, 2.0));
  
  cout << "not equal : " << (E1 != E2) << endl;
  cout << "less      : " << (E1  < E2) << endl;
  
  cout << E1 << endl;
}
