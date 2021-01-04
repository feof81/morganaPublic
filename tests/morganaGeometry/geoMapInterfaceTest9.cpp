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

#include "typesInterface.hpp"
#include "geoMapInterface.hpp"

int main(int argc, char *argv[])
{
  geoMapInterface<linearHexa> interface;
  
  cout << "refShape   : " << interface.getReferenceShapes() << endl;
  cout << "refGeo     : " << interface.getReferenceGeometry() << endl;
  cout << "dim        : " << interface.getNDim() << endl;
  cout << "numFaces   : " << interface.getNumFaces()    << " " << interface.getNumGeoTypes(FACE) << endl;
  cout << "numEdges   : " << interface.getNumEdges()    << " " << interface.getNumGeoTypes(EDGE) << endl;
  cout << "numVertices: " << interface.getNumVertices() << " " << interface.getNumGeoTypes(VERTEX) << endl;

  
  point3d Y(0.5,0.5,0.5);
  sVect<point3d> points(8);
  
  points(1) = point3d(1.0, 1.0, 1.0);
  points(2) = point3d(2.0, 1.0, 1.0);
  points(3) = point3d(2.0, 2.0, 1.0);
  points(4) = point3d(1.0, 2.0, 1.0);
  
  points(5) = point3d(1.0, 1.0, 2.0);
  points(6) = point3d(2.0, 1.0, 2.0);
  points(7) = point3d(2.0, 2.0, 2.0);
  points(8) = point3d(1.0, 2.0, 2.0);
  
  //Test eval
  cout << "eval:     " << interface.getPosition(points,Y) << endl;
  cout << "gradient: " << interface.getGradient(points,Y) << endl;
  cout << "degree:   " << interface.getDegree() << endl;
}

/*
refShape   : 5
refGeo     : 4
dim        : 3
numFaces   : 6 6
numEdges   : 12 12
numVertices: 8 8
eval:     1.5 1.5 1.5

gradient: Tensor components
1 0 0 
0 1 0 
0 0 1 


degree:   3
*/
