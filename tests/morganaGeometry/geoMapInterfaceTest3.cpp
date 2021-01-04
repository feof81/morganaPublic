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
  geoMapInterface<linearTriangle> interface;
  
  cout << "refShape   : " << interface.getReferenceShapes() << endl;
  cout << "refGeo     : " << interface.getReferenceGeometry() << endl;
  cout << "dim        : " << interface.getNDim() << endl;
  cout << "numFaces   : " << interface.getNumFaces()    << " " << interface.getNumGeoTypes(FACE) << endl;
  cout << "numEdges   : " << interface.getNumEdges()    << " " << interface.getNumGeoTypes(EDGE) << endl;
  cout << "numVertices: " << interface.getNumVertices() << " " << interface.getNumGeoTypes(VERTEX) << endl;

  
  point3d Y(0.0,0.5,0.0);
  sVect<point3d> points(3);
  
  points(1) = point3d(1.0,1.0,0.0);
  points(2) = point3d(1.0,2.0,0.0);
  points(3) = point3d(2.0,1.0,0.0);
  
  //Test eval
  cout << "eval:     " << interface.getPosition(points,Y) << endl;
  cout << "gradient: " << interface.getGradient(points,Y) << endl;
  cout << "degree:   " << interface.getDegree() << endl;
}

/*
refShape   : 3
refGeo     : 3
dim        : 2
numFaces   : 1 1
numEdges   : 3 3
numVertices: 3 3
eval:     1.5 1 0

gradient: Tensor components
0 1 0 
1 0 0 
0 0 0 


degree:   1 
 */
