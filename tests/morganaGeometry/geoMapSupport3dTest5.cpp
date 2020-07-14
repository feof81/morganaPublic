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

#include "geoMapSupport3d.hpp"

int main(int argc, char *argv[])
{
  //Geo support tor tretra
  geoMapSupport3d<linearTetra> supporter;
  
  point3d P, Y;
  sVect<point3d> points(4);
  
  points(1) = point3d(0.0,0.0,0.0);
  points(2) = point3d(1.0,0.0,0.0);
  points(3) = point3d(1.0,1.0,0.0);
  points(4) = point3d(0.5,0.5,1.0);
  
  supporter.setPoints(points);
  
  P.set(0.5, 0.5, -1.0);
  cout << "found: " << supporter.projection(P,Y) << endl;
  cout << "Y    : " << Y;
  cout << "Q    : " << supporter.getPosition(points,Y) << endl;
  
  P.set(0.5, -1.0, 0.0);
  cout << "found: " << supporter.projection(P,Y) << endl;
  cout << "Y    : " << Y;
  cout << "Q    : " << supporter.getPosition(points,Y) << endl;
  
  P.set(0.5, 0.5, 2.0);
  cout << "found: " << supporter.projection(P,Y) << endl;
  cout << "Y    : " << Y;
  cout << "Q    : " << supporter.getPosition(points,Y) << endl;
}

/*
 found: 1
Y    : 0 0.5 0
Q    : 0.5 0.5 0

found: 1
Y    : 0.5 0 0
Q    : 0.5 0 0

found: 1
Y    : 0 0 1
Q    : 0.5 0.5 1
 */
