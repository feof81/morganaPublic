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
#include "pVect.hpp"
#include "pMapItem.h"


/*! Run with one processor */
int main(int argc, char *argv[])
{
  //Data
  sVect<point3d> data;
  data.push_back(point3d(1.0,0.0,0.0));
  data.push_back(point3d(0.0,1.0,0.0));
  data.push_back(point3d(0.0,0.0,1.0));
  
  //Map
  pMap<pMapItem> map;
  map.push_back(pMapItem(0,3));
  map.push_back(pMapItem(0,5));
  map.push_back(pMapItem(0,7));
  
  //pVect
  pVect<point3d,pMapItem> trunk(map,data);
  cout << trunk << endl;
  cout << "size: " << trunk.sizeL() << endl;
  
  point3d P(1.0,1.0,1.0);
  trunk.getG(5) = P;
  cout << trunk.getL(2) << endl;
}

/*
1
 map:  pid: 0 lid: 1 gid: 3
 data: 1 0 0

2
 map:  pid: 0 lid: 2 gid: 5
 data: 0 1 0

3
 map:  pid: 0 lid: 3 gid: 7
 data: 0 0 1


size: 3
1 1 1
*/
