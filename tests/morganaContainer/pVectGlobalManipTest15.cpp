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
#include "pMapItem.h"
#include "pVect.hpp"
#include "pVectManip.hpp"
#include "pVectComm.hpp"
#include "pVectGlobalManip.hpp"

using namespace Teuchos;


//! Run with two processors
int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  RCP<communicator> world(new communicator);
  
  assert(world->size() == 2);
  
  typedef pVect<point3d,pMapItem> PVECT;
  PVECT sVector, rVector;    
  pMapItem item;
  
  pVectGlobalManip<point3d,pMapItem> vectorManip(world);
  
  if(world->rank() == 0)
  {
    item.setPid(0); item.setLid(1); item.setGid(5); sVector.push_back(item,point3d( 1.0,-2.0, 0.0),false);
    item.setPid(0); item.setLid(2); item.setGid(6); sVector.push_back(item,point3d(-1.0,-1.0,-1.0),false);
    item.setPid(0); item.setLid(3); item.setGid(1); sVector.push_back(item,point3d( 1.0, 2.0,-1.0),false);
    item.setPid(0); item.setLid(4); item.setGid(2); sVector.push_back(item,point3d( 2.0, 3.0,-2.0),false);
    
    sVector.updateFinder();
  }
  
  if(world->rank() == 1)
  {
    item.setPid(1); item.setLid(1); item.setGid(1); sVector.push_back(item,point3d(-0.5, 1.0, 0.0),false);
    item.setPid(1); item.setLid(2); item.setGid(2); sVector.push_back(item,point3d(-0.5, 1.0, 0.0),false);
    item.setPid(1); item.setLid(3); item.setGid(3); sVector.push_back(item,point3d( 0.0, 0.0, 1.0),false);
    item.setPid(1); item.setLid(4); item.setGid(4); sVector.push_back(item,point3d(-1.0,-1.0,-1.0),false);
    
    sVector.updateFinder();
  }
  
  std::plus<Real> op;
  vectorManip.allReduce(sVector,op);
  
  world->barrier();
  if(world->rank() == 0)
  { cout << sVector << endl; }
  sleep(0.2);
  
  world->barrier();
  if(world->rank() == 1)
  { cout << sVector << endl; }
  sleep(0.2);
}

/*
1
 map:  pid: 0 lid: 1 gid: 5
 data: 1 -2 0

2
 map:  pid: 0 lid: 2 gid: 6
 data: -1 -1 -1

3
 map:  pid: 0 lid: 3 gid: 1
 data: 0.5 3 -1

4
 map:  pid: 0 lid: 4 gid: 2
 data: 1.5 4 -2


1
 map:  pid: 1 lid: 1 gid: 1
 data: 0.5 3 -1

2
 map:  pid: 1 lid: 2 gid: 2
 data: 1.5 4 -2

3
 map:  pid: 1 lid: 3 gid: 3
 data: 0 0 1

4
 map:  pid: 1 lid: 4 gid: 4
 data: -1 -1 -1

*/

