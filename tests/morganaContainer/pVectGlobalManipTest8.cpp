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
#include "pMapItemShare.h"
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
  
  typedef pVect<point3d,pMapItemShare> PVECT;
  PVECT sVector, rVector;    
  pMapItemShare item;
  
  
  if(world->rank() == 0)
  {
    item.setPid(0); item.setLid(1); item.setGid(5); sVector.push_back(item,point3d( 1.0,-2.0, 0.0),false);
    item.setPid(0); item.setLid(2); item.setGid(6); sVector.push_back(item,point3d(-1.0,-1.0,-1.0),false);
    
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
  
  point3d pointMax, pointMin;
  
  pVectGlobalManip<point3d,pMapItemShare> vectorManip(world);
  vectorManip.buildGlobalNumbering(sVector);
  vectorManip.dataMinMax(sVector,pointMin,pointMax);
  
  
  world->barrier();
  if(world->rank() == 0)
  {
    cout << "Min: " << endl;
    cout << pointMin << endl << endl;
    
    cout << "Max: " << endl;
    cout << pointMax << endl;
  }
  
  world->barrier();
  if(world->rank() == 1)
  {
    cout << "Min: " << endl;
    cout << pointMin << endl << endl;
    
    cout << "Max: " << endl;
    cout << pointMax << endl;
  }
  
  
  pVectComm<point3d,pMapItemShare> vectCommunicator(world);
  vectCommunicator.vectorData(sVector,pointMin,pointMax);
  
  
  world->barrier();
  if(world->rank() == 0)
  { cout << sVector << endl; }
  sleep(1);
  
  world->barrier();
  if(world->rank() == 1)
  { cout << sVector << endl; }
  sleep(1);
}


/*
Min: 
-1 -1 -1


Max: 
1 -2 0

1
 map:  pid: 1 lid: 1 gid: 2 shared: 1 owned 1
 data: -0.5 1 0

2
 map:  pid: 1 lid: 2 gid: 2 shared: 1 owned 0
 data: -0.5 1 0

3
 map:  pid: 1 lid: 3 gid: 1 shared: 1 owned 1
 data: -1 -1 -1

4
 map:  pid: 0 lid: 4 gid: 1 shared: 1 owned 0
 data: -1 -1 -1
*/


/*
Min: 
-1 -1 -1


Max: 
1 -2 0

1
 map:  pid: 0 lid: 1 gid: 4 shared: 0 owned 1
 data: 1 -2 0

2
 map:  pid: 1 lid: 2 gid: 3 shared: 0 owned 1
 data: 0 0 1
*/
