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

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "pMap.hpp"
#include "pMapItemShare.h"
#include "pMpiOptimization.hpp"
#include "pMapComm.hpp"
#include "pMapGlobalManip.h"

using namespace Teuchos;


//! Run with three processors
int main(int argc, char *argv[])
{
  //Comunication stuff
  environment  env(argc,argv);
  RCP<communicator> world(new communicator);
  
  assert(world->size() == 3);
  
  
  //Common data
  typedef pMap<pMapItemShare> OBJ;
  
  pMapItemShare  item;
  RCP<OBJ>  sMap(new OBJ);
  RCP<OBJ>  rMap(new OBJ);
  
  if(world->rank() == 1)
  {  
    item.setLid(1); item.setGid(1);  item.setOwned(true); item.setShared(false);  sMap->push_back(item);
    item.setLid(2); item.setGid(2);  item.setOwned(true); item.setShared(false);  sMap->push_back(item);
    item.setLid(3); item.setGid(3);  item.setOwned(true); item.setShared(false);  sMap->push_back(item);
    item.setLid(4); item.setGid(4);  item.setOwned(true); item.setShared(true);   sMap->push_back(item);
    item.setLid(5); item.setGid(5);  item.setOwned(true); item.setShared(true);   sMap->push_back(item);
  }
  
  if(world->rank() == 2)
  {
    item.setLid(1); item.setGid(4);  item.setOwned(false); item.setShared(true);  sMap->push_back(item);
    item.setLid(2); item.setGid(5);  item.setOwned(false); item.setShared(true);  sMap->push_back(item);
    item.setLid(3); item.setGid(6);  item.setOwned(true);  item.setShared(false); sMap->push_back(item);
    item.setLid(4); item.setGid(7);  item.setOwned(true);  item.setShared(false); sMap->push_back(item);
    item.setLid(5); item.setGid(8);  item.setOwned(true);  item.setShared(false); sMap->push_back(item);
    item.setLid(6); item.setGid(9);  item.setOwned(true);  item.setShared(true);  sMap->push_back(item);
    item.setLid(7); item.setGid(10); item.setOwned(true);  item.setShared(true);  sMap->push_back(item);
  }
  
  if(world->rank() == 0)
  {
    item.setLid(1); item.setGid(9);  item.setOwned(false); item.setShared(true);  sMap->push_back(item);
    item.setLid(2); item.setGid(10); item.setOwned(false); item.setShared(true);  sMap->push_back(item);
    item.setLid(3); item.setGid(11); item.setOwned(true);  item.setShared(false); sMap->push_back(item);
    item.setLid(4); item.setGid(12); item.setOwned(true);  item.setShared(false); sMap->push_back(item);
    item.setLid(5); item.setGid(13); item.setOwned(true);  item.setShared(false); sMap->push_back(item);
    item.setLid(6); item.setGid(14); item.setOwned(true);  item.setShared(false); sMap->push_back(item);
    item.setLid(7); item.setGid(15); item.setOwned(true);  item.setShared(false); sMap->push_back(item);
  }
  
  
  pMapGlobalManip<pMapItemShare> manipulator(world);
  
  UInt sizeG   = manipulator.sizeG(sMap);
  UInt sharedL = manipulator.sharedL(sMap);
  UInt sharedG = manipulator.sharedG(sMap);
  UInt ownedL  = manipulator.ownedL(sMap);
  
  world->barrier();
  
  if(world->rank() == 0)
  {
    cout << "SizeG:     " << sizeG << endl;
    cout << "SharedL:   " << sharedL << endl;
    cout << "SharedG:   " << sharedG << endl;
    cout << "OwnedL:    " << ownedL << endl;
  }
  
  world->barrier();
  if(world->rank() == 1)
  {
    cout << "SizeG:     " << sizeG << endl;
    cout << "SharedL:   " << sharedL << endl;
    cout << "SharedG:   " << sharedG << endl;
    cout << "OwnedL:    " << ownedL << endl;
  }
  sleep(1);
  
  world->barrier();
  if(world->rank() == 2)
  {
    cout << "SizeG:     " << sizeG << endl;
    cout << "SharedL:   " << sharedL << endl;
    cout << "SharedG:   " << sharedG << endl;
    cout << "OwnedL:    " << ownedL << endl;
  }
  sleep(1);
}

/* pid 0
SizeG:     15
SharedL:   2
SharedG:   4
OwnedL:    5
*/


/* pid 1
SizeG:     15
SharedL:   2
SharedG:   4
OwnedL:    5
*/


/* pid 2
SizeG:     15
SharedL:   4
SharedG:   4
OwnedL:    5
*/
