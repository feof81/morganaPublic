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
#include "pMapItem.h"
#include "pMpiOptimization.hpp"
#include "pMapComm.hpp"

using namespace Teuchos;


/*! Run with two processors */
int main(int argc, char *argv[])
{
  //Comunication stuff
  environment  env(argc,argv);
  RCP<communicator> world(new communicator);
  
  assert(world->size() == 2);
  
  //Common data
  typedef pMap<pMapItem> OBJ;
  pMapItem  item;
  RCP<OBJ>  sMap(new OBJ);
    
  
  if(world->rank() == 0)
  {
    item.setLid(1); item.setGid(7);  sMap->push_back(item);
    item.setLid(2); item.setGid(8);  sMap->push_back(item);
    item.setLid(3); item.setGid(9);  sMap->push_back(item);
    item.setLid(4); item.setGid(10); sMap->push_back(item);
    item.setLid(5); item.setGid(11); sMap->push_back(item);
    item.setLid(6); item.setGid(12); sMap->push_back(item);
    
    pMapComm<pMapItem> mapComm(world);
    mapComm.vectorNormal(sMap);
  }
  
  if(world->rank() == 1)
  {
    item.setLid(1); item.setGid(1); sMap->push_back(item);
    item.setLid(2); item.setGid(2); sMap->push_back(item);
    item.setLid(3); item.setGid(3); sMap->push_back(item);
    item.setLid(4); item.setGid(4); sMap->push_back(item);
    item.setLid(5); item.setGid(5); sMap->push_back(item);
    item.setLid(6); item.setGid(6); sMap->push_back(item);
    
    pMapComm<pMapItem> mapComm(world);
    mapComm.vectorNormal(sMap);
  }

  world->barrier();
  
  if(world->rank() == 0)
  {
    cout << "Processor: " << world->rank() << endl;
    cout << *sMap << endl;
  }
  
  world->barrier();
  
  if(world->rank() == 1)
  {
    cout << "Processor: " << world->rank() << endl;
    cout << *sMap << endl;
  }
}

/*
Processor: 0
pid: 0 lid: 1 gid: 1
pid: 0 lid: 2 gid: 2
pid: 0 lid: 3 gid: 3
pid: 0 lid: 4 gid: 4
pid: 0 lid: 5 gid: 5
pid: 0 lid: 6 gid: 6
*/

/*
Processor: 1
pid: 0 lid: 1 gid: 7
pid: 0 lid: 2 gid: 8
pid: 0 lid: 3 gid: 9
pid: 0 lid: 4 gid: 10
pid: 0 lid: 5 gid: 11
pid: 0 lid: 6 gid: 12
*/
