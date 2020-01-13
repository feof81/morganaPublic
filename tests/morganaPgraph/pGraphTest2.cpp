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

#include "pMapItem.h"
#include "pMap.hpp"

#include "pGraphItem.h"
#include "pGraph.hpp"
#include "traitsMpiOptimization.hpp"

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;

// mpirun -np 2 ./bin/morgana

/*! Run with two processors */
int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  
  pMapItem   pItem;
  pGraphItem grphItem(3);
  
  pMap<pMapItem> colMap;
  pGraph<pGraphItem,pMapItem,pMapItem> sendGraph, recvGraph;
  
  if(world.rank() == 0)
  {
    pItem.setLid(1); pItem.setGid(21); colMap.push_back(pItem);
    pItem.setLid(2); pItem.setGid(22); colMap.push_back(pItem);
    pItem.setLid(3); pItem.setGid(23); colMap.push_back(pItem);
    pItem.setLid(4); pItem.setGid(24); colMap.push_back(pItem);
  
    pItem.setLid(1); pItem.setGid(11); grphItem(1) = 1; grphItem(2) = 2; grphItem(3) = 3; sendGraph.push_back(pItem,grphItem);
    pItem.setLid(2); pItem.setGid(12); grphItem(1) = 2; grphItem(2) = 3; grphItem(3) = 4; sendGraph.push_back(pItem,grphItem);
  
    sendGraph.updateFinder();
    sendGraph.setColMap(colMap);
    
    world.send(1,0,sendGraph);
    world.recv(1,1,recvGraph);
  }
  else
  {
    pItem.setLid(1); pItem.setGid(121); colMap.push_back(pItem);
    pItem.setLid(2); pItem.setGid(122); colMap.push_back(pItem);
    pItem.setLid(3); pItem.setGid(123); colMap.push_back(pItem);
    pItem.setLid(4); pItem.setGid(124); colMap.push_back(pItem);
  
    pItem.setLid(1); pItem.setGid(111); grphItem(1) = 1; grphItem(2) = 2; grphItem(3) = 3; sendGraph.push_back(pItem,grphItem);
    pItem.setLid(2); pItem.setGid(112); grphItem(1) = 2; grphItem(2) = 3; grphItem(3) = 4; sendGraph.push_back(pItem,grphItem);
  
    sendGraph.updateFinder();
    sendGraph.setColMap(colMap);
    
    world.recv(0,0,recvGraph);
    world.send(0,1,sendGraph);
  }
  
  
  world.barrier();
  if(world.rank() == 0)
  { cout << recvGraph << endl; }
  sleep(0.2);
  
  world.barrier();
  if(world.rank() == 1)
  { cout << recvGraph << endl; }
  sleep(0.2);
}

/* Pid 0
ROW
1
 map:  pid: 0 lid: 1 gid: 111
 data: Num Connected  : 3
Connected Id's : 1 2 3 

2
 map:  pid: 0 lid: 2 gid: 112
 data: Num Connected  : 3
Connected Id's : 2 3 4 

COL
pid: 0 lid: 1 gid: 121
pid: 0 lid: 2 gid: 122
pid: 0 lid: 3 gid: 123
pid: 0 lid: 4 gid: 124
*/

/* Pid 1
ROW
1
 map:  pid: 0 lid: 1 gid: 11
 data: Num Connected  : 3
Connected Id's : 1 2 3 

2
 map:  pid: 0 lid: 2 gid: 12
 data: Num Connected  : 3
Connected Id's : 2 3 4 

COL
pid: 0 lid: 1 gid: 21
pid: 0 lid: 2 gid: 22
pid: 0 lid: 3 gid: 23
pid: 0 lid: 4 gid: 24
*/
