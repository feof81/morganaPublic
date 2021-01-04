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

#include "pVect.hpp"
#include "pVectComm.hpp"

#include "pGraphItem.h"
#include "pGraph.hpp"
#include "pGraphManip.hpp"
#include "pGraphComm.hpp"
#include "pGraphGlobalManip.hpp"

#include "traitsMpiOptimization.hpp"

using namespace Teuchos;


//! Run with two processors
int main(int argc, char *argv[])
{
  typedef pGraphItem ITEM;
  typedef pMapItem   ROWMAP;
  typedef pMapItem   COLMAP;
  
  environment  env(argc,argv);
  RCP<communicator> world(new communicator);
  
  assert(world->size() == 2);
  
  ITEM grphItem(3);
  ROWMAP pItem;
  
  pMap<COLMAP> colMap;
  pGraph<ITEM,ROWMAP,COLMAP>  grafo;
  
  if(world->rank() == 0)
  {
    pItem.setLid(1); pItem.setGid(1); colMap.push_back(pItem);
    pItem.setLid(2); pItem.setGid(2); colMap.push_back(pItem);
    pItem.setLid(3); pItem.setGid(3); colMap.push_back(pItem);
    pItem.setLid(4); pItem.setGid(4); colMap.push_back(pItem);
    pItem.setLid(5); pItem.setGid(5); colMap.push_back(pItem);
    pItem.setLid(6); pItem.setGid(6); colMap.push_back(pItem);
  
    pItem.setLid(1); pItem.setGid(0); grphItem(1) = 1; grphItem(2) = 5; grphItem(3) = 6; grafo.push_back(pItem,grphItem);
    pItem.setLid(2); pItem.setGid(0); grphItem(1) = 1; grphItem(2) = 2; grphItem(3) = 5; grafo.push_back(pItem,grphItem);
    pItem.setLid(3); pItem.setGid(0); grphItem(1) = 2; grphItem(2) = 4; grphItem(3) = 5; grafo.push_back(pItem,grphItem);
    pItem.setLid(4); pItem.setGid(0); grphItem(1) = 2; grphItem(2) = 3; grphItem(3) = 4; grafo.push_back(pItem,grphItem);
  
    grafo.setColMap(colMap);
    grafo.updateRowFinder();
    grafo.updateColFinder();
  }
  else
  {
    pItem.setLid(1); pItem.setGid(2); colMap.push_back(pItem);
    pItem.setLid(2); pItem.setGid(3); colMap.push_back(pItem);
    pItem.setLid(3); pItem.setGid(7); colMap.push_back(pItem);
    pItem.setLid(4); pItem.setGid(8); colMap.push_back(pItem);
    pItem.setLid(5); pItem.setGid(4); colMap.push_back(pItem);
    pItem.setLid(6); pItem.setGid(5); colMap.push_back(pItem);
    
    pItem.setLid(1); pItem.setGid(0); grphItem(1) = 1; grphItem(2) = 5; grphItem(3) = 6; grafo.push_back(pItem,grphItem);
    pItem.setLid(2); pItem.setGid(0); grphItem(1) = 1; grphItem(2) = 2; grphItem(3) = 5; grafo.push_back(pItem,grphItem);
    pItem.setLid(3); pItem.setGid(0); grphItem(1) = 2; grphItem(2) = 4; grphItem(3) = 5; grafo.push_back(pItem,grphItem);
    pItem.setLid(4); pItem.setGid(0); grphItem(1) = 2; grphItem(2) = 3; grphItem(3) = 4; grafo.push_back(pItem,grphItem);
    
    grafo.setColMap(colMap);
    grafo.updateRowFinder();
    grafo.updateColFinder();
  }
  
  grafo.pushToGlobal();
  
  pVectGlobalManip<ITEM,ROWMAP> fixer(world);
  fixer.buildGlobalNumbering(grafo);
  
  
  world->barrier();
  if(world->rank() == 0)
  { cout << grafo << endl;}
  sleep(0.2);
  
  world->barrier();
  if(world->rank() == 1)
  { cout << grafo << endl;}
  sleep(0.2);
}

/* Pid 0
ROW
1
 map:  pid: 0 lid: 1 gid: 2
 data: Num Connected  : 3
Connected Id's : 1 5 6 

2
 map:  pid: 0 lid: 2 gid: 1
 data: Num Connected  : 3
Connected Id's : 1 2 5 

3
 map:  pid: 0 lid: 3 gid: 4
 data: Num Connected  : 3
Connected Id's : 2 4 5 

4
 map:  pid: 0 lid: 4 gid: 3
 data: Num Connected  : 3
Connected Id's : 2 3 4 

COL
pid: 0 lid: 1 gid: 1
pid: 0 lid: 2 gid: 2
pid: 0 lid: 3 gid: 3
pid: 0 lid: 4 gid: 4
pid: 0 lid: 5 gid: 5
pid: 0 lid: 6 gid: 6
*/

/* Pid1
ROW
1
 map:  pid: 1 lid: 1 gid: 4
 data: Num Connected  : 3
Connected Id's : 2 4 5 

2
 map:  pid: 1 lid: 2 gid: 3
 data: Num Connected  : 3
Connected Id's : 2 3 4 

3
 map:  pid: 1 lid: 3 gid: 5
 data: Num Connected  : 3
Connected Id's : 3 8 4 

4
 map:  pid: 1 lid: 4 gid: 6
 data: Num Connected  : 3
Connected Id's : 3 7 8 

COL
pid: 0 lid: 1 gid: 2
pid: 0 lid: 2 gid: 3
pid: 0 lid: 3 gid: 7
pid: 0 lid: 4 gid: 8
pid: 0 lid: 5 gid: 4
pid: 0 lid: 6 gid: 5
*/
