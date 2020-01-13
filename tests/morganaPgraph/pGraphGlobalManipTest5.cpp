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
  typedef pGraphItem     ITEM;
  typedef pMapItem       ROWMAP;
  typedef pMapItemShare  COLMAP;
  
  environment  env(argc,argv);
  RCP<communicator> world(new communicator);
  
  assert(world->size() == 2);
  
  ITEM grphItem(3);
  ROWMAP pItem;
  COLMAP pShare;
  
  pMap<COLMAP> colMap;
  pGraph<ITEM,ROWMAP,COLMAP>  grafo;
  pGraphGlobalManip<ITEM,ROWMAP,COLMAP> fixer(world);
  
  if(world->rank() == 0)
  {
    pItem.setPid(world->rank());
    pShare.setPid(world->rank());
    
    pShare.setLid(1); pShare.setGid(1); pShare.setShared(false); pShare.setOwned(true);  colMap.push_back(pShare);
    pShare.setLid(2); pShare.setGid(2); pShare.setShared(true);  pShare.setOwned(true);  colMap.push_back(pShare);
    pShare.setLid(3); pShare.setGid(3); pShare.setShared(true);  pShare.setOwned(false); colMap.push_back(pShare);
    pShare.setLid(4); pShare.setGid(4); pShare.setShared(true);  pShare.setOwned(false); colMap.push_back(pShare);
    pShare.setLid(5); pShare.setGid(5); pShare.setShared(true);  pShare.setOwned(true);  colMap.push_back(pShare);
    pShare.setLid(6); pShare.setGid(6); pShare.setShared(false); pShare.setOwned(true);  colMap.push_back(pShare);
  
    pItem.setLid(1); pItem.setGid(2); grphItem(1) = 1; grphItem(2) = 5; grphItem(3) = 6; grafo.push_back(pItem,grphItem);
    pItem.setLid(2); pItem.setGid(1); grphItem(1) = 1; grphItem(2) = 2; grphItem(3) = 5; grafo.push_back(pItem,grphItem);
    pItem.setLid(3); pItem.setGid(4); grphItem(1) = 2; grphItem(2) = 4; grphItem(3) = 5; grafo.push_back(pItem,grphItem);
    pItem.setLid(4); pItem.setGid(3); grphItem(1) = 2; grphItem(2) = 3; grphItem(3) = 4; grafo.push_back(pItem,grphItem);
  
    grafo.setColMap(colMap);
    grafo.updateRowFinder();
    grafo.updateColFinder();
  }
  else
  {
    pItem.setPid(world->rank());
    pShare.setPid(world->rank());
    
    pShare.setLid(1); pShare.setGid(2); pShare.setShared(true);  pShare.setOwned(false); colMap.push_back(pShare);
    pShare.setLid(2); pShare.setGid(3); pShare.setShared(true);  pShare.setOwned(true);  colMap.push_back(pShare);
    pShare.setLid(3); pShare.setGid(7); pShare.setShared(false); pShare.setOwned(true);  colMap.push_back(pShare);
    pShare.setLid(4); pShare.setGid(8); pShare.setShared(false); pShare.setOwned(true);  colMap.push_back(pShare);
    pShare.setLid(5); pShare.setGid(4); pShare.setShared(true);  pShare.setOwned(true);  colMap.push_back(pShare);
    pShare.setLid(6); pShare.setGid(5); pShare.setShared(true);  pShare.setOwned(false); colMap.push_back(pShare);
    
    pItem.setLid(1); pItem.setGid(4); grphItem(1) = 1; grphItem(2) = 5; grphItem(3) = 6; grafo.push_back(pItem,grphItem);
    pItem.setLid(2); pItem.setGid(3); grphItem(1) = 1; grphItem(2) = 2; grphItem(3) = 5; grafo.push_back(pItem,grphItem);
    pItem.setLid(3); pItem.setGid(5); grphItem(1) = 2; grphItem(2) = 4; grphItem(3) = 5; grafo.push_back(pItem,grphItem);
    pItem.setLid(4); pItem.setGid(6); grphItem(1) = 2; grphItem(2) = 3; grphItem(3) = 4; grafo.push_back(pItem,grphItem);
    
    grafo.setColMap(colMap);
    grafo.updateRowFinder();
    grafo.updateColFinder();
  }
  
  
  fixer.overlap(grafo);
  
  
  world->barrier();
  if(world->rank() == 0)
  { cout << grafo << endl; }  
  
  world->barrier();
  if(world->rank() == 1)
  { cout << grafo << endl; }
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

5
 map:  pid: 0 lid: 5 gid: 5
 data: Num Connected  : 3
Connected Id's : 3 7 4 

6
 map:  pid: 0 lid: 6 gid: 6
 data: Num Connected  : 3
Connected Id's : 3 8 7 

COL
pid: 0 lid: 1 gid: 1 shared: 1 owned 1
pid: 0 lid: 2 gid: 2 shared: 1 owned 1
pid: 0 lid: 3 gid: 3 shared: 1 owned 0
pid: 0 lid: 4 gid: 4 shared: 1 owned 0
pid: 0 lid: 5 gid: 5 shared: 1 owned 1
pid: 0 lid: 6 gid: 6 shared: 1 owned 1
pid: 0 lid: 7 gid: 8 shared: 1 owned 0
pid: 0 lid: 8 gid: 7 shared: 1 owned 0
*/


/* Pid 1
ROW
1
 map:  pid: 1 lid: 1 gid: 4
 data: Num Connected  : 3
Connected Id's : 1 5 6 

2
 map:  pid: 1 lid: 2 gid: 3
 data: Num Connected  : 3
Connected Id's : 1 2 5 

3
 map:  pid: 1 lid: 3 gid: 5
 data: Num Connected  : 3
Connected Id's : 2 4 5 

4
 map:  pid: 1 lid: 4 gid: 6
 data: Num Connected  : 3
Connected Id's : 2 3 4 

5
 map:  pid: 1 lid: 5 gid: 1
 data: Num Connected  : 3
Connected Id's : 7 1 6 

6
 map:  pid: 1 lid: 6 gid: 2
 data: Num Connected  : 3
Connected Id's : 7 6 8 

COL
pid: 1 lid: 1 gid: 2 shared: 1 owned 0
pid: 1 lid: 2 gid: 3 shared: 1 owned 1
pid: 1 lid: 3 gid: 7 shared: 1 owned 1
pid: 1 lid: 4 gid: 8 shared: 1 owned 1
pid: 1 lid: 5 gid: 4 shared: 1 owned 1
pid: 1 lid: 6 gid: 5 shared: 1 owned 0
pid: 1 lid: 7 gid: 1 shared: 1 owned 0
pid: 1 lid: 8 gid: 6 shared: 1 owned 0
*/
