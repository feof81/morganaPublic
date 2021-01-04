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


//! Test with graph items
int main(int argc, char *argv[])
{
  typedef pGraphItem ITEM;
  typedef pMapItem   ROWMAP;
  typedef pMapItem   COLMAP;
  typedef pVect<ITEM,ROWMAP> PVECT;
  typedef pGraph<ITEM,ROWMAP,COLMAP> PGRAPH;
  
  ITEM grphItem(3);
  ROWMAP pItem;
  
  pMap<COLMAP> colMap;
  PGRAPH       grafo;
  
  pItem.setLid(1); pItem.setGid(1); colMap.push_back(pItem);
  pItem.setLid(2); pItem.setGid(2); colMap.push_back(pItem);
  pItem.setLid(3); pItem.setGid(3); colMap.push_back(pItem);
  pItem.setLid(4); pItem.setGid(4); colMap.push_back(pItem);
  pItem.setLid(5); pItem.setGid(5); colMap.push_back(pItem);
  pItem.setLid(6); pItem.setGid(6); colMap.push_back(pItem);
  
  pItem.setLid(1); pItem.setGid(1); grphItem(1) = 1; grphItem(2) = 5; grphItem(3) = 6; grafo.push_back(pItem,grphItem);
  pItem.setLid(2); pItem.setGid(2); grphItem(1) = 1; grphItem(2) = 2; grphItem(3) = 5; grafo.push_back(pItem,grphItem);
  pItem.setLid(3); pItem.setGid(3); grphItem(1) = 2; grphItem(2) = 4; grphItem(3) = 5; grafo.push_back(pItem,grphItem);
  pItem.setLid(4); pItem.setGid(4); grphItem(1) = 2; grphItem(2) = 3; grphItem(3) = 4; grafo.push_back(pItem,grphItem);
  
  grafo.setColMap(colMap);
  grafo.updateRowFinder();
  grafo.updateColFinder();
  
  //Max-Min determination
  set<ITEM> lista;
  
  for(UInt i=1; i <= grafo.rowSize(); ++i)
  {
    lista.insert(grafo.getItemL(i));
  }
  
  ITEM minItem = *lista.begin();
  ITEM maxItem;
 
  set<ITEM>::iterator iter;
  for(iter = lista.begin(); iter != lista.end(); ++iter)
  { maxItem = *iter; }
  
  assert(minItem < maxItem);
  
  cout << "min: " << minItem << endl;
  cout << "max: " << maxItem << endl;
  
  
  //Ordering test
  ITEM itemA(3), itemB(3);
  
  itemA(1) = 1; itemA(2) = 2; itemA(3) = 5;
  itemB(1) = 1; itemB(2) = 3; itemB(3) = 2;
  
  itemA.updateSorting();
  itemB.updateSorting();
  
  cout << "itemA" << endl;
  itemA.printSorted();
  
  cout << "itemB" << endl;
  itemB.printSorted();
  
  cout << "test: " << (itemB < itemA) << endl << endl;
  
  
  //Subclass testing
  UInt segments = 4;
  
  dataSegmentationUtility<ITEM> utility(segments,minItem,maxItem);
  
  cout << "Segments" << endl;
  cout << utility << endl;
  
  
  //Segmentation
  sVect<PVECT> targetVects;
  
  pVectManip<ITEM,ROWMAP> segmenter(grafo);
  segmenter.segmentationData(targetVects,segments,minItem,maxItem);
  
  
  for(UInt i=1; i <= targetVects.size(); ++i)
  {
    cout << "SEGMENT: " << i << endl;
    cout << targetVects(i) << endl << endl;
  }
}

/*
min: Num Connected  : 3
Connected Id's : 1 2 5 

max: Num Connected  : 3
Connected Id's : 2 4 5 

itemA
ordered ids
1 2 5 
itemB
ordered ids
1 2 3 
test: 1

Segments
Num Connected  : 3
Connected Id's : 1 2 5 

Num Connected  : 3
Connected Id's : 1 3 5 

Num Connected  : 3
Connected Id's : 2 4 5 


SEGMENT: 1


SEGMENT: 2
1
 map:  pid: 0 lid: 2 gid: 2
 data: Num Connected  : 3
Connected Id's : 1 2 5 



SEGMENT: 3


SEGMENT: 4
1
 map:  pid: 0 lid: 1 gid: 1
 data: Num Connected  : 3
Connected Id's : 1 5 6 

2
 map:  pid: 0 lid: 3 gid: 3
 data: Num Connected  : 3
Connected Id's : 2 4 5 

3
 map:  pid: 0 lid: 4 gid: 4
 data: Num Connected  : 3
Connected Id's : 2 3 4
*/
