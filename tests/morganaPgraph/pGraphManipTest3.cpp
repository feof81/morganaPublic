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
#include "pGraphManip.hpp"
#include "traitsMpiOptimization.hpp"


using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


/*! Run with one processor */
int main(int argc, char *argv[])
{
  typedef pGraphItem ITEM;
  typedef pMapItem   ROWMAP;
  typedef pMapItem   COLMAP;
  
  pMapItem   pItem;
  pGraphItem grphItem(3);
  
  pMap<pMapItem> colMap;
  pGraph<ITEM,ROWMAP,COLMAP> grafo;
  
  pItem.setLid(4); pItem.setGid(24); colMap.push_back(pItem);
  pItem.setLid(3); pItem.setGid(23); colMap.push_back(pItem);
  pItem.setLid(2); pItem.setGid(22); colMap.push_back(pItem);
  pItem.setLid(1); pItem.setGid(21); colMap.push_back(pItem);
  
  pItem.setLid(2); pItem.setGid(12); grphItem(1) = 2; grphItem(2) = 3; grphItem(3) = 4; grafo.push_back(pItem,grphItem);
  pItem.setLid(1); pItem.setGid(11); grphItem(1) = 1; grphItem(2) = 2; grphItem(3) = 3; grafo.push_back(pItem,grphItem);
  
  grafo.setColMap(colMap);
  
  //Finder test
  pGraphManip<ITEM,ROWMAP,COLMAP> manipulator(grafo);
  manipulator.buildFinder();
  
  manipulator.setIndexing();
  
  cout << "max Row: " << manipulator.getMaxRowGid() << endl;
  cout << "max Col: " << manipulator.getMaxColGid() << endl;
  cout << "check: " << manipulator.checkCol() << endl;
}

/*
max Row: 12
max Col: 24
check: 1
*/
