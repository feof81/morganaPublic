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
  //Typedefs---------------------------------------------------------
  typedef pGraphItem ITEM;
  typedef pMapItem   ROWMAP;
  typedef pMapItem   COLMAP;
  
  //Comm dev---------------------------------------------------------
  boost::mpi::environment  env(argc, argv);
  boost::mpi::communicator newCommDev;
  
  bool commActive = (newCommDev.rank() < 2);
  boost::mpi::communicator oldCommDev = newCommDev.split(commActive);
  
  assert(newCommDev.size() == 4);
  
  //Alloc------------------------------------------------------------
  ITEM grphItem(3);
  ROWMAP pItem;
  
  pGraph<ITEM,ROWMAP,COLMAP> newGraph;
  pGraphGlobalManip<ITEM,ROWMAP,COLMAP> fixer;
  
  //Testing----------------------------------------------------------
  if(commActive)
  {
    pMap<COLMAP> colMap;    
    pGraph<ITEM,ROWMAP,COLMAP> oldGraph;
    
    if(oldCommDev.rank() == 0)
    {
      pItem.setPid(oldCommDev.rank());
      
      pItem.setLid(1); pItem.setGid(1); pItem.setGid(1); colMap.push_back(pItem);
      pItem.setLid(2); pItem.setGid(2); colMap.push_back(pItem);
      pItem.setLid(3); pItem.setGid(3); colMap.push_back(pItem);
      pItem.setLid(4); pItem.setGid(4); colMap.push_back(pItem);
      pItem.setLid(5); pItem.setGid(5); colMap.push_back(pItem);
      pItem.setLid(6); pItem.setGid(6); colMap.push_back(pItem);
    
      pItem.setLid(1); pItem.setGid(1); grphItem(1) = 1; grphItem(2) = 5; grphItem(3) = 6; oldGraph.push_back(pItem,grphItem);
      pItem.setLid(2); pItem.setGid(2); grphItem(1) = 1; grphItem(2) = 2; grphItem(3) = 5; oldGraph.push_back(pItem,grphItem);
      pItem.setLid(3); pItem.setGid(3); grphItem(1) = 2; grphItem(2) = 4; grphItem(3) = 5; oldGraph.push_back(pItem,grphItem);
      pItem.setLid(4); pItem.setGid(4); grphItem(1) = 2; grphItem(2) = 3; grphItem(3) = 4; oldGraph.push_back(pItem,grphItem);
    
      oldGraph.setColMap(colMap);
      oldGraph.updateRowFinder();
      oldGraph.updateColFinder();
    }
    
    if(oldCommDev.rank() == 1)
    {
      pItem.setPid(oldCommDev.rank());
      
      pItem.setLid(1); pItem.setGid(2); colMap.push_back(pItem);
      pItem.setLid(2); pItem.setGid(3); colMap.push_back(pItem);
      pItem.setLid(3); pItem.setGid(7); colMap.push_back(pItem);
      pItem.setLid(4); pItem.setGid(8); colMap.push_back(pItem);
      pItem.setLid(5); pItem.setGid(4); colMap.push_back(pItem);
      pItem.setLid(6); pItem.setGid(5); colMap.push_back(pItem);
      
      pItem.setLid(1); pItem.setGid(3); grphItem(1) = 1; grphItem(2) = 5; grphItem(3) = 6; oldGraph.push_back(pItem,grphItem);
      pItem.setLid(2); pItem.setGid(4); grphItem(1) = 1; grphItem(2) = 2; grphItem(3) = 5; oldGraph.push_back(pItem,grphItem);
      pItem.setLid(3); pItem.setGid(5); grphItem(1) = 2; grphItem(2) = 4; grphItem(3) = 5; oldGraph.push_back(pItem,grphItem);
      pItem.setLid(4); pItem.setGid(6); grphItem(1) = 2; grphItem(2) = 3; grphItem(3) = 4; oldGraph.push_back(pItem,grphItem);
      
      oldGraph.setColMap(colMap);
      oldGraph.updateRowFinder();
      oldGraph.updateColFinder();
    }
    
    fixer.expandCommunicator(commActive,
                             oldCommDev,
                             oldGraph,
                             newCommDev,
                             newGraph);
  }
  
  //Printout---------------------------------------------------------
  for(UInt pid=0; pid < newCommDev.size(); pid++)
  {
    if(pid == newCommDev.rank())
    {
      cout << "Process : " << pid << endl;
      cout << newGraph << endl;
    }
    
    newCommDev.barrier();
  }
}
