/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"

#include <cmath>
#include <iostream>

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "pMap.hpp"
#include "pMapItem.h"
#include "pMpiOptimization.hpp"
#include "pMapComm.hpp"
#include "pMapGlobalManip.h"

using namespace Teuchos;



//! Run with two processors
int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  
  //Alloc
  typedef pMapItem      MAPITEM;
  typedef pMap<MAPITEM> PMAP;
  pMapItem item;
  PMAP mappa;
  
  //The Morgana map
  if(world.rank() == 0)
  {
    item.setLid(1); item.setGid(2);  mappa.push_back(item);
    item.setLid(2); item.setGid(4);  mappa.push_back(item);
  }
  else
  {
    item.setLid(1); item.setGid(1);  mappa.push_back(item);
    item.setLid(2); item.setGid(3);  mappa.push_back(item);
  }
  

  //The Tpetra map
  typedef int                                                     ORDINALTYPE;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType            PLATFORM;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  NODE;
  typedef const Tpetra::Map<ORDINALTYPE,ORDINALTYPE,NODE>         TPETRA_MAP;
  
  PLATFORM & tpetraPlatform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<TPETRA_MAP> tpetraMap;
  
  pMapGlobalManip<MAPITEM> manipulator(world);
  manipulator.exportTpetraMap(mappa, tpetraMap);
  
  cout << tpetraMap->getNodeElementList() << endl;
}
