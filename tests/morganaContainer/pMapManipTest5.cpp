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
#include "pMapManip.hpp"

using namespace Teuchos;


/*! Run with one processor */
int main(int argc, char *argv[])
{
  typedef pMap<pMapItem> MAP;
  
  //The map
  pMapItem  item;
  RCP<MAP> map(new MAP);
    
  item.setLid(0); item.setGid(2);  map->push_back(item); //Build the map with 6 items
  item.setLid(0); item.setGid(4);  map->push_back(item);
  item.setLid(0); item.setGid(6);  map->push_back(item);
  item.setLid(0); item.setGid(8);  map->push_back(item);
  item.setLid(0); item.setGid(10); map->push_back(item);
  item.setLid(0); item.setGid(15); map->push_back(item);
  
  
  //Segmenter
  pMapManip<pMapItem> doctor;
  doctor.setNormalIndexing(map); //Set the lids
  
  cout << *map << endl;
}
