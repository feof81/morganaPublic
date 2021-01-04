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

using namespace std;
using namespace boost::mpi;


/*! Run with two processors */
int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  
  if(world.rank() == 0)
  {
    pMapItem Item; //Create a map to be set
    pMap<pMapItem> sMap;
    
    Item.setLid(1); Item.setGid(1); //Insert the first map item
    sMap.push_back(Item);
    
    Item.setLid(2); Item.setGid(2); //Insert the second map item
    sMap.push_back(Item);
    
    world.send(1,0,sMap); //Send the map
  }
  
  if(world.rank() == 1)
  {
    pMap<pMapItem> rMap; //Create a map to be received
    
    cout << "Rank 1 receiving" << endl; //Receive the map
    world.recv(0,0,rMap);
    
    cout << "Process: " << world.rank() << " Size: " << world.size() << endl;
    cout << rMap << endl;

  }
}

/*
Rank 1 reciving
Rank 0 sending
Process: 1 Size: 2
pid: 0 lid: 1 gid: 1
pid: 0 lid: 2 gid: 2
*/
