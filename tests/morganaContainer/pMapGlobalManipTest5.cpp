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
#include "pMapItemShare.h"
#include "pMpiOptimization.hpp"
#include "pMapComm.hpp"
#include "pMapGlobalManip.h"

using namespace Teuchos;



//! Run with three processors
int main(int argc, char *argv[])
{
  //Comunication stuff
  environment  env(argc,argv);
  RCP<communicator> world(new communicator);
  
  assert(world->size() == 3);
  
  //Common data
  typedef pMap<pMapItemShare>    OBJ;
  typedef pMap<pMapItemSendRecv> SENDRECV;
  
  pMapItemShare  item;
  RCP<OBJ>  sMap(new OBJ);
  RCP<OBJ>  rMap(new OBJ);
  
  if(world->rank() == 1)
  {  
    item.setLid(1); item.setGid(1);  item.setOwned(true); item.setShared(false);  sMap->push_back(item);
    item.setLid(2); item.setGid(2);  item.setOwned(true); item.setShared(false);  sMap->push_back(item);
    item.setLid(3); item.setGid(3);  item.setOwned(true); item.setShared(false);  sMap->push_back(item);
    item.setLid(4); item.setGid(4);  item.setOwned(true); item.setShared(true);   sMap->push_back(item);
    item.setLid(5); item.setGid(5);  item.setOwned(true); item.setShared(true);   sMap->push_back(item);
  }
  
  if(world->rank() == 2)
  {
    item.setLid(1); item.setGid(4);  item.setOwned(false); item.setShared(true);  sMap->push_back(item);
    item.setLid(2); item.setGid(5);  item.setOwned(false); item.setShared(true);  sMap->push_back(item);
    item.setLid(3); item.setGid(6);  item.setOwned(true);  item.setShared(false); sMap->push_back(item);
    item.setLid(4); item.setGid(7);  item.setOwned(true);  item.setShared(false); sMap->push_back(item);
    item.setLid(5); item.setGid(8);  item.setOwned(true);  item.setShared(false); sMap->push_back(item);
    item.setLid(6); item.setGid(9);  item.setOwned(true);  item.setShared(true);  sMap->push_back(item);
    item.setLid(7); item.setGid(10); item.setOwned(true);  item.setShared(true);  sMap->push_back(item);
  }
  
  if(world->rank() == 0)
  {
    item.setLid(1); item.setGid(9);  item.setOwned(false); item.setShared(true);  sMap->push_back(item);
    item.setLid(2); item.setGid(10); item.setOwned(false); item.setShared(true);  sMap->push_back(item);
    item.setLid(3); item.setGid(11); item.setOwned(true);  item.setShared(false); sMap->push_back(item);
    item.setLid(4); item.setGid(12); item.setOwned(true);  item.setShared(false); sMap->push_back(item);
    item.setLid(5); item.setGid(13); item.setOwned(true);  item.setShared(false); sMap->push_back(item);
    item.setLid(6); item.setGid(14); item.setOwned(true);  item.setShared(false); sMap->push_back(item);
    item.setLid(7); item.setGid(15); item.setOwned(true);  item.setShared(false); sMap->push_back(item);
  }
  
  
  pMapGlobalManip<pMapItemShare> manipulator(world);
  
  RCP<SENDRECV> mapSend(new SENDRECV);
  RCP<SENDRECV> mapRecv(new SENDRECV);
  
  manipulator.createSendRecvMap(sMap,mapSend,mapRecv);
  
  world->barrier();
  if(world->rank() == 0)
  {
    cout << "Sending Map" << endl;
    cout << *mapSend << endl;
    
    cout << "Reciving Map" << endl;
    cout << *mapRecv << endl;
  }
  sleep(1.0);
  
  world->barrier();
  if(world->rank() == 1)
  {
    cout << "Sending Map" << endl;
    cout << *mapSend << endl;
    
    cout << "Reciving Map" << endl;
    cout << *mapRecv << endl;
  }
  sleep(1.0);
  
  world->barrier();
  if(world->rank() == 2)
  {
    cout << "Sending Map" << endl;
    cout << *mapSend << endl;
    
    cout << "Reciving Map" << endl;
    cout << *mapRecv << endl;
  }
}

/* pid 0
Sending Map

Reciving Map
pid: 0 lid: 1 gid: 9 sid: 2 rid 0
pid: 0 lid: 2 gid: 10 sid: 2 rid 0
*/

/* pid 1
Sending Map
pid: 1 lid: 4 gid: 4 sid: 1 rid 2
pid: 1 lid: 5 gid: 5 sid: 1 rid 2

Reciving Map
*/

/* pid 2
Sending Map
pid: 2 lid: 6 gid: 9 sid: 2 rid 0
pid: 2 lid: 7 gid: 10 sid: 2 rid 0

Reciving Map
pid: 2 lid: 1 gid: 4 sid: 1 rid 2
pid: 2 lid: 2 gid: 5 sid: 1 rid 2
*/

