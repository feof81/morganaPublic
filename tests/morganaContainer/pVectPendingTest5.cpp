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

#include "pVect.hpp"
#include "pVectManip.hpp"
#include "pMapItem.h"
#include "pVectComm.hpp"

using namespace Teuchos;


//! Run with two processors
int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  RCP<communicator> world(new communicator);
  
  assert(world->size() == 2);
  
  typedef pVect<Real,pMapItem>     PVECT;
  typedef pVectComm<Real,pMapItem> PVCOMM;
  typedef PVCOMM::PVPS             PVPS;
  
  PVECT    sVector, rVector;    
  pMapItem item;
  
  pVectComm<Real,pMapItem> vectorComm(world);
  pMap<pMapItemSendRecv> sendingMap;
  pMap<pMapItemSendRecv> receivingMap;
  
  if(world->rank() == 0)
  {
    item.setPid(0); item.setLid(1); item.setGid(5); sVector.push_back(item,5.0,false);
    item.setPid(0); item.setLid(2); item.setGid(6); sVector.push_back(item,6.0,false);
    
    sVector.updateFinder();
    
    sendingMap.push_back(pMapItemSendRecv(1,5,0,1));
    receivingMap.push_back(pMapItemSendRecv(2,2,1,0));
    receivingMap.push_back(pMapItemSendRecv(3,3,1,0));
  }
  
  if(world->rank() == 1)
  {
    item.setPid(1); item.setLid(1); item.setGid(1); sVector.push_back(item,1.0,false);
    item.setPid(1); item.setLid(2); item.setGid(2); sVector.push_back(item,2.0,false);
    item.setPid(1); item.setLid(3); item.setGid(3); sVector.push_back(item,3.0,false);
    item.setPid(1); item.setLid(4); item.setGid(4); sVector.push_back(item,4.0,false);
    
    sVector.updateFinder();
    
    sendingMap.push_back(pMapItemSendRecv(2,2,1,0));
    sendingMap.push_back(pMapItemSendRecv(3,3,1,0));
    receivingMap.push_back(pMapItemSendRecv(1,5,0,1));
  }
  
  //Send and recv
  sVect<PVECT> commSegments;
  sVect<PVPS>  sendPvps;
  sVect<PVPS>  recvPvps;
  
  vectorComm.sendRecvI(sendingMap,
                       receivingMap,
                       sVector,
                       rVector,
                       commSegments,
                       sendPvps,
                       recvPvps);
  
  vectorComm.sendRecvO(sendingMap,
                       receivingMap,
                       sVector,
                       rVector,
                       commSegments,
                       sendPvps,
                       recvPvps);
  
  world->barrier();
  
  if(world->rank() == 0)
  { cout << rVector << endl; }
  
  
  world->barrier();
  
  if(world->rank() == 1)
  { cout << rVector << endl; }
}


/*
1
 map:  pid: 1 lid: 1 gid: 2
 data: 2
2
 map:  pid: 1 lid: 2 gid: 3
 data: 3
*/


/*
1
 map:  pid: 0 lid: 1 gid: 5
 data: 5
*/
