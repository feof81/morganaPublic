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
  
  typedef pVect<Real,pMapItem> PVECT;
  PVECT    sVector, rVector;    
  pMapItem item;
  
  pVectComm<Real,pMapItem> vectorComm(world);
  
  assert(world->size() == 2);
  
  if(world->rank() == 0)
  {
    item.setPid(0); item.setLid(1); item.setGid(5); sVector.push_back(item,5.0,false);
    item.setPid(0); item.setLid(2); item.setGid(6); sVector.push_back(item,6.0,false);
    
    UInt sid = 0, rid = 1;
    vectorComm.sendRecv(sid,rid,sVector,rVector);
  }
  
  if(world->rank() == 1)
  {
    item.setPid(1); item.setLid(1); item.setGid(1); sVector.push_back(item,1.0,false);
    item.setPid(1); item.setLid(2); item.setGid(2); sVector.push_back(item,2.0,false);
    item.setPid(1); item.setLid(3); item.setGid(3); sVector.push_back(item,3.0,false);
    item.setPid(1); item.setLid(4); item.setGid(4); sVector.push_back(item,4.0,false);
    
    UInt sid = 0, rid = 1;
    vectorComm.sendRecv(sid,rid,sVector,rVector);
    cout << rVector << endl;
  }
}

/* Pid 1
1
 map:  pid: 0 lid: 1 gid: 5
 data: 5
2
 map:  pid: 0 lid: 2 gid: 6
 data: 6
*/

