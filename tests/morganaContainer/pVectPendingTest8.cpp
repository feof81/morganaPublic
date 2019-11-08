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
  
  typedef pVect<Real,pMapItem>     PVECT;
  typedef pVectComm<Real,pMapItem> PVCOMM;
  typedef PVCOMM::PVRS             PVRS;
  PVECT    sVector, rVector;    
  pMapItem item;
  
  PVCOMM vectorComm(world);
  PVRS   pvrs;
  
  assert(world->size() == 2);
  
  UInt sid = 0, rid = 1;
  UInt numLines = 10;
  
  if(world->rank() == 0)
  {
    item.setPid(0); 
    
    for(UInt i=1; i <= numLines; ++i)
    {
      item.setLid(i);
      item.setGid(i);
      sVector.push_back(item,Real(i),false);
    }
  }
 
  //Startup of the recursive comm channel
  sVect<request> reqsRR;
  
  vectorComm.sendRR(sid,
                    rid,
                    sVector,
                    pvrs,
                    reqsRR);
  
  vectorComm.recvRR(sid,
                    rid,
                    sVector,
                    pvrs);
  
  wait_all(reqsRR.begin(),reqsRR.end());
  
  //Begin of the comm
  vectorComm.sendRI(sid,
                    rid,
                    sVector,
                    pvrs);
  
  vectorComm.recvRI(sid,
                    rid,
                    sVector,
                    pvrs);
  
  //End of the comm
  vectorComm.sendRO(sid,
                    rid,
                    sVector,
                    pvrs);
  
  vectorComm.recvRO(sid,
                    rid,
                    sVector,
                    pvrs);

  //Printout
  if(world->rank() == 1)
  { cout << sVector << endl; }
}


