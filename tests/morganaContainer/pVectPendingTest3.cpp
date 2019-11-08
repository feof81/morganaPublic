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
  //Mpi--------------------------------------------------------------
  environment  env(argc,argv);
  RCP<communicator> world(new communicator);
  
  //Alloc------------------------------------------------------------
  typedef pVect<point3d,pMapItem>     PVECT;
  typedef pVectComm<point3d,pMapItem> PVCOMM;
  typedef PVCOMM::PVPS                PVPS;
  
  PVECT    sVector, rVector;    
  pMapItem item;
  time_t start, end;
  
  PVCOMM vectorComm(world);
  PVPS   commData;
  
  assert(world->size() == 2);
  
  //Init-------------------------------------------------------------
  if(world->rank() == 0)
  {
    for(UInt k=1; k <= 20000000; ++k)
    {
      item.setPid(0);
      item.setLid(k);
      item.setGid(4+k);
      sVector.push_back(item,point3d(4+k,0.0,0.0),false);
    }
  }
  
  if(world->rank() == 1)
  {
    item.setPid(1); item.setLid(1); item.setGid(1); sVector.push_back(item,point3d(1.0,0.0,0.0),false);
    item.setPid(1); item.setLid(2); item.setGid(2); sVector.push_back(item,point3d(2.0,0.0,0.0),false);
    item.setPid(1); item.setLid(3); item.setGid(3); sVector.push_back(item,point3d(3.0,0.0,0.0),false);
    item.setPid(1); item.setLid(4); item.setGid(4); sVector.push_back(item,point3d(4.0,0.0,0.0),false);
  }
  
  //Comm-------------------------------------------------------------
  if(world->rank() == 0) {cout << "Comm "; time(&start); cout << endl;}
  
  UInt sid = 0, rid = 1;
  
  vectorComm.sendRecvI(sid,
                       rid,
                       sVector,
                       rVector,
                       commData);
  
  if(world->rank() == 0) {time(&end); cout << "Comm done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Computing--------------------------------------------------------
  if(world->rank() == 0) {cout << "Comput "; time(&start); cout << endl;}
  
  sVect<UInt> numeriPrimi;
  bool isPrime;
  
  for(UInt k=1; k <= 100000; k++)
  {
    isPrime = true;
    
    for(UInt j=2; j < k; ++j)
    {
      if((k % j) == 0)
      {
        isPrime = false;
        break;
      }
    }
    
    if(isPrime) { numeriPrimi.push_back(k); }
  }
  
  if(world->rank() == 0) {time(&end); cout << "Comput done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Comm-------------------------------------------------------------
  if(world->rank() == 0) {cout << "Comm "; time(&start); cout << endl;}
  
  vectorComm.sendRecvO(sid,
                       rid,
                       sVector,
                       rVector,
                       commData);
  
  if(world->rank() == 0) {time(&end); cout << "Comm done (" << difftime(end, start) << " s)" << endl << endl;}
}

