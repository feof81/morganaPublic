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
#include "staticVector.hpp"

using namespace Teuchos;


//! Run with two processors
int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  RCP<communicator> world(new communicator);
  
  typedef staticVector<10> DOFTYPE;
  
  time_t start, end;
  UInt numItems = 1e6;
  UInt numSend  = 100;
  sVect<DOFTYPE> dummy(numItems);
  
  assert(world->size() == 2);
  
  if(world->rank() == 0)
  {
    //Comm Data
    UInt sid = 0, rid = 1;
    
    dummy(1)[0] = 1;
    dummy(2)[0] = 2;
    
    //Skeleton
    boost::mpi::skeleton_proxy<sVect<DOFTYPE> > skeleton = boost::mpi::skeleton(dummy);
    
    sVect<request> reqs(1);
    reqs(1) = world->isend(rid,1,skeleton);
    wait_all(reqs.begin(),reqs.end());
    
    boost::mpi::content conSVect = boost::mpi::get_content(dummy);
    
    //Contenent
    cout << "Send time "; time(&start); cout << endl;

    for(UInt k=1; k <= numSend; ++k)
    {
      sVect<request> reqs(1);
      reqs(1) = world->isend(rid,1,conSVect);
      wait_all(reqs.begin(),reqs.end());
    }
    
    time(&end); cout << "Send Time (" << difftime(end, start) << " s)" << endl << endl;
  }
  
  if(world->rank() == 1)
  {
    //Comm Data    
    UInt sid = 0, rid = 1;
    
    //Skeleton
    boost::mpi::skeleton_proxy<sVect<DOFTYPE> > skeleton = boost::mpi::skeleton(dummy);
    
    sVect<request> reqs(1);
    reqs(1) = world->irecv(sid,1,skeleton);
    wait_all(reqs.begin(),reqs.end());
    
    boost::mpi::content conRVect = boost::mpi::get_content(dummy);
    
    //Contenent
    cout << "Recv time "; time(&start); cout << endl;

    for(UInt k=1; k <= numSend; ++k)
    {
      sVect<request> reqs(1);
      reqs(1) = world->irecv(sid,1,conRVect);
      wait_all(reqs.begin(),reqs.end());
    }
    
    time(&end); cout << "Recv Time (" << difftime(end, start) << " s)" << endl << endl;
    
    //cout << dummy << endl;
  }
}



