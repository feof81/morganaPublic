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

#include "traitsMpiOptimization.hpp"
#include "point3d.h"
#include "simpleFormats.hpp"
#include "pVect.hpp"
#include "pMapItem.h"


using namespace std;
using namespace boost::mpi;


//! To be run with two processors: mpirun -np 2 ./bin/morgana
int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  pVect<point3d,pMapItem> pVectSend, pVectRecv;
  
  if(world.rank() == 0)
  {
    point3d P1(2.0,0.0,0.0), P2(0.0,2.0,0.0);
    
    pVectSend.push_back(P1,pMapItem(1,3));
    pVectSend.push_back(P2,pMapItem(2,4));
    pVectSend.updateFinder();
    
    
    world.send(1,0,pVectSend);
    world.recv(1,1,pVectRecv);
  }
  else
  {
    point3d P1(1.0,0.0,0.0), P2(0.0,1.0,0.0);
    
    pVectSend.push_back(P1,pMapItem(1,1));
    pVectSend.push_back(P2,pMapItem(2,2));
    pVectSend.updateFinder();
    
    world.recv(0,0,pVectRecv);
    world.send(0,1,pVectSend);
  }
  
  
  world.barrier();
  if(world.rank() == 0)
  {
    cout << "Process: " << world.rank() << " Size: " << world.size() << endl;
    cout << pVectRecv << endl;
  }
  
  world.barrier();
  if(world.rank() == 1)
  {
    cout << "Process: " << world.rank() << " Size: " << world.size() << endl;
    cout << pVectRecv << endl;
  }
}


/* Pid 0
Process: 0 Size: 2
1
 map:  pid: 0 lid: 1 gid: 1
 data: 1 0 0

2
 map:  pid: 0 lid: 2 gid: 2
 data: 0 1 0
*/


/* Pid 1
Process: 1 Size: 2
1
 map:  pid: 0 lid: 1 gid: 3
 data: 2 0 0

2
 map:  pid: 0 lid: 2 gid: 4
 data: 0 2 0
*/
