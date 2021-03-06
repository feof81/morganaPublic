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

#include "typesInterface.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;

int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  stateVector recV;
  
  //Compute
  if(world.rank() == 0)
  {
    request reqs[2];
    stateVector V(3);
    
    V(1) = 1.0;
    V(2) = 2.0;
    V(3) = 3.0;
    
    reqs[0] = world.isend(1,0,V);
    reqs[1] = world.irecv(1,1,recV);
    
    wait_all(reqs,reqs+2);
  }
  else
  {
    request reqs[2];
    stateVector V(3);
    
    V(1) = 3.0;
    V(2) = 2.0;
    V(3) = 1.0;
    
    reqs[0] = world.isend(0,1,V);
    reqs[1] = world.irecv(0,0,recV);
    
    wait_all(reqs,reqs+2);
  }
  
  
  //Print
  if(world.rank() == 0)
  {
    cout << "Process: " << world.rank() << " Size: " << world.size() << endl;
    cout << recV << endl;
  }
  
  world.barrier();
  
  
  if(world.rank() == 1)
  {
    cout << "Process: " << world.rank() << " Size: " << world.size() << endl;
    cout << recV << endl;
  }
}
