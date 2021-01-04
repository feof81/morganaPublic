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

#include "komplex.h"
#include "typesInterface.hpp"
#include "traitsMpiOptimization.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;

int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  komplex recC;
  
  //Compute
  if(world.rank() == 0)
  {
    request reqs[2];
    komplex C(1.0,0.0);
    
    reqs[0] = world.isend(1,0,C);
    reqs[1] = world.irecv(1,1,recC);
    
    wait_all(reqs,reqs+2);
  }
  else
  {
    request reqs[2];
    komplex C(0.0,1.0);
    
    reqs[0] = world.isend(0,1,C);
    reqs[1] = world.irecv(0,0,recC);
    
    wait_all(reqs,reqs+2);
  }
  
  //Printout
  if(world.rank() == 0)
  {
    cout << "Process: " << world.rank() << " Size: " << world.size() << endl;
    cout << recC << endl;
  }
  
  world.barrier();
  
  if(world.rank() == 1)
  {
    cout << "Process: " << world.rank() << " Size: " << world.size() << endl;
    cout << recC << endl;
  }
}
