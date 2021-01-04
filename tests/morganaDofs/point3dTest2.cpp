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

#include "point3d.h"
#include "typesInterface.hpp"
#include "traitsMpiOptimization.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  point3d recP;
  
  //Pass
  if(world.rank() == 0)
  {
    point3d P(1.0,0.0,0.0);
    
    world.send(1,0,P);
    world.recv(1,1,recP);
  }
  else
  {
    point3d P(0.0,0.0,1.0);
    
    world.recv(0,0,recP);
    world.send(0,1,P);
  }
  
  //Printout
  if(world.rank() == 0)
  {
    cout << "Process: " << world.rank() << " Size: " << world.size() << endl;
    cout << recP << endl;
  }
  
  world.barrier();
  
  if(world.rank() == 1)
  {
    cout << "Process: " << world.rank() << " Size: " << world.size() << endl;
    cout << recP << endl;
  }
}
