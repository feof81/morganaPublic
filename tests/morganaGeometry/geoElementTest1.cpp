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

#include "geoShapes.h"
#include "geoElement.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;

int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  cout << "I am: " << world.rank() << endl;
  
  if(world.rank() == 0)
  {
    geoElement<linearTetra> tet(true), recTet(true);
    
    tet.getGeoId() = 1;
    tet(1) = 1;
    tet(2) = 3;
    tet(3) = 5;
    
    cout << "Rank 0 sending" << endl;
    world.send(1,0,tet);
    cout << "Rank 0 reciving" << endl;
    world.recv(1,1,recTet);
    
    //cout << "Process: " << world.rank() << " Size: " << world.size() << endl;
    //cout << recTet << endl;
  }
  else
  {
    geoElement<linearTetra> tet(true), recTet(true);
    
    tet.getGeoId() = 2;
    tet(1) = 2;
    tet(2) = 4;
    tet(3) = 5;
    
    cout << "Rank 1 reciving" << endl;
    world.recv(0,0,recTet);
    cout << "Rank 1 sending" << endl;
    world.send(0,1,tet);
    
    cout << "Process: " << world.rank() << " Size: " << world.size() << endl;
    cout << recTet << endl;
  }
}


/* PID 0
I am: 0
Rank 0 sending
I am: 1
Rank 1 reciving
Rank 0 reciving
Rank 1 sending
Process: 0 Size: 2
GeoId          : 2
Num Connected  : 4
Connected Id's : 2 4 5 0
*/


/* PID 1
I am: 1
Rank 1 reciving
I am: 0
Rank 0 sending
Rank 1 sending
Rank 0 reciving
Process: 1 Size: 2
GeoId          : 1
Num Connected  : 4
Connected Id's : 1 3 5 0
*/

