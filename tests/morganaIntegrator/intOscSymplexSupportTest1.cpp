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

#include "intOscSymplexSupport.h"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;


int main(int argc, char *argv[])
{
  point3d   K( 3.0,  3.0, 0.0);  
  point3d  P1( 0.0,  0.0, 0.0);
  point3d  P2( 1.0,  0.0, 0.0);
  point3d  P3( 0.0,  1.0, 0.0);
  point3d  P4( 0.0,  0.0, 1.0);
  
  cout << "S0 1d: " << intOscSymplexSupport::intS0_1d(K,P1,P2);
  cout << "S0 2d: " << intOscSymplexSupport::intS0_2d(K,P1,P2,P3);
  cout << "S0 3d: " << intOscSymplexSupport::intS0_3d(K,P1,P2,P3,P4) << endl;
}

