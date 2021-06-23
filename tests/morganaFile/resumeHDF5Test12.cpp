/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "resumeHDF5.h"


using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);  
  
  sVect<UInt> indices;
  sVect<Real> serialized;
  sVect<sVect<Real> > inArray(3), outArray;
  
  if(world.rank() == 0)
  {
    inArray(1).resize(1);
    inArray(2).resize(2);
    inArray(3).resize(3);
  
    inArray(1)(1) = 1.0;
    inArray(2)(1) = 2.0; inArray(2)(2) = 3.0;
    inArray(3)(1) = 4.0; inArray(3)(2) = 5.0; inArray(3)(3) = 6.0;
  }
  else
  {
    inArray(3).resize(1);
    inArray(2).resize(2);
    inArray(1).resize(3);
    
    inArray(3)(1) = 1.0;
    inArray(2)(1) = 2.0; inArray(2)(2) = 3.0;
    inArray(1)(1) = 4.0; inArray(1)(2) = 5.0; inArray(1)(3) = 6.0;
  }
  
  resumeHDF5 resumer(world);
  resumer.printToFile("array",inArray);
  resumer.loadFromFile("array",outArray);
  
  cout << outArray << endl;
}
