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
  
  pMap<pMapItem> map, rMap;
  resumeHDF5 resumer(world);
   
  if(world.rank() == 0)
  {
    map.push_back(pMapItem(1,2,0));
    map.push_back(pMapItem(2,4,0));
    map.push_back(pMapItem(3,5,0));
    map.push_back(pMapItem(4,3,0));
  }
  
  if(world.rank() == 1)
  {
    map.push_back(pMapItem(1,1,1));
  }
  
  map.bufferLids();
  
  resumer.printToFile("map",map);
  resumer.loadFromFile("map",rMap);
  
  if(world.rank() == 0)
  { 
    cout <<  map << endl;
    cout << rMap << endl;
  }
}
