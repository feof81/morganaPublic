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
  typedef linearTriangle   GEOSHAPE;
  typedef geoElement<GEOSHAPE> ITEM;
  typedef pMapItem           ROWMAP;
  typedef pMapItem           COLMAP;
  
  environment  env(argc,argv);
  RCP<communicator> world(new communicator);
  
  assert(world->size() == 2);
  
  ITEM grphItem(true);
  ROWMAP pItem;
  
  pMap<COLMAP> colMap;
  pGraph<ITEM,ROWMAP,COLMAP> grafo, rGraph;
  resumeHDF5 resumer(world);
  
  if(world->rank() == 0)
  {
    pItem.setPid(0);
    
    pItem.setLid(1); pItem.setGid(1); colMap.push_back(pItem);
    pItem.setLid(2); pItem.setGid(2); colMap.push_back(pItem);
    pItem.setLid(3); pItem.setGid(3); colMap.push_back(pItem);
    pItem.setLid(4); pItem.setGid(4); colMap.push_back(pItem);
    pItem.setLid(5); pItem.setGid(5); colMap.push_back(pItem);
    pItem.setLid(6); pItem.setGid(6); colMap.push_back(pItem);
  
    pItem.setLid(1); pItem.setGid(1); grphItem(1) = 1; grphItem(2) = 5; grphItem(3) = 6; grphItem.setGeoId(1); grafo.push_back(pItem,grphItem);
    pItem.setLid(2); pItem.setGid(2); grphItem(1) = 1; grphItem(2) = 2; grphItem(3) = 5; grphItem.setGeoId(1); grafo.push_back(pItem,grphItem);
    pItem.setLid(3); pItem.setGid(3); grphItem(1) = 2; grphItem(2) = 4; grphItem(3) = 5; grphItem.setGeoId(1); grafo.push_back(pItem,grphItem);
    pItem.setLid(4); pItem.setGid(4); grphItem(1) = 2; grphItem(2) = 3; grphItem(3) = 4; grphItem.setGeoId(1); grafo.push_back(pItem,grphItem);
  
    grafo.setColMap(colMap);
    grafo.updateRowFinder();
    grafo.updateColFinder();
  }
  else
  {
    pItem.setPid(1);
    
    pItem.setLid(1); pItem.setGid(2); colMap.push_back(pItem);
    pItem.setLid(2); pItem.setGid(3); colMap.push_back(pItem);
    pItem.setLid(3); pItem.setGid(7); colMap.push_back(pItem);
    pItem.setLid(4); pItem.setGid(8); colMap.push_back(pItem);
    pItem.setLid(5); pItem.setGid(4); colMap.push_back(pItem);
    pItem.setLid(6); pItem.setGid(5); colMap.push_back(pItem);
    
    pItem.setLid(1); pItem.setGid(5); grphItem(1) = 1; grphItem(2) = 5; grphItem(3) = 6; grphItem.setGeoId(2); grafo.push_back(pItem,grphItem);
    pItem.setLid(2); pItem.setGid(6); grphItem(1) = 1; grphItem(2) = 2; grphItem(3) = 5; grphItem.setGeoId(2); grafo.push_back(pItem,grphItem);
    pItem.setLid(3); pItem.setGid(7); grphItem(1) = 2; grphItem(2) = 4; grphItem(3) = 5; grphItem.setGeoId(2); grafo.push_back(pItem,grphItem);
    pItem.setLid(4); pItem.setGid(8); grphItem(1) = 2; grphItem(2) = 3; grphItem(3) = 4; grphItem.setGeoId(2); grafo.push_back(pItem,grphItem);
    
    grafo.setColMap(colMap);
    grafo.updateRowFinder();
    grafo.updateColFinder();
  }
  
  grafo.pushToGlobal();
  
  resumer.printToFile("pGraph",grafo);
  resumer.loadFromFile("pGraph",rGraph);
  
  if(world->rank() == 1)
  {
    cout <<  grafo << endl;
    cout << rGraph << endl;
  }
}
