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

#include "pMapItem.h"
#include "pMapItemShare.h"

#include "pGraph.hpp"

#include "geoShapes.h"
#include "mesh2d.hpp"
#include "mesh2dGlobalManip.hpp"
#include "meshDoctor2d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;

int main(int argc, char *argv[])
{
  typedef linearQuad              GEOSHAPE2D;
  typedef geoElement<GEOSHAPE2D>  GEOELEMENT2D;
  typedef pMapItemShare ELMAP;
  typedef pMapItem      NODEMAP;
  
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  mesh2d<GEOSHAPE2D,ELMAP,NODEMAP> grid2d;
  
  
  if(world.rank() == 0)
  {
    //Maps
    NODEMAP nodeMapItem;
    ELMAP   elMapItem;
    
    nodeMapItem.setPid(0);
    elMapItem.setPid(0);
    
    //Nodes2d
    pVect<point3d,NODEMAP> nodes2d;
    
    nodes2d.reserve(9);
    nodeMapItem.setLid(1); nodeMapItem.setGid(2);
    nodes2d.push_back(point3d(0.0, 0.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(2); nodeMapItem.setGid(3);
    nodes2d.push_back(point3d(1.0, 0.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(3); nodeMapItem.setGid(5);
    nodes2d.push_back(point3d(0.0, 1.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(4); nodeMapItem.setGid(6);
    nodes2d.push_back(point3d(1.0, 1.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(5); nodeMapItem.setGid(8);
    nodes2d.push_back(point3d(0.0, 2.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(6); nodeMapItem.setGid(9);
    nodes2d.push_back(point3d(1.0, 2.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(7); nodeMapItem.setGid(7);
    nodes2d.push_back(point3d(2.0, 1.0, 0.0),nodeMapItem);

    nodeMapItem.setLid(8); nodeMapItem.setGid(10);
    nodes2d.push_back(point3d(2.0, 2.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(9); nodeMapItem.setGid(1);
    nodes2d.push_back(point3d(0.1, 0.1, 0.1),nodeMapItem);
    
    
    //Elements2d
    GEOELEMENT2D quad(true);
    pGraph<GEOELEMENT2D,ELMAP,NODEMAP> elList2d;
    
    elList2d.reserve(3);
    quad.setGeoId(1);
    
    elMapItem.setLid(1); elMapItem.setGid(1); elMapItem.setShared(false); elMapItem.setOwned(true);
    quad(1) = 1; quad(2) = 2; quad(3) = 4; quad(4) = 3;
    elList2d.push_back(quad,elMapItem);
    
    elMapItem.setLid(2); elMapItem.setGid(2); elMapItem.setShared(false); elMapItem.setOwned(true);
    quad(1) = 3; quad(2) = 4; quad(3) = 6; quad(4) = 5;
    elList2d.push_back(quad,elMapItem);
    
    elMapItem.setLid(3); elMapItem.setGid(4); elMapItem.setShared(true); elMapItem.setOwned(false);
    quad(1) = 4; quad(2) = 7; quad(3) = 8; quad(4) = 6;
    elList2d.push_back(quad,elMapItem);
    
    elList2d.updateRowFinder();
    elList2d.updateColFinder();
    
    
    //The grid2d
    grid2d.setNodes(nodes2d);
    grid2d.setElements(elList2d);
  }
  else
  {
    //Maps
    NODEMAP nodeMapItem;
    ELMAP   elMapItem;
    
    nodeMapItem.setPid(1);
    elMapItem.setPid(1);
    
    //Nodes2d
    pVect<point3d,NODEMAP> nodes2d;
    
    nodes2d.reserve(6);
    nodeMapItem.setLid(1); nodeMapItem.setGid(3);
    nodes2d.push_back(point3d(1.0, 0.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(2); nodeMapItem.setGid(4);
    nodes2d.push_back(point3d(2.0, 0.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(3); nodeMapItem.setGid(6);
    nodes2d.push_back(point3d(1.0, 1.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(4); nodeMapItem.setGid(7);
    nodes2d.push_back(point3d(2.0, 1.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(5); nodeMapItem.setGid(9);
    nodes2d.push_back(point3d(1.0, 2.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(6); nodeMapItem.setGid(10);
    nodes2d.push_back(point3d(2.0, 2.0, 0.0),nodeMapItem);
    
    
    //Elements2d
    GEOELEMENT2D quad(true);
    pGraph<GEOELEMENT2D,ELMAP,NODEMAP> elList2d;
    
    elList2d.reserve(2);
    quad.setGeoId(1);
    
    elMapItem.setLid(1); elMapItem.setGid(3); elMapItem.setShared(false); elMapItem.setOwned(true);
    quad(1) = 1; quad(2) = 2; quad(3) = 4; quad(4) = 3;
    elList2d.push_back(quad,elMapItem);
    
    elMapItem.setLid(2); elMapItem.setGid(4); elMapItem.setShared(true); elMapItem.setOwned(true);
    quad(1) = 3; quad(2) = 4; quad(3) = 6; quad(4) = 5;
    elList2d.push_back(quad,elMapItem);
    
    elList2d.updateRowFinder();
    elList2d.updateColFinder();
    
    
    //The grid2d
    grid2d.setNodes(nodes2d);
    grid2d.setElements(elList2d);
  }
  
  meshDoctor2d<GEOSHAPE2D,ELMAP,NODEMAP> doctor2d(world);
  bool checkIds = doctor2d.checkGeoIds(grid2d);
  UInt numIds   = doctor2d.countGeoIds(grid2d);
  
  assert(checkIds);
  doctor2d.removeUnusedPoints(grid2d);
  
  
  world.barrier();
  if(world.rank() == 0)
  {
    cout << "Num ids: " << numIds << endl;
    cout << "Nodes" << endl;
    cout << grid2d.getNodes() << endl << endl;
    
    cout << "Elements" << endl;
    cout << grid2d.getElements() << endl;
  }
  sleep(1);
  
  world.barrier();
  if(world.rank() == 1)
  {
    cout << "Num ids: " << numIds << endl;
    cout << "Nodes" << endl;
    cout << grid2d.getNodes() << endl << endl;
    
    cout << "Elements" << endl;
    cout << grid2d.getElements() << endl;
  }
  sleep(1);
}


/* PID 0 
GeoIds checking: 1
GeoIds checking: 1
Num GeoIds: 2
1
 map:  pid: 0 lid: 1 gid: 1
 data: 0 0 0

2
 map:  pid: 0 lid: 2 gid: 2
 data: 1 0 0

3
 map:  pid: 0 lid: 3 gid: 3
 data: 1 1 0


ROW
1
Num GeoIds: 2
 map:  pid: 0 lid: 1 gid: 1 shared: 0 owned 1
 data: 
GeoId          : 1
Num Connected  : 2
Connected Id's : 1 2 


2
 map:  pid: 0 lid: 2 gid: 2 shared: 0 owned 1
 data: 
GeoId          : 1
Num Connected  : 2
Connected Id's : 2 3 


COL
pid: 0 lid: 1 gid: 1
pid: 0 lid: 2 gid: 2
pid: 0 lid: 3 gid: 3 
 */


/* PID 1
GeoIds checking: 1
GeoIds checking: 1
Num GeoIds: 2
1
Num GeoIds: 2
 map:  pid: 1 lid: 1 gid: 1
 data: 0 0 0

2
 map:  pid: 1 lid: 2 gid: 3
 data: 1 1 0

3
 map:  pid: 1 lid: 3 gid: 4
 data: 0 1 0


ROW
1
 map:  pid: 1 lid: 1 gid: 3 shared: 0 owned 1
 data: 
GeoId          : 2
Num Connected  : 2
Connected Id's : 2 3 


2
 map:  pid: 1 lid: 2 gid: 4 shared: 0 owned 1
 data: 
GeoId          : 2
Num Connected  : 2
Connected Id's : 3 1 


COL
pid: 1 lid: 1 gid: 1
pid: 1 lid: 2 gid: 3
pid: 1 lid: 3 gid: 4 
 */
