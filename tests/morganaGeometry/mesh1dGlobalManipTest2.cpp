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
#include "mesh1d.hpp"
#include "mesh1dGlobalManip.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;

int main(int argc, char *argv[])
{
  typedef geoElement<linearLine> GEOELEMENT;
  typedef pMapItem   ELMAP;
  typedef pMapItem   NODEMAP;
  
  ELMAP   elMap;
  NODEMAP nodeMap;
  
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  mesh1d<linearLine,ELMAP,NODEMAP> grid1d;
  
  if(world.rank() == 0)
  {
    //Nodes
    pVect<point3d,NODEMAP> nodes;
    
    nodes.reserve(2);
    
    nodeMap.setLid(1); nodeMap.setGid(1); nodeMap.setPid(0);
    nodes.push_back(point3d(0.0, 0.0, 0.0),nodeMap);
    nodeMap.setLid(2); nodeMap.setGid(2); nodeMap.setPid(0);
    nodes.push_back(point3d(1.0, 0.0, 0.0),nodeMap);
    
    //Elements
    GEOELEMENT li(true);
    pGraph<GEOELEMENT,ELMAP,NODEMAP> elList;
    
    elList.reserve(2);
    nodeMap.setLid(1); nodeMap.setGid(1); nodeMap.setPid(0);
    li.setGeoId(1); li(1) = 1; li(2) = 2;
    elList.push_back(li, nodeMap);
    nodeMap.setLid(2); nodeMap.setGid(2); nodeMap.setPid(0);
    li.setGeoId(1); li(1) = 2; li(2) = 3;
    elList.push_back(li, nodeMap);
    
    elList.colIsLocal() = false;
    
    //The grid2d
    grid1d.setNodes(nodes);
    grid1d.setElements(elList);
    grid1d.setMeshStandard(STDL);
  }
  else
  {
    //Nodes    
    pVect<point3d,NODEMAP> nodes;
    
    nodes.reserve(2);
    
    elMap.setLid(1); elMap.setGid(3); elMap.setPid(1);
    nodes.push_back(point3d(1.0, 1.0, 0.0),elMap);
    elMap.setLid(2); elMap.setGid(4); elMap.setPid(1);
    nodes.push_back(point3d(0.0, 1.0, 0.0),elMap);
    
    //Elements
    GEOELEMENT li(true);
    pGraph<GEOELEMENT,ELMAP,NODEMAP> elList;
    
    elList.reserve(2);
    
    nodeMap.setLid(1); nodeMap.setGid(3); nodeMap.setPid(1);
    li.setGeoId(2); li(1) = 3; li(2) = 4;
    elList.push_back(li, nodeMap);   
    nodeMap.setLid(2); nodeMap.setGid(4); nodeMap.setPid(1);
    li.setGeoId(2); li(1) = 4; li(2) = 1;
    elList.push_back(li, nodeMap);
    
    elList.colIsLocal() = false;
    
    //The grid2d
    grid1d.setNodes(nodes);
    grid1d.setElements(elList);
    grid1d.setMeshStandard(STDL);
  }
  
  grid1d.transferMap();
  
  //Mehs manipulator
  mesh1dGlobalManip<linearLine,ELMAP,NODEMAP> lineManipulator(world);
  
  lineManipulator.meshPartition(grid1d);
  bool flag = lineManipulator.check(grid1d);
  
  
  world.barrier();
  if(world.rank() == 0)
  {
    cout << "Check: " << flag << endl;
    cout << grid1d.getNodes() << endl;
    cout << grid1d.getElements() << endl;
  }
  sleep(1);
  
  world.barrier();
  if(world.rank() == 1)
  {
    cout << "Check: " << flag << endl;
    cout << grid1d.getNodes() << endl;
    cout << grid1d.getElements() << endl;
  }
  sleep(1);
}


/* PID 0
 Check: 1
1
 map:  pid: 0 lid: 1 gid: 1
Check: 1
 data: 0 0 0

2
 map:  pid: 0 lid: 2 gid: 2
 data: 1 0 0

3
 map:  pid: 0 lid: 3 gid: 3
 data: 1 1 0


ROW
1
 map:  pid: 0 lid: 1 gid: 1
 data: 
GeoId          : 1
Num Connected  : 2
Connected Id's : 1 2 


2
 map:  pid: 0 lid: 2 gid: 2
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
 Check: 1
Check: 1
1
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
 map:  pid: 1 lid: 1 gid: 3
 data: 
GeoId          : 2
Num Connected  : 2
Connected Id's : 2 3 


2
 map:  pid: 1 lid: 2 gid: 4
 data: 
GeoId          : 2
Num Connected  : 2
Connected Id's : 3 1 


COL
pid: 1 lid: 1 gid: 1
pid: 1 lid: 2 gid: 3
pid: 1 lid: 3 gid: 4
*/
