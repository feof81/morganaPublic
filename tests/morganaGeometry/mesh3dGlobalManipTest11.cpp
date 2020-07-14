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
#include "mesh3d.hpp"

#include "mesh3dGlobalManip.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;

int main(int argc, char *argv[])
{
  typedef pMapItem             ELMAP;
  typedef pMapItemShare        NODEMAP;
  typedef linearHexa           GEOSHAPE;
  typedef geoElement<GEOSHAPE> GEOELEMENT;
  
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  
  GEOELEMENT                                 tet(GEOELEMENT::numPoints);
  pVect<point3d,NODEMAP>                     nodes;
  pGraph<GEOELEMENT,ELMAP,NODEMAP>           elList;
  mesh3d<GEOSHAPE,ELMAP,NODEMAP>             grid3d;
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP>  gridManipulator(world);
  
  
  if(world.rank() == 0)
  {
    //Nodes
    nodes.reserve(16);
    nodes.push_back(point3d(0.0, 0.0, 0.0),pMapItemShare(1,1,0,false,true));
    nodes.push_back(point3d(1.0, 0.0, 0.0),pMapItemShare(2,2,0,true,true));
    nodes.push_back(point3d(2.0, 0.0, 0.0),pMapItemShare(3,3,0,true,true));
    nodes.push_back(point3d(3.0, 0.0, 0.0),pMapItemShare(4,4,0,true,false));
    
    nodes.push_back(point3d(0.0, 1.0, 0.0),pMapItemShare(5,6,0,false,true));
    nodes.push_back(point3d(1.0, 1.0, 0.0),pMapItemShare(6,7,0,true,true));
    nodes.push_back(point3d(2.0, 1.0, 0.0),pMapItemShare(7,8,0,true,true));
    nodes.push_back(point3d(3.0, 1.0, 0.0),pMapItemShare(8,9,0,true,false));
    
    nodes.push_back(point3d(0.0, 0.0, 1.0),pMapItemShare(9, 11,0,false,true));
    nodes.push_back(point3d(1.0, 0.0, 1.0),pMapItemShare(10,12,0,true,true));
    nodes.push_back(point3d(2.0, 0.0, 1.0),pMapItemShare(11,13,0,true,true));
    nodes.push_back(point3d(3.0, 0.0, 1.0),pMapItemShare(12,14,0,true,false));
    
    nodes.push_back(point3d(0.0, 1.0, 1.0),pMapItemShare(13,16,0,false,true));
    nodes.push_back(point3d(1.0, 1.0, 1.0),pMapItemShare(14,17,0,true,true));
    nodes.push_back(point3d(2.0, 1.0, 1.0),pMapItemShare(15,18,0,true,true));
    nodes.push_back(point3d(3.0, 1.0, 1.0),pMapItemShare(16,19,0,true,false));
    nodes.updateFinder();
    
    //Elements  
    elList.reserve(3);
    
    tet.setGeoId(1);
    tet(1) = 1; tet(2) = 2; tet(3) = 6; tet(4) = 5; tet(5) = 9; tet(6) = 10; tet(7) = 14; tet(8) = 13;
    elList.push_back(tet, pMapItem(1,1,0));
    
    tet.setGeoId(1);
    tet(1) = 2; tet(2) = 3; tet(3) = 7; tet(4) = 6; tet(5) = 10; tet(6) = 11; tet(7) = 15; tet(8) = 14;
    elList.push_back(tet, pMapItem(2,3,0));
    
    tet.setGeoId(2);
    tet(1) = 3; tet(2) = 4; tet(3) = 8; tet(4) = 7; tet(5) = 11; tet(6) = 12; tet(7) = 16; tet(8) = 15;
    elList.push_back(tet, pMapItem(3,2,0));
    
    //The grid3d
    grid3d.setElements(elList);
    grid3d.setNodes(nodes);
  }
  else
  {
    //Nodes
    nodes.reserve(16);
    nodes.push_back(point3d(1.0, 0.0, 0.0),pMapItemShare(1,2,1,true,false));
    nodes.push_back(point3d(2.0, 0.0, 0.0),pMapItemShare(2,3,1,true,false));
    nodes.push_back(point3d(3.0, 0.0, 0.0),pMapItemShare(3,4,1,true,true));
    nodes.push_back(point3d(4.0, 0.0, 0.0),pMapItemShare(4,5,1,false,true));
    
    nodes.push_back(point3d(1.0, 1.0, 0.0),pMapItemShare(5,7,1,true,false));
    nodes.push_back(point3d(2.0, 1.0, 0.0),pMapItemShare(6,8,1,true,false));
    nodes.push_back(point3d(3.0, 1.0, 0.0),pMapItemShare(7,9,1,true,true));
    nodes.push_back(point3d(4.0, 1.0, 0.0),pMapItemShare(8,10,1,false,true));
    
    nodes.push_back(point3d(1.0, 0.0, 1.0),pMapItemShare(9, 12,1,true,false));
    nodes.push_back(point3d(2.0, 0.0, 1.0),pMapItemShare(10,13,1,true,false));
    nodes.push_back(point3d(3.0, 0.0, 1.0),pMapItemShare(11,14,1,true,true));
    nodes.push_back(point3d(4.0, 0.0, 1.0),pMapItemShare(12,15,1,false,true));
    
    nodes.push_back(point3d(1.0, 1.0, 1.0),pMapItemShare(13,17,1,true,false));
    nodes.push_back(point3d(2.0, 1.0, 1.0),pMapItemShare(14,18,1,true,false));
    nodes.push_back(point3d(3.0, 1.0, 1.0),pMapItemShare(15,19,1,true,true));
    nodes.push_back(point3d(4.0, 1.0, 1.0),pMapItemShare(16,20,1,false,true));
    nodes.updateFinder();
    
    //Elements  
    elList.reserve(3);
    
    tet.setGeoId(1);
    tet(1) = 1; tet(2) = 2; tet(3) = 6; tet(4) = 5; tet(5) = 9; tet(6) = 10; tet(7) = 14; tet(8) = 13;
    elList.push_back(tet, pMapItem(1,3,1));
    
    tet.setGeoId(1);
    tet(1) = 2; tet(2) = 3; tet(3) = 7; tet(4) = 6; tet(5) = 10; tet(6) = 11; tet(7) = 15; tet(8) = 14;
    elList.push_back(tet, pMapItem(2,2,1));
    
    tet.setGeoId(2);
    tet(1) = 3; tet(2) = 4; tet(3) = 8; tet(4) = 7; tet(5) = 11; tet(6) = 12; tet(7) = 16; tet(8) = 15;
    elList.push_back(tet, pMapItem(3,4,1));
    
    //The grid3d
    grid3d.setElements(elList);
    grid3d.setNodes(nodes);
  }
  
  grid3d.setMeshStandard(STDA);

  
  //Mesh overlapping
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> balancer(world);  
  balancer.destroyOverlap(grid3d);
  
  
  world.barrier();
  if(world.rank() == 0)
  {
    cout << grid3d.getElements() << endl;
    cout << grid3d.getNodes() << endl;
  }
  sleep(1);
  
  world.barrier();
  if(world.rank() == 1)
  {
    cout << grid3d.getElements() << endl;
    cout << grid3d.getNodes() << endl;
  }
  sleep(1);
}


/*PID 0
ROW
1
 map:  pid: 0 lid: 1 gid: 1
 data: 
GeoId          : 1
Num Connected  : 8
Connected Id's : 1 2 5 4 7 8 11 10 


2
 map:  pid: 0 lid: 2 gid: 3
 data: 
GeoId          : 1
Num Connected  : 8
Connected Id's : 2 3 6 5 8 9 12 11 


COL
pid: 0 lid: 1 gid: 1 shared: 0 owned 1
pid: 0 lid: 2 gid: 2 shared: 0 owned 1
pid: 0 lid: 3 gid: 3 shared: 1 owned 1
pid: 0 lid: 4 gid: 6 shared: 0 owned 1
pid: 0 lid: 5 gid: 7 shared: 0 owned 1
pid: 0 lid: 6 gid: 8 shared: 1 owned 1
pid: 0 lid: 7 gid: 11 shared: 0 owned 1
pid: 0 lid: 8 gid: 12 shared: 0 owned 1
pid: 0 lid: 9 gid: 13 shared: 1 owned 1
pid: 0 lid: 10 gid: 16 shared: 0 owned 1
pid: 0 lid: 11 gid: 17 shared: 0 owned 1
pid: 0 lid: 12 gid: 18 shared: 1 owned 1


1
 map:  pid: 0 lid: 1 gid: 1 shared: 0 owned 1
 data: 0 0 0

2
 map:  pid: 0 lid: 2 gid: 2 shared: 0 owned 1
 data: 1 0 0

3
 map:  pid: 0 lid: 3 gid: 3 shared: 1 owned 1
 data: 2 0 0

4
 map:  pid: 0 lid: 4 gid: 6 shared: 0 owned 1
 data: 0 1 0

5
 map:  pid: 0 lid: 5 gid: 7 shared: 0 owned 1
 data: 1 1 0

6
 map:  pid: 0 lid: 6 gid: 8 shared: 1 owned 1
 data: 2 1 0

7
 map:  pid: 0 lid: 7 gid: 11 shared: 0 owned 1
 data: 0 0 1

8
 map:  pid: 0 lid: 8 gid: 12 shared: 0 owned 1
 data: 1 0 1

9
 map:  pid: 0 lid: 9 gid: 13 shared: 1 owned 1
 data: 2 0 1

10
 map:  pid: 0 lid: 10 gid: 16 shared: 0 owned 1
 data: 0 1 1

11
 map:  pid: 0 lid: 11 gid: 17 shared: 0 owned 1
 data: 1 1 1

12
 map:  pid: 0 lid: 12 gid: 18 shared: 1 owned 1
 data: 2 1 1 
 */


/*PID 1
ROW
1
 map:  pid: 1 lid: 1 gid: 2
 data: 
GeoId          : 1
Num Connected  : 8
Connected Id's : 1 2 5 4 7 8 11 10 


2
 map:  pid: 1 lid: 2 gid: 4
 data: 
GeoId          : 2
Num Connected  : 8
Connected Id's : 2 3 6 5 8 9 12 11 


COL
pid: 1 lid: 1 gid: 3 shared: 1 owned 0
pid: 1 lid: 2 gid: 4 shared: 0 owned 1
pid: 1 lid: 3 gid: 5 shared: 0 owned 1
pid: 1 lid: 4 gid: 8 shared: 1 owned 0
pid: 1 lid: 5 gid: 9 shared: 0 owned 1
pid: 1 lid: 6 gid: 10 shared: 0 owned 1
pid: 1 lid: 7 gid: 13 shared: 1 owned 0
pid: 1 lid: 8 gid: 14 shared: 0 owned 1
pid: 1 lid: 9 gid: 15 shared: 0 owned 1
pid: 1 lid: 10 gid: 18 shared: 1 owned 0
pid: 1 lid: 11 gid: 19 shared: 0 owned 1
pid: 1 lid: 12 gid: 20 shared: 0 owned 1


1
 map:  pid: 1 lid: 1 gid: 3 shared: 1 owned 0
 data: 2 0 0

2
 map:  pid: 1 lid: 2 gid: 4 shared: 0 owned 1
 data: 3 0 0

3
 map:  pid: 1 lid: 3 gid: 5 shared: 0 owned 1
 data: 4 0 0

4
 map:  pid: 1 lid: 4 gid: 8 shared: 1 owned 0
 data: 2 1 0

5
 map:  pid: 1 lid: 5 gid: 9 shared: 0 owned 1
 data: 3 1 0

6
 map:  pid: 1 lid: 6 gid: 10 shared: 0 owned 1
 data: 4 1 0

7
 map:  pid: 1 lid: 7 gid: 13 shared: 1 owned 0
 data: 2 0 1

8
 map:  pid: 1 lid: 8 gid: 14 shared: 0 owned 1
 data: 3 0 1

9
 map:  pid: 1 lid: 9 gid: 15 shared: 0 owned 1
 data: 4 0 1

10
 map:  pid: 1 lid: 10 gid: 18 shared: 1 owned 0
 data: 2 1 1

11
 map:  pid: 1 lid: 11 gid: 19 shared: 0 owned 1
 data: 3 1 1

12
 map:  pid: 1 lid: 12 gid: 20 shared: 0 owned 1
 data: 4 1 1 
 */
