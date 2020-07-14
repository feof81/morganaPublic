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
#include "meshDoctor3d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;

int main(int argc, char *argv[])
{
  typedef linearTetra            GEOSHAPE3D;
  typedef geoElement<GEOSHAPE3D> GEOELEMENT3D;
  typedef pMapItemShare          ELMAP;
  typedef pMapItem               NODEMAP;
  
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  mesh3d<GEOSHAPE3D,ELMAP,NODEMAP> grid3d;
  
  
  if(world.rank() == 0)
  {
    //Nodes 3d
    pVect<point3d,NODEMAP> nodes3d;
    
    nodes3d.reserve(5);
    nodes3d.push_back(point3d(0.3, 0.3, 0.3),pMapItem(1,1)); nodes3d.getMapL(1).setPid(0);
    nodes3d.push_back(point3d(0.0, 0.0, 0.0),pMapItem(2,2)); nodes3d.getMapL(2).setPid(0);
    nodes3d.push_back(point3d(1.0, 0.0, 0.0),pMapItem(3,3)); nodes3d.getMapL(3).setPid(0);
    nodes3d.push_back(point3d(1.0, 1.0, 0.0),pMapItem(4,4)); nodes3d.getMapL(4).setPid(0);
    nodes3d.push_back(point3d(0.0, 1.0, 0.0),pMapItem(5,5)); nodes3d.getMapL(5).setPid(0);
    nodes3d.push_back(point3d(0.5, 0.5, 1.0),pMapItem(6,6)); nodes3d.getMapL(6).setPid(0);
    nodes3d.updateFinder();
    
    //Elements 3d
    GEOELEMENT3D tet(true);
    pGraph<GEOELEMENT3D,ELMAP,NODEMAP> elList3d;
    
    elList3d.reserve(1);
    tet.setGeoId(1); tet(1) = 2; tet(2) = 3; tet(3) = 4; tet(4) = 6;
    elList3d.push_back(tet, pMapItemShare(1,1,false,true)); elList3d.getRowMapL(1).setPid(0);
    
    //The grid3d
    grid3d.setNodes(nodes3d);
    grid3d.setElements(elList3d);
  }
  else
  {
    //Nodes 3d
    pVect<point3d,NODEMAP> nodes3d;
    
    nodes3d.reserve(5);
    nodes3d.push_back(point3d(0.4, 0.4, 0.4),pMapItem(1,1)); nodes3d.getMapL(1).setPid(1);
    nodes3d.push_back(point3d(0.0, 0.0, 0.0),pMapItem(2,2)); nodes3d.getMapL(2).setPid(1);
    nodes3d.push_back(point3d(1.0, 0.0, 0.0),pMapItem(3,3)); nodes3d.getMapL(3).setPid(1);
    nodes3d.push_back(point3d(1.0, 1.0, 0.0),pMapItem(4,4)); nodes3d.getMapL(4).setPid(1);
    nodes3d.push_back(point3d(0.0, 1.0, 0.0),pMapItem(5,5)); nodes3d.getMapL(5).setPid(1);
    nodes3d.push_back(point3d(0.5, 0.5, 1.0),pMapItem(6,6)); nodes3d.getMapL(6).setPid(1);
    nodes3d.updateFinder();
    
    //Elements 3d
    GEOELEMENT3D tet(true);
    pGraph<GEOELEMENT3D,ELMAP,NODEMAP> elList3d;
    
    elList3d.reserve(1);
    tet.setGeoId(2); tet(1) = 2; tet(2) = 4; tet(3) = 5; tet(4) = 6;
    elList3d.push_back(tet, pMapItemShare(1,2,false,true)); elList3d.getRowMapL(1).setPid(1);
    
    //The grid3d
    grid3d.setNodes(nodes3d);
    grid3d.setElements(elList3d);
  }
  
  grid3d.transferMap();
  
  meshDoctor3d<linearTetra,ELMAP,NODEMAP> doctor(world);
  doctor.removeUnusedPoints(grid3d);
  
  
  world.barrier();
  if(world.rank() == 0)
  {
    cout << "Nodes" << endl;
    cout << grid3d.getNodes() << endl << endl;
    
    cout << "Elements" << endl;
    cout << grid3d.getElements() << endl;
  }
  sleep(1);
  
  world.barrier();
  if(world.rank() == 1)
  {
    cout << "Nodes" << endl;
    cout << grid3d.getNodes() << endl << endl;
    
    cout << "Elements" << endl;
    cout << grid3d.getElements() << endl;
  }
  sleep(1);
}


/* PID 0 
Nodes
1
 map:  pid: 0 lid: 1 gid: 1
 data: 0 0 0

2
 map:  pid: 0 lid: 2 gid: 2
 data: 1 0 0

3
 map:  pid: 0 lid: 3 gid: 3
 data: 1 1 0

4
 map:  pid: 0 lid: 4 gid: 5
 data: 0.5 0.5 1



Elements
ROW
1
 map:  pid: 0 lid: 1 gid: 1 shared: 0 owned 1
 data: 
GeoId          : 1
Num Connected  : 4
Connected Id's : 1 2 3 4 


COL
pid: 0 lid: 1 gid: 1
pid: 0 lid: 2 gid: 2
pid: 0 lid: 3 gid: 3
pid: 0 lid: 4 gid: 5 
 */


/* PID 1
Nodes
1
 map:  pid: 1 lid: 1 gid: 1
 data: 0 0 0

2
 map:  pid: 1 lid: 2 gid: 3
 data: 1 1 0

3
 map:  pid: 1 lid: 3 gid: 4
 data: 0 1 0

4
 map:  pid: 1 lid: 4 gid: 5
 data: 0.5 0.5 1



Elements
ROW
1
 map:  pid: 1 lid: 1 gid: 2 shared: 0 owned 1
 data: 
GeoId          : 2
Num Connected  : 4
Connected Id's : 1 2 3 4 


COL
pid: 1 lid: 1 gid: 1
pid: 1 lid: 2 gid: 3
pid: 1 lid: 3 gid: 4
pid: 1 lid: 4 gid: 5 
 */

