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

#include "meshDoctor3d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;

int main(int argc, char *argv[])
{
  typedef linearTetra          GEOSHAPE;
  typedef geoElement<GEOSHAPE> GEOELEMENT;
  typedef pMapItemShare        ELMAP;
  typedef pMapItem             NODEMAP;
  
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  mesh3d<GEOSHAPE,ELMAP,NODEMAP> grid3d;
  
  if(world.rank() == 0)
  {
    //Nodes
    pVect<point3d,NODEMAP> nodes;
    
    nodes.reserve(5);
    nodes.push_back(point3d(0.0, 0.0, 0.0),pMapItem(1,1));
    nodes.push_back(point3d(1.0, 0.0, 0.0),pMapItem(2,2));
    nodes.push_back(point3d(1.0, 1.0, 0.0),pMapItem(3,3));
    nodes.push_back(point3d(0.0, 1.0, 0.0),pMapItem(4,4));
    nodes.push_back(point3d(0.5, 0.5, 1.0),pMapItem(5,5));
    nodes.updateFinder();
    
    //Elements
    GEOELEMENT tet(true);
    pGraph<GEOELEMENT,ELMAP,NODEMAP> elList;
    
    elList.reserve(2);
    
    tet.setGeoId(1); tet(1) = 1; tet(2) = 2; tet(3) = 3; tet(4) = 5;
    elList.push_back(tet, pMapItemShare(1,1,true,true));
    tet.setGeoId(2); tet(1) = 1; tet(2) = 3; tet(3) = 4; tet(4) = 5;
    elList.push_back(tet, pMapItemShare(2,2,true,false));
    
    //The grid3d
    grid3d.setNodes(nodes);
    grid3d.setElements(elList);
  }
  else
  {
    //Nodes
    pVect<point3d,NODEMAP> nodes;
    
    nodes.reserve(5);
    nodes.push_back(point3d(0.0, 0.0, 0.0),pMapItem(1,1));
    nodes.push_back(point3d(1.0, 0.0, 0.0),pMapItem(2,2));
    nodes.push_back(point3d(1.0, 1.0, 0.0),pMapItem(3,3));
    nodes.push_back(point3d(0.0, 1.0, 0.0),pMapItem(4,4));
    nodes.push_back(point3d(0.5, 0.5, 1.0),pMapItem(5,5));
    nodes.updateFinder();
    
    //Elements
    GEOELEMENT tet(true);
    pGraph<GEOELEMENT,ELMAP,NODEMAP> elList;
    
    elList.reserve(2);
    
    tet.setGeoId(2); tet(1) = 1; tet(2) = 2; tet(3) = 3; tet(4) = 5;
    elList.push_back(tet, pMapItemShare(1,1,true,false));
    tet.setGeoId(3); tet(1) = 1; tet(2) = 3; tet(3) = 4; tet(4) = 5;
    elList.push_back(tet, pMapItemShare(2,2,true,true));
    
    //The grid3d
    grid3d.setNodes(nodes);
    grid3d.setElements(elList);
  }
  
  
  //Mesh doctor
  meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP> doctor(world);
  bool flag = doctor.checkGeoIds(grid3d);
  
  cout << " flag: " << flag << endl;
}
