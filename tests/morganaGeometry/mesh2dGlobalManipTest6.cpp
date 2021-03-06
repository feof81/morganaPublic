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

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;

int main(int argc, char *argv[])
{
  typedef linearTriangle       GEOSHAPE;
  typedef geoElement<GEOSHAPE> GEOELEMENT;
  typedef pMapItemShare        ELMAP;
  typedef pMapItemShare        NODEMAP;
  
  environment  env(argc,argv);
  communicator world;
  
  mesh2d<GEOSHAPE,ELMAP,NODEMAP> grid2d;
  assert(world.size() == 2);
  
  sVect<UInt> elWeights(3);
  
  if(world.rank() == 0)
  {
    //Nodes
    pVect<point3d,NODEMAP> nodes;
    
    nodes.reserve(4);
    nodes.push_back(point3d(0.0, 0.0, 0.0),pMapItemShare(1,1,false,true));
    nodes.push_back(point3d(1.0, 0.0, 0.0),pMapItemShare(2,2,false,true));
    nodes.push_back(point3d(2.0, 0.0, 0.0),pMapItemShare(3,3,false,true));
    nodes.push_back(point3d(3.0, 0.0, 0.0),pMapItemShare(4,4,false,true));
    nodes.updateFinder();
    
    //Elements
    GEOELEMENT tri(true);
    pGraph<GEOELEMENT,ELMAP,NODEMAP> elList;
    
    elList.reserve(3);
    
    tri.setGeoId(1); tri(1) = 1; tri(2) = 2; tri(3) = 6;
    elList.push_back(tri, pMapItemShare(1,1,false,true));
    
    tri.setGeoId(1); tri(1) = 1; tri(2) = 6; tri(3) = 5;
    elList.push_back(tri, pMapItemShare(2,2,false,true));
    
    tri.setGeoId(1); tri(1) = 2; tri(2) = 7; tri(3) = 6;
    elList.push_back(tri, pMapItemShare(3,5,false,true));
    
    elList.colIsLocal() = false;
    
    //Weights 
    elWeights(1) = 1.0;
    elWeights(2) = 1.0;
    elWeights(3) = 10.0;
    
    //The grid2d
    grid2d.setNodes(nodes);
    grid2d.setElements(elList);
    grid2d.setMeshStandard(STDL);
  }
  else
  {
    //Nodes
    pVect<point3d,NODEMAP> nodes;
    
    nodes.reserve(4);
    nodes.push_back(point3d(0.0, 1.0, 0.0),pMapItemShare(1,5,false,true));
    nodes.push_back(point3d(1.0, 1.0, 0.0),pMapItemShare(2,6,false,true));
    nodes.push_back(point3d(2.0, 1.0, 0.0),pMapItemShare(3,7,false,true));
    nodes.push_back(point3d(3.0, 1.0, 0.0),pMapItemShare(4,8,false,true));
    nodes.updateFinder();
    
    //Elements
    GEOELEMENT tri(true);
    pGraph<GEOELEMENT,ELMAP,NODEMAP> elList;
    
    elList.reserve(3);
    
    tri.setGeoId(1); tri(1) = 2; tri(2) = 3; tri(3) = 7;
    elList.push_back(tri, pMapItemShare(1,6,false,true));
    
    tri.setGeoId(1); tri(1) = 3; tri(2) = 8; tri(3) = 7;
    elList.push_back(tri, pMapItemShare(2,3,false,true));
    
    tri.setGeoId(1); tri(1) = 3; tri(2) = 4; tri(3) = 8;
    elList.push_back(tri, pMapItemShare(3,4,false,true));
    
    elList.colIsLocal() = false;
    
    //Weights 
    elWeights(1) = 10.0;
    elWeights(2) = 10.0;
    elWeights(3) = 10.0;
    
    //The grid2d
    grid2d.setNodes(nodes);
    grid2d.setElements(elList);
    grid2d.setMeshStandard(STDL);
  }
  
  grid2d.transferMap();
  
  mesh2dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> manipulator(world);
  manipulator.meshPartition(grid2d);
  
  bool flag = manipulator.check(grid2d);
  assert(flag);
  
  manipulator.meshBalancing(grid2d,elWeights);
  
  
  
  world.barrier();
  if(world.rank() == 0)
  {
    cout << grid2d.getElements() << endl;
    cout << grid2d.getNodes() << endl;
  }
  sleep(1);
  
  world.barrier();
  if(world.rank() == 1)
  {
    cout << grid2d.getElements() << endl;
    cout << grid2d.getNodes() << endl;
  }
  sleep(1);
}
