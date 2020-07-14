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

#include "pMapItem.h"
#include "pMapItemShare.h"

#include "geoShapes.h"
#include "meshInit3d.hpp"

#include "connect2d.hpp"
#include "connect3d.hpp"
#include "geoOctTree.hpp"


// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 1);
  
  typedef linearTetra    GEOSHAPE3D;
  typedef linearTriangle GEOSHAPE2D;
  typedef linearLine     GEOSHAPE1D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  typedef mesh3d<GEOSHAPE3D,ELMAP,NODEMAP>         MESH3D;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>         MESH2D;
  typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>         MESH1D;
  typedef connect3d<GEOSHAPE3D,ELMAP,NODEMAP>      CONNECT3D;
  typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>      CONNECT2D;
  
  
  //! Mesh Init------------------------------------------------------
  string meshFile = "./geometries/cubes3d/testCubeE.msh";
  
  meshInit3d<GEOSHAPE3D,ELMAP,NODEMAP> init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  
  //! Search 3d------------------------------------------------------
  typedef geoOctTree<MESH2D,CONNECT2D> OCTTREE;
  
  UInt size = 0;
  bool test = true;
  sVect<point3d> elemNodes;
  sVect<UInt>    elements;
  std::set<UInt> elemList;
  
  OCTTREE octTree(grid2d,connectGrid2d);
  octTree.localInit(1.5);
  
  for(UInt i=1; i <= grid2d->getNumElements(); ++i)
  {
    elemNodes = grid2d->getElementNodesL(i);
    elements  = octTree.getMatchingElements( (elemNodes(1) + elemNodes(2) + elemNodes(3)) / 3.0);
    
    size = std::max(size, UInt(elements.size()));
    elemList.clear();
    
    for(UInt k=1; k <= elements.size(); ++k)
    { elemList.insert(elements(k)); }
    
    test = test && (elemList.count(i) == 1);
  }
  
  cout << "Test-1: " << test << endl;
  test = true;
  
  for(UInt i=1; i <= grid2d->getNumElements(); ++i)
  {
    elemNodes = grid2d->getElementNodesL(i);
    elements  = octTree.getMatchingElements(elemNodes(1));
    
    elemList.clear();
    
    for(UInt k=1; k <= elements.size(); ++k)
    { elemList.insert(elements(k)); }
    
    test = test && (elemList.count(i) == 1);
  }
  
  cout << "Test-2: " << test << endl;
  
  
  cout << "Outpoint: " << octTree.getMatchingElements(point3d(-1.0,-1.0,-1.0)).size() << endl;
  cout << "Size: " << size << endl;
}
