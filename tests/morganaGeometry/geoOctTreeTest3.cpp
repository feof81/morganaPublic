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
#include "meshInit2d.hpp"

#include "connect1d.hpp"
#include "connect2d.hpp"
#include "geoOctTree.hpp"


using namespace std;
using namespace boost::mpi;
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef linearTriangle GEOSHAPE2D;
  typedef linearLine     GEOSHAPE1D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>     MESH2D;
  typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>     MESH1D;
  typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>  CONNECT2D;
  typedef connect1d<GEOSHAPE1D,ELMAP,NODEMAP>  CONNECT1D;
  typedef searchData<ELMAP>  SEARCHDATA;
  
  //Mesh Init
  UInt pid = world.rank();
  time_t start, end;
  string meshFile  = "./geometries/rectangles2d/rectangleD.msh";
  
  meshInit2d<GEOSHAPE2D,ELMAP,NODEMAP> init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH1D>           grid1d = init.getGrid1d();
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT1D> connectGrid1d = init.getConnectGrid1d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  
  //! Search 2d------------------------------------------------------
  typedef geoOctTree<MESH2D,CONNECT2D> OCTTREE;
  
  UInt size = 0;
  bool test = true;
  sVect<point3d> elemNodes;
  sVect<UInt>    elements;
  std::set<UInt> elemList;
  
  OCTTREE octTree(grid2d,connectGrid2d);
  octTree.localInit(3.0);
  
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
