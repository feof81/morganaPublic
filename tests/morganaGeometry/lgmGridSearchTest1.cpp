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
#include "lgmGridSearch.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  //Communicator-----------------------------------------------------
  environment  env(argc,argv);
  communicator world;
  
  //Typedefs---------------------------------------------------------
  typedef linearTetra    GEOSHAPE3D;
  typedef linearTriangle GEOSHAPE2D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  typedef mesh3d<GEOSHAPE3D,ELMAP,NODEMAP>     MESH3D;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>     MESH2D;
  typedef connect3d<GEOSHAPE3D,ELMAP,NODEMAP>  CONNECT3D;
  typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>  CONNECT2D;
  typedef meshInit3d<GEOSHAPE3D,ELMAP,NODEMAP> MESHINIT;
  typedef lgmGridSearch<MESH2D,MESH3D>         GRIDSEARCH;
  typedef GRIDSEARCH::INVECT  INVECT;
  typedef GRIDSEARCH::OUTVECT OUTVECT;
  typedef GRIDSEARCH::INDATA  INDATA;
  
  //Grids------------------------------------------------------------
  string meshFile = "./tests/morganaMeshes/mignon3dA.msh";
  
  MESHINIT init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  //Gen input vector-------------------------------------------------
  INVECT inVect;
  INDATA inData;
  
  for(UInt el=1; el <= grid2d->getNumElements(); ++el)
  {
    inData.setPoints(grid2d->getElementNodesL(el));
    inVect.push_back(grid2d->getElements().getMapL(el),inData);
  }
  
  inVect.updateFinder();
  
  sVect<point3d> inLocCoords;
  inLocCoords.push_back(point3d(0.25,0.25,0.0));
  
  //Test-------------------------------------------------------------
  GRIDSEARCH gridMatcher(grid2d,
                         connectGrid2d,
                         grid3d,
                         connectGrid3d);
  
  OUTVECT outVect = gridMatcher.localSearch(inLocCoords,inVect);
  
  cout << outVect << endl;
}