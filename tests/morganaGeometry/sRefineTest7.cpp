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
#include "printMesh.hpp"
#include "printMesh2dHDF5.hpp"
#include "printMesh3dHDF5.hpp"

#include "sRefine3d.hpp"

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
  string meshFile = "./tests/morganaMeshes/mignon3dA.msh";
  
  meshInit3d<GEOSHAPE3D,ELMAP,NODEMAP> init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  
  //! Refinement-----------------------------------------------------
  sVect<UInt> elList(2);
  elList(1) = 1;
  elList(2) = 2;
  
  sRefine3d<GEOSHAPE3D,ELMAP,NODEMAP> refiner;
  refiner.upload(grid3d, grid2d);
  refiner.setRefinementParams(0.2, 2, elList);
  
  refiner.refineLeb(50);
  
  refiner.download(grid3d, grid2d);
  
  //! Printing-------------------------------------------------------
  grid3d->transferMap();
  
  cout << grid3d->getNodes() << endl;
  cout << grid3d->getElements() << endl;
  
  //Grid2d
  grid2d->transferMap();
  
  printMesh3dHDF5<GEOSHAPE3D,ELMAP> printer3d(world);
  printer3d.printGeoIds("griglia3d",0,grid3d);
  
  printMesh2dHDF5<GEOSHAPE2D,ELMAP> printer2d(world);
  printer2d.printGeoIds("griglia2d",0,grid2d);
}
