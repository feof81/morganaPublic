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

#include "printMesh.hpp"
#include "dual2d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef linearTriangle GEOSHAPE2D;
  typedef linearLine     GEOSHAPE1D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>  MESH2D;
  typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>  MESH1D;
  
  
  string meshFile = "./tests/morganaMeshes/mignon2dC.msh";
  //string meshFile = "./geometries/cubes3d/testCubeA.msh";
  
  meshInit2d<GEOSHAPE2D,ELMAP,NODEMAP> init(world);
  init.gmMesh_to_stdA(meshFile);
  
  //Dual mesh
  dual2d<GEOSHAPE2D,ELMAP,NODEMAP> dualGrid(world);
  dualGrid.setMesh2d(init.getGrid2d());
  dualGrid.setConnect2d(init.getConnectGrid2d());
  
  dualGrid.buildDualMesh();
  dualGrid.buildFiniteVolumeData();
}
