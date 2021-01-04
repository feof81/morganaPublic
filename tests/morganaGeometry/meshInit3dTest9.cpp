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

#include "printMesh.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using Teuchos::RCP;
using Teuchos::rcp;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef linearTetra    GEOSHAPE;
  typedef linearTriangle GEOSHAPE2D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  typedef mesh3d<GEOSHAPE,ELMAP,NODEMAP>    MESH3D;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>  MESH2D;
  
  
  //Mesh Init
  string meshFile = "../geometries/biCube3d/biCubeA.msh";
  
  meshInit3d<GEOSHAPE,ELMAP,NODEMAP> init(world);
  init.gmMesh_to_stdB(meshFile);
  
  
  //Local printing 3d
  RCP<MESH3D> grid3d = init.getGrid3d(); 
  string local3d     = "init3dTest9_localMesh3d_" + num2str(world.rank()) + "pid.inp";
  UInt pid           = world.rank();
  
  printMesh<GEOSHAPE,ELMAP,NODEMAP> localPrinter3d;
  localPrinter3d.paraviewLocal(local3d, pid, *grid3d);
  
  
  //Local printing 2d
  RCP<MESH2D> grid2d = init.getGrid2d();
  string local2d = "init3dTest9_localMesh2d_" + num2str(world.rank()) + "pid.inp";
  
  printMesh<GEOSHAPE2D,ELMAP,NODEMAP> localPrinter2d;
  localPrinter2d.paraviewLocal(local2d, pid, *grid2d);
  
  
  //Global printing 3d
  UInt printPid = 0;
  
  mesh3d<GEOSHAPE,ELMAP,NODEMAP>            globalGrid3d;
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gathering(world);
  gathering.gather(printPid,globalGrid3d,*grid3d);
  
  printMesh<GEOSHAPE,ELMAP,NODEMAP> meshPrinter;
  meshPrinter.paraviewSerial("init3dTest9_globalMesh.inp",globalGrid3d);
}
