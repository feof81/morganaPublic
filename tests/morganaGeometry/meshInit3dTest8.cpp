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
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef linearHexa     GEOSHAPE;
  typedef linearQuad     GEOSHAPE2D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  typedef mesh3d<GEOSHAPE,ELMAP,NODEMAP>    MESH3D;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>  MESH2D;
  
  
  string meshFile  = "../geometries/hexa3d/hexa3dA.neu";
  string colorFile = "../geometries/hexa3d/hexa3dA_color.neu";
  
  meshInit3d<GEOSHAPE,ELMAP,NODEMAP> init(world);
  init.neutral_to_stdB(meshFile, colorFile);
  
  
  //Local printing 3d
  RCP<MESH3D> grid3d = init.getGrid3d(); 
  string local3d     = "init3dTest8_localMesh3d_" + num2str(world.rank()) + "pid.inp";
  UInt pid           = world.rank();
  
  printMesh<GEOSHAPE,ELMAP,NODEMAP> localPrinter3d;
  localPrinter3d.paraviewLocal(local3d, pid, *grid3d);
  
  
  //Local printing 2d
  RCP<MESH2D> grid2d = init.getGrid2d();
  string local2d = "init3dTest8_localMesh2d_" + num2str(world.rank()) + "pid.inp";
  
  printMesh<GEOSHAPE2D,ELMAP,NODEMAP> localPrinter2d;
  localPrinter2d.paraviewLocal(local2d, pid, *grid2d);
  
  
  //Global printing 3d
  UInt printPid = 0;
  
  mesh3d<GEOSHAPE,ELMAP,NODEMAP>            globalGrid3d;
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gathering3d(world);
  gathering3d.gather(printPid,globalGrid3d,*grid3d);
  
  printMesh<GEOSHAPE,ELMAP,NODEMAP> meshPrinter3d;
  meshPrinter3d.paraviewSerial("init3dTest8_globalMesh3d.inp",globalGrid3d);
  
  
  //Global printing 2d
  mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>            globalGrid2d;
  mesh2dGlobalManip<GEOSHAPE2D,ELMAP,NODEMAP> gathering2d(world);
  gathering2d.gather(printPid,globalGrid2d,*grid2d);
  
  printMesh<GEOSHAPE2D,ELMAP,NODEMAP> meshPrinter2d;
  meshPrinter2d.paraviewSerial("init3dTest8_globalMesh2d.inp",globalGrid2d);
}
