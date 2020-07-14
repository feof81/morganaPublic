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
#include "printMesh2dHDF5.hpp"
#include "printMesh3dHDF5.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  //Typedefs---------------------------------------------------------
  typedef linearTetra     GEOSHAPE;
  typedef linearTriangle  GEOSHAPE2D;
  typedef pMapItemShare   PMAPTYPE;
  
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> INIT;  
  typedef INIT::MESH3D     MESH3D;
  typedef INIT::MESH2D     MESH2D;
  typedef INIT::CONNECT3D  CONNECT3D;
  typedef INIT::CONNECT2D  CONNECT2D;
  
  
  //Mesh loading-----------------------------------------------------
  string meshFile = "./geometries/biCube3d/biCubeC.msh";
  //string meshFile = "./geometries/cubes3d/testCubeB.msh";
  
  INIT init(world);
  init.gmMesh_to_stdA(meshFile);
  
  //Meshes links-----------------------------------------------------
  Teuchos::RCP<MESH3D>           grid3d = init.getGrid3d();
  Teuchos::RCP<MESH2D>           grid2d = init.getGrid2d();
  Teuchos::RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  Teuchos::RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  //Mesh rebalance---------------------------------------------------
  set<UInt> activeGeoIds;
  activeGeoIds.insert(1);
  //activeGeoIds.insert(2);
  
  init.reinit_stdA(grid2d,
                   grid3d,
                   activeGeoIds);
  
  //Meshes links-----------------------------------------------------
  grid3d        = init.getGrid3d();
  grid2d        = init.getGrid2d();
  connectGrid3d = init.getConnectGrid3d();
  connectGrid2d = init.getConnectGrid2d();
  
  //Printer----------------------------------------------------------
  string local2d = "localMesh2d_" + num2str(world.rank()) + "pid.inp";
  string local3d = "localMesh3d_" + num2str(world.rank()) + "pid.inp";
  UInt       pid = world.rank();
  
  printMesh<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> localPrinter2d;
  localPrinter2d.paraviewLocal(local2d, pid, *grid2d);
  
  printMesh<GEOSHAPE,PMAPTYPE,PMAPTYPE> localPrinter3d;
  localPrinter3d.paraviewLocal(local3d, pid, *grid3d);
  
  printMesh3dHDF5<GEOSHAPE,PMAPTYPE> printer3d(world);
  printer3d.printGeoIds("grid3d",0,grid3d);
}
