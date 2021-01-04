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
#include "meshDoctor2d.hpp"
#include "meshDoctor3d.hpp"

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
  typedef linearTetra     GEOSHAPE3D;
  typedef linearTriangle  GEOSHAPE2D;
  typedef pMapItemShare   PMAPTYPE;
  
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> INIT;  
  typedef INIT::MESH3D     MESH3D;
  typedef INIT::MESH2D     MESH2D;
  typedef INIT::CONNECT3D  CONNECT3D;
  typedef INIT::CONNECT2D  CONNECT2D;
  
  typedef meshDoctor2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> DOCTOR2D;
  typedef meshDoctor3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> DOCTOR3D;
  
  
  //Mesh loading-----------------------------------------------------
  string meshFile = "./geometries/biCube3d/biCubeC.msh";
  //string meshFile = "./geometries/cubes3d/testCubeB.msh";
  
  INIT init(world);
  init.gmMesh_to_stdA(meshFile);
  
  //Meshes links-----------------------------------------------------
  Teuchos::RCP<MESH3D>           oldGrid3d = init.getGrid3d();
  Teuchos::RCP<MESH2D>           oldGrid2d = init.getGrid2d();
  Teuchos::RCP<CONNECT3D> oldConnectGrid3d = init.getConnectGrid3d();
  Teuchos::RCP<CONNECT2D> oldConnectGrid2d = init.getConnectGrid2d();
  
  //Mesh refinement--------------------------------------------------
  Teuchos::RCP<DOCTOR2D> doctor2d(new DOCTOR2D(world));
  Teuchos::RCP<DOCTOR3D> doctor3d(new DOCTOR3D(world));
  
  std::set<UInt> activeGeoIds = doctor3d->getGeoIds(oldGrid3d);
  
  Teuchos::RCP<MESH3D> newMesh3d(new MESH3D);
  doctor3d->refineUniform(newMesh3d,
                          oldGrid3d,  
                          oldConnectGrid3d);
  
  Teuchos::RCP<MESH2D> newMesh2d(new MESH2D);
  doctor2d->refineUniform(newMesh2d,
                          oldGrid2d,  
                          oldConnectGrid2d);
  
  init.reinit_stdA(newMesh2d,
                   newMesh3d,
                   activeGeoIds);
  
  //Printing---------------------------------------------------------
  typedef printMesh2dHDF5<GEOSHAPE2D,PMAPTYPE> PRINTER2D;
  typedef printMesh3dHDF5<GEOSHAPE3D,PMAPTYPE> PRINTER3D;
  
  PRINTER2D printer2d(world);
  PRINTER3D printer3d(world);
  
  string string2d = "newGrid2d";
  string string3d = "newGrid3d";
  
  printer2d.print(string2d,0,newMesh2d);
  printer3d.print(string3d,0,newMesh3d);
}
