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

#include "fePr3d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"
#include "feStaticFieldPrinter3d.hpp"

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
  
  typedef fePr3d<0,PMAPTYPE>  FETYPE;
  typedef Real                DOFTYPE;
  typedef feStaticField3d<FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_geoIdMode>            FIELD;
  typedef feStaticField3dGlobalManip<FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_geoIdMode> FIELD_MANIP;
  
  typedef FIELD::OPTIONS FIELD_OPTIONS;
  
  
  //Mesh loading-----------------------------------------------------
  string meshFile = "./geometries/biCube3d/biCubeC.msh";
  //string meshFile = "./geometries/cubes3d/testCubeB.msh";
  
  INIT init(world);
  init.gmMesh_to_stdA(meshFile);
  
  //Meshes links-----------------------------------------------------
  Teuchos::RCP<MESH3D>       oldGrid3d = init.getGrid3d();
  Teuchos::RCP<MESH2D>       oldGrid2d = init.getGrid2d();
  Teuchos::RCP<CONNECT3D> oldConnect3d = init.getConnectGrid3d();
  Teuchos::RCP<CONNECT2D> oldConnect2d = init.getConnectGrid2d();
  
  //Init the field---------------------------------------------------
  FIELD_OPTIONS fieldOptions;
  fieldOptions.addGeoId(1);
  
  FIELD field;
  field.setCommunicator(world);
  field.setGeometry(oldGrid3d,oldConnect3d);
  field.setOptions(fieldOptions);
  field.startup();
  
  sArray<string> exprString(1,1);
  exprString(1,1) = "x";
  
  FIELD_MANIP fieldManip(world);
  fieldManip.initilize(exprString,field);
  
  //Mesh rebalance---------------------------------------------------
  set<UInt> activeGeoIds;
  activeGeoIds.insert(1);
  //activeGeoIds.insert(2);
  
  init.reinit_stdA(oldGrid2d,
                   oldGrid3d,
                   activeGeoIds);

  Teuchos::RCP<MESH3D>       newGrid3d = init.getGrid3d();
  Teuchos::RCP<MESH2D>       newGrid2d = init.getGrid2d();
  Teuchos::RCP<CONNECT3D> newConnect3d = init.getConnectGrid3d();
  Teuchos::RCP<CONNECT2D> newConnect2d = init.getConnectGrid2d();
  
  //Field rebalance--------------------------------------------------
  fieldManip.changeGrid(*newGrid3d,*newConnect3d,fieldOptions,field);
  
  //Field printer----------------------------------------------------
  feStaticFieldPrinter3d<FIELD> printer(world);
  string flag = "FieldBalanced";
  printer.printHDF5_elements(flag,0,field);
}
