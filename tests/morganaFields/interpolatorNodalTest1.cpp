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
#include "meshInit3d.hpp"

#include "fePr2d.hpp"
#include "fePr3d.hpp"

#include "feStaticField2d.hpp"
#include "feStaticField3d.hpp"

#include "feStaticField2dGlobalManip.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "feStaticFieldPrinter2d.hpp"
#include "feStaticFieldPrinter3d.hpp"

#include "interpolatorNodal.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  UInt pid = world.rank();
  
  
  //Mesh Init--------------------------------------------------------------------------------------
  typedef pMapItemShare      PMAPTYPE;
  typedef point3d            DOFTYPE;
  typedef fePr3d<1,PMAPTYPE> SOURCE_FETYPE;
  typedef fePr2d<1,PMAPTYPE> TARGET_FETYPE;
  typedef typename SOURCE_FETYPE::GEOSHAPE SOURCE_GEOSHAPE;
  typedef typename TARGET_FETYPE::GEOSHAPE TARGET_GEOSHAPE;
  
  string sourceMeshFile = "./geometries/cubes3d/testCubeA.msh";
  string targetMeshFile = "./geometries/rectangles2d/rectangleAbis.msh";
  
  meshInit3d<SOURCE_GEOSHAPE,PMAPTYPE,PMAPTYPE> sourceInit(world);
  sourceInit.gmMesh_to_stdB(sourceMeshFile, false);
  
  meshInit2d<TARGET_GEOSHAPE,PMAPTYPE,PMAPTYPE> targetInit(world);
  targetInit.gmMesh_to_stdB(targetMeshFile, false);
  
  
  //Download data----------------------------------------------------------------------------------
  typedef meshInit3d<SOURCE_GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH3D     MESH3D;
  typedef meshInit3d<SOURCE_GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT3D  CONNECT3D;
  
  typedef meshInit2d<TARGET_GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH2D     MESH2D;
  typedef meshInit2d<TARGET_GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT2D  CONNECT2D;
  
  RCP<MESH3D>           grid3d = sourceInit.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = sourceInit.getConnectGrid3d();
  
  RCP<MESH2D>           grid2d = targetInit.getGrid2d();
  RCP<CONNECT2D> connectGrid2d = targetInit.getConnectGrid2d();
  
  
  //Set the source field---------------------------------------------------------------------------
  typedef feStaticField3d<SOURCE_FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_allMode> SOURCE_FIELD;
  
  dofMapStatic3d_options sourceOptions;
  
  SOURCE_FIELD sourceField;
  sourceField.setCommunicator(world);
  sourceField.setGeometry(grid3d,connectGrid3d);
  sourceField.setOptions(sourceOptions);
  sourceField.startup();
  
  
  //Source field sourceInitialization--------------------------------------------------------------------
  feStaticField3dGlobalManip<SOURCE_FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_allMode> manipulator(world);
  
  sArray<string> functions(3,1);
  functions(1,1) = "x";
  functions(2,1) = "y";
  functions(3,1) = "z";
  
  manipulator.initilize(functions,sourceField);
  
  
  //Target field-----------------------------------------------------------------------------------
  typedef feStaticField2d<TARGET_FETYPE,DOFTYPE,dms2d_vectMajor,dms2d_allMode> TARGET_FIELD;
  
  dofMapStatic2d_options targetOptions;
  
  TARGET_FIELD targetField;
  targetField.setCommunicator(world);
  targetField.setGeometry(grid2d,connectGrid2d);
  targetField.setOptions(targetOptions);
  targetField.startup();
  
  
  //Interpolator-----------------------------------------------------------------------------------
  interpolatorNodal<SOURCE_FIELD,TARGET_FIELD> interpolator(world);
  interpolator.setMesh(grid2d,grid3d,connectGrid2d,connectGrid3d);
  interpolator.localInit(2.0);
  interpolator.globalInit();
  interpolator.findDofs(targetField);
  interpolator.exchangeData(sourceField,targetField);
    
  //Source field printing
  feStaticFieldPrinter3d<SOURCE_FIELD> sourcePrinter(world);
  sourcePrinter.printHDF5_nodes("sourceField",pid,sourceField);
  
  //Target field printing
  feStaticFieldPrinter2d<TARGET_FIELD> targetPrinter(world);
  targetPrinter.printHDF5_nodes("targetField",pid,targetField);
}
