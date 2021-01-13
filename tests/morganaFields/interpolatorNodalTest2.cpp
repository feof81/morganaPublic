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
  typedef fePr3d<1,PMAPTYPE> TARGET_FETYPE;
  typedef typename SOURCE_FETYPE::GEOSHAPE SOURCE_GEOSHAPE;
  typedef typename TARGET_FETYPE::GEOSHAPE TARGET_GEOSHAPE;
  
  string sourceMeshFile = "./geometries/cadInterface/sphereCylinder05_hole.msh";
  string targetMeshFile = "./geometries/cadInterface/sphereCylinder05.msh";
  
  meshInit3d<SOURCE_GEOSHAPE,PMAPTYPE,PMAPTYPE> sourceInit(world);
  sourceInit.gmMesh_to_stdB(sourceMeshFile);
  
  meshInit3d<TARGET_GEOSHAPE,PMAPTYPE,PMAPTYPE> targetInit(world);
  targetInit.gmMesh_to_stdB(targetMeshFile);
  
  //Download data----------------------------------------------------------------------------------
  typedef meshInit3d<SOURCE_GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH3D     MESH3D;
  typedef meshInit3d<SOURCE_GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT3D  CONNECT3D;
  
  typedef meshInit3d<SOURCE_GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH2D     MESH2D;
  typedef meshInit3d<SOURCE_GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT2D  CONNECT2D;
  
  RCP<MESH2D>       sourceGrid2d = sourceInit.getGrid2d();
  RCP<MESH3D>       sourceGrid3d = sourceInit.getGrid3d();
  RCP<CONNECT2D> sourceConnect2d = sourceInit.getConnectGrid2d();
  RCP<CONNECT3D> sourceConnect3d = sourceInit.getConnectGrid3d();
  
  RCP<MESH2D>       targetGrid2d = targetInit.getGrid2d();
  RCP<MESH3D>       targetGrid3d = targetInit.getGrid3d();
  RCP<CONNECT2D> targetConnect2d = targetInit.getConnectGrid2d();
  RCP<CONNECT3D> targetConnect3d = targetInit.getConnectGrid3d();
  
  //Set the source field---------------------------------------------------------------------------
  typedef feStaticField3d<SOURCE_FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_allMode> SOURCE_FIELD;
  
  dofMapStatic3d_options sourceOptions;
  
  SOURCE_FIELD sourceField;
  sourceField.setCommunicator(world);
  sourceField.setGeometry(sourceGrid3d,sourceConnect3d);
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
  typedef feStaticField3d<TARGET_FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_allMode> TARGET_FIELD;
  
  dofMapStatic3d_options targetOptions;
  
  TARGET_FIELD targetField;
  targetField.setCommunicator(world);
  targetField.setGeometry(targetGrid3d,targetConnect3d);
  targetField.setOptions(targetOptions);
  targetField.startup();
  
  //Interpolator-----------------------------------------------------------------------------------
  time_t start, end;
  
  interpolatorNodal<SOURCE_FIELD,TARGET_FIELD> interpolator(world);
  interpolator.setMesh(sourceGrid2d,
                       sourceGrid3d,
                       sourceConnect2d,
                       sourceConnect3d);
  
  if(world.rank() == 0) {cout << "Local Init" << " "; time(&start); cout << endl;}
  interpolator.localInit(2.0);
  if(world.rank() == 0) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  if(world.rank() == 0) {cout << "Global Init" << " "; time(&start); cout << endl;}
  interpolator.globalInit();
  if(world.rank() == 0) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  if(world.rank() == 0) {cout << "Find Dofs" << " "; time(&start); cout << endl;}
  interpolator.findDofs(targetField);
  if(world.rank() == 0) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}

  if(world.rank() == 0) {cout << "Exchange Data" << " "; time(&start); cout << endl;}
  interpolator.exchangeData(sourceField,targetField);
  if(world.rank() == 0) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
    
  //Source field printing
  feStaticFieldPrinter3d<SOURCE_FIELD> sourcePrinter(world);
  sourcePrinter.printHDF5_nodes("sourceField",pid,sourceField);
  
  //Target field printing
  feStaticFieldPrinter3d<TARGET_FIELD> targetPrinter(world);
  targetPrinter.printHDF5_nodes("targetField",pid,targetField);
}
