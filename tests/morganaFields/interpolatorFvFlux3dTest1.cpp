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
#include "meshInit3d.hpp"

#include "feFvFlux3d.hpp"
#include "fePr3d.hpp"

#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "interpolatorFvFlux3d.hpp"
#include "traitsFvEl3d.hpp"


// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  UInt pid = world.rank();
  assert(world.size() == 2);
  
  
  //Mesh Init--------------------------------------------------------------------------------------
  typedef pMapItemShare  PMAPTYPE;
  typedef point3d        DOFTYPE;
  
  typedef fePr3d<1,PMAPTYPE>               SOURCE_FETYPE;
  typedef typename SOURCE_FETYPE::GEOSHAPE SOURCE_GEOSHAPE;
  
  string sourceMeshFile = "../tests/morganaMeshes/mignon3dB.msh";
  
  meshInit3d<SOURCE_GEOSHAPE,PMAPTYPE,PMAPTYPE> sourceInit(world);
  sourceInit.gmMesh_to_stdB(sourceMeshFile, false);

  
  //Download data----------------------------------------------------------------------------------
  typedef meshInit3d<SOURCE_GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH3D     MESH3D;
  typedef meshInit3d<SOURCE_GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT3D  CONNECT3D;
  
  RCP<MESH3D>           grid3d = sourceInit.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = sourceInit.getConnectGrid3d();
  
  
  //Set the source field---------------------------------------------------------------------------
  typedef feStaticField3d<SOURCE_FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_geoIdMode> SOURCE_FIELD;
  
  dofMapStatic3d_options sourceOptions;
  sourceOptions.addGeoId(1);
  //sourceOptions.addGeoId(2);
  
  SOURCE_FIELD sourceField;
  sourceField.setCommunicator(world);
  sourceField.setGeometry(grid3d,connectGrid3d);
  sourceField.setOptions(sourceOptions);
  sourceField.startup();
  
  
  //Source field sourceInitialization--------------------------------------------------------------------
  feStaticField3dGlobalManip<SOURCE_FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_geoIdMode> manipulator(world);
  
  sArray<string> functions(3,1);
  functions(1,1) = "1";
  functions(2,1) = "0";
  functions(3,1) = "0";
  
  manipulator.initilize(functions,sourceField);
  
  
  //Target field-----------------------------------------------------------------------------------
  typedef interpolatorFvFlux3d<SOURCE_FIELD> INTERPOLATOR;
  typedef INTERPOLATOR::TARGETFIELD TARGET_FIELD;
  
  dofMapStatic3d_options targetOptions;
  targetOptions.addGeoId(1);
  //targetOptions.addGeoId(2);
  
  TARGET_FIELD targetField;
  targetField.setCommunicator(world);
  targetField.setGeometry(grid3d,connectGrid3d);
  targetField.setOptions(targetOptions);
  targetField.startup();
  
  
  //Interpolator-----------------------------------------------------------------------------------
  INTERPOLATOR  interpolator(world);
  interpolator.setCommDev(world);
  interpolator.interpolate(sourceField,targetField);
  
  if(world.rank() == 0)
  {
    cout << targetField.getDofVect() << endl;
  }
}
