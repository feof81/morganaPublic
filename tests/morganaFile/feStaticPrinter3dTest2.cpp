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

#include "feQr3d.hpp"
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
  
  UInt pid = world.rank();
  
  //Mesh Init
  typedef pMapItemShare      PMAPTYPE;
  typedef point3d            DOFTYPE;
  typedef feQr3d<1,PMAPTYPE> FETYPE;
  typedef typename FETYPE::GEOSHAPE GEOSHAPE;
  
  string meshFile  = "./tests/morganaMeshes/mignonHexaB.unv";
  string colorFile = "./tests/morganaMeshes/mignonHexaB_color.unv";
  
  meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.femap_to_stdA(meshFile, colorFile);
  
  
  //Download data
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH3D     MESH3D;
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT3D  CONNECT3D;
  
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  
  //Set the field
  dofMapStatic3d_options mapOptions;
  mapOptions.addGeoId(1);
  mapOptions.addGeoId(2);
  
  typedef feStaticField3d<FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_geoIdMode> FIELD;
  
  FIELD field3d;
  field3d.setCommunicator(world);
  field3d.setGeometry(grid3d,connectGrid3d);
  field3d.setOptions(mapOptions);
  field3d.startup();
  
  
  //Field initialization
  feStaticField3dGlobalManip<FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_geoIdMode> manipulator(world);
  
  sArray<string> functions(3,1);
  functions(1,1) = "x";
  functions(2,1) = "y";
  functions(3,1) = "z";
  
  manipulator.initilize(functions,field3d);
  
    
  //Mesh printing
  feStaticFieldPrinter3d<FIELD> printer(world);
  printer.printHDF5_nodes("testField",pid,field3d);
}
