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
  
  //Mesh Init
  typedef pMapItem           PMAPTYPE;
  typedef point3d            DOFTYPE;
  typedef feQr3d<1,PMAPTYPE> FETYPE;
  typedef typename FETYPE::GEOSHAPE GEOSHAPE;
  
  string meshFile  = "../tests/morganaMeshes/mignonHexaB.unv";
  string colorFile = "../tests/morganaMeshes/mignonHexaB_color.unv";
  
  meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.femap_to_stdA(meshFile, colorFile, false);
  
  
  //Download data
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH3D     MESH3D;
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT3D  CONNECT3D;
  
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  
  //Set the field
  dofMapStatic3d_options mapOptions;
  mapOptions.addGeoId(1);
  
  feStaticField3d<FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_geoIdMode> field3d;
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
  
    
  //Mesh gathering
  RCP<MESH3D>  newGrid3d = rcp(new MESH3D);
  mesh3dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> gathering(world);
  gathering.gather(1,newGrid3d,grid3d);
  
  RCP<CONNECT3D> newConnectGrid3d = rcp(new CONNECT3D(world));
  newConnectGrid3d->setMesh3d(newGrid3d);
  newConnectGrid3d->buildConnectivity();  
  
  
  //Change the grid3d 
  manipulator.changeGrid(*newGrid3d, *newConnectGrid3d, mapOptions, field3d);
  
  
  //Printout of the dofs
  if(pid == 1)
  {
    cout << "Dofs" << endl;
    cout << field3d.getDofVect() << endl;
    
  }
}
