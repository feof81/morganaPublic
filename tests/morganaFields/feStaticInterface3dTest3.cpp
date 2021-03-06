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

#include "fePr3d.hpp"
#include "feStaticField3d.hpp"

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
  
  typedef linearTetra    GEOSHAPE;
  typedef pMapItemShare  PMAPTYPE;
  
  string meshFile = "../tests/morganaMeshes/mignon3dC.msh";
  
  meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  
  //Download data
  typedef point3d            DOFTYPE;
  typedef fePr3d<1,PMAPTYPE> FETYPE;
  
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH3D     MESH3D;
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT3D  CONNECT3D;
  
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  
  //Set the field
  dofMapStatic3d_options mapOptions;
  mapOptions.addGeoId(1);
  mapOptions.addGeoId(2);
  
  feStaticField3d<FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_geoIdMode> field3d;
  field3d.setCommunicator(world);
  field3d.setGeometry(grid3d,connectGrid3d);
  field3d.setOptions(mapOptions);
  
  field3d.startup();
  
  for(UInt i=1; i <= field3d.size(); ++i)
  {
    field3d.setDofL(i, grid3d->getNodeL(i) );
  }
  
  
  //Static interface
  feStaticInterface3d<FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_geoIdMode> interfaceClone3d;
  interfaceClone3d.clone(field3d.getDofMapper());
  
  feStaticInterface3d<FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_geoIdMode> interface3d(world,*grid3d,*connectGrid3d);
  interface3d.setOptions(mapOptions);
  interface3d.startup();
  
  
  if(world.rank() == 0)
  {
    UInt numBasis = interface3d.getNumBasis();
    sVect<UInt> indices(numBasis);
    sVect<point3d> val(numBasis);
    
    UInt elL = 1;
    UInt I=1, J=1;
    point3d Y(0.1,0.1,0.1);
    
    interface3d.feIndicesLL(elL,I,J,indices);
    interface3d.feEvalL(elL,I,J,Y,val);
    
    cout << "Original Indices: " << endl;
    cout << indices << endl;
    
    cout << "Original Eval: " << endl;
    cout << val << endl;
    
    interfaceClone3d.feIndicesLL(elL,I,J,indices);
    interfaceClone3d.feEvalL(elL,I,J,Y,val);
    
    cout << "Clone Indices: " << endl;
    cout << indices << endl;
    
    cout << "Clone Eval: " << endl;
    cout << val << endl;
  }
}

