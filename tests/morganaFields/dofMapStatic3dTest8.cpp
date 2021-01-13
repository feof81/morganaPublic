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
#include "dofMapStatic3d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef linearTetra    GEOSHAPE;
  typedef pMapItemShare  PMAPTYPE;
  
  string meshFile = "../tests/morganaMeshes/mignon3dC.msh";
  //string meshFile = "./geometries/cubes3d/testCubeD.msh";
  
  meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile,false);
  
  
  //Map test 
  typedef point3d            DOFTYPE;
  typedef fePr3d<2,PMAPTYPE> FETYPE;
  
  dofMapStatic3d_options mapOptions;
  mapOptions.addGeoId(1);
  
  dofMapStatic3d<FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_allMode> dofMap, dofMapCopy;
  
  dofMap.setCommunicator(world);
  dofMap.setGeometry(init.getGrid3d(), init.getConnectGrid3d());
  dofMap.setOptions(mapOptions);
  dofMap.startup();
  
  dofMapCopy = dofMap;
  
  feStaticDofCard3d card;
  card.setGeoType(VERTEX);
  card.setLevel(1);
  card.setLocalId(4);
  card.setLocalElId(1);
  
  UInt I = 1, J = 1;
  
  
  world.barrier();
  if(world.rank() == 1)
  {
    cout << "Dof Map " << endl;
    cout << dofMap.getDofMap() << endl << endl;
    
    cout << "Dof Map Copy " << endl;
    cout << dofMapCopy.getDofMap() << endl << endl;
    
    cout << "Original " << endl;
    cout << "dofL  : " << dofMap.mapDofL(card) << endl;
    cout << "dofG  : " << dofMap.mapDofG(card) << endl;
    cout << "listL : " << dofMap.mapListL(I,J,card) << endl;
    cout << "listG : " << dofMap.mapListG(I,J,card) << endl;
    
    cout << "Copy " << endl;
    cout << "dofL  : " << dofMapCopy.mapDofL(card) << endl;
    cout << "dofG  : " << dofMapCopy.mapDofG(card) << endl;
    cout << "listL : " << dofMapCopy.mapListL(I,J,card) << endl;
    cout << "listG : " << dofMapCopy.mapListG(I,J,card) << endl;
  }
}
