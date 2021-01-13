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
#include "search3dA.hpp"

#include "feOscCr3d.hpp"
#include "feOscCr3d_extern.h"
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
  
  Real K = 2.5;
  
  typedef linearTriangle GEOSHAPE2D;
  typedef linearTetra    GEOSHAPE3D;
  typedef pMapItemShare  PMAPTYPE;
  
  string meshFile = "./tests/morganaMeshes/mignon3dC.msh";
  
  meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Download data 
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>::MESH2D     MESH2D;
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>::MESH3D     MESH3D;
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>::CONNECT2D  CONNECT2D;
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>::CONNECT3D  CONNECT3D;
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //Build the feCards
  typedef feOscCr3d_extern<1,PMAPTYPE>::FECARDS FECARDS;
  
  feOscCr3d_extern<1,PMAPTYPE> feExtern;
  feExtern.setGeometry(grid3d,connectGrid3d);
  feExtern.setCommDev(world);
  FECARDS feCards = feExtern.buildFeCards(K);
  
  if(pid == 0)
  {
    for(UInt el=1; el <= grid2d->getNumElements(); ++el)
    {
      cout << "Nodes" << endl;
      cout << grid2d->getElementNodesL(el) << endl;
      
      cout << feCards(el) << endl;
    }
  }
}
