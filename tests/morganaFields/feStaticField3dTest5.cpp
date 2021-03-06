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
  
  typedef linearTriangle GEOSHAPE2D;
  typedef linearTetra    GEOSHAPE3D;
  typedef pMapItemShare  PMAPTYPE;
  
  string meshFile = "./tests/morganaMeshes/mignon3dC.msh";
  
  meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Download data
  typedef komplex            DOFTYPE;
  typedef fePr3d<1,PMAPTYPE> FETYPE;
  
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>::MESH2D     MESH2D;
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>::MESH3D     MESH3D;
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>::CONNECT2D  CONNECT2D;
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>::CONNECT3D  CONNECT3D;
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
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
  
  point3d P;
  
  for(UInt i=1; i <= field3d.size(); ++i)
  {
    P = grid3d->getNodeL(i);
    field3d.setDofL(i, komplex(P.getX(), P.getY()) );
  }
  
  //Search 3d Init
  typedef search3dA<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>   SEARCH3D;
  typedef searchData<PMAPTYPE>                   SEARCHDATA;
  
  SEARCH3D search(world);
  search.setMesh(grid2d,grid3d,connectGrid2d,connectGrid3d);
  
  search.localInit(3.0);
  search.globalInit();
  
  //Global search point
  pVect<point3d,PMAPTYPE>    Pvect;
  pVect<SEARCHDATA,PMAPTYPE> dataVect;
  
  point3d Y;
  komplex V;
  
  P.setX(0.3); P.setY(0.4); P.setZ(0.0);
  Pvect.push_back(P,PMAPTYPE(1,1,pid,false,false));
  Pvect.updateFinder();
  
  search.findGlobal(Pvect,dataVect);
  
  UInt elL = dataVect(1).getElMap().getLid();
         Y = dataVect(1).getLocCoord();
  
  //Evaluate
  if(pid == dataVect(1).getElMap().getPid())
  {
    field3d.evalL(elL,Y,V);
    cout << "pid " << pid << " eval " << V;
  }
}

