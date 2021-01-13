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
#include "search1dA.hpp"

#include "fePr1d.hpp"
#include "feStaticField1d.hpp"

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
  
  typedef linearLine      GEOSHAPE1D;
  typedef linearTriangle  GEOSHAPE2D;
  typedef pMapItemShare   PMAPTYPE;  
  
  string meshFile = "../tests/morganaMeshes/mignon2dB.msh";
  
  meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  
  //Download data
  typedef point3d            DOFTYPE;
  typedef fePr1d<1,PMAPTYPE> FETYPE;
  
  typedef meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>::MESH1D     MESH1D;
  typedef meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>::CONNECT1D  CONNECT1D;
  
  RCP<MESH1D>           grid1d = init.getGrid1d();
  RCP<CONNECT1D> connectGrid1d = init.getConnectGrid1d();
  
  
  //Set the field
  dofMapStatic1d_options mapOptions;
  mapOptions.addGeoId(1);
  mapOptions.addGeoId(2);
  
  feStaticField1d<FETYPE,DOFTYPE,dms1d_vectMajor,dms1d_geoIdMode> field1d;
  field1d.setCommunicator(world);
  field1d.setGeometry(grid1d,connectGrid1d);
  field1d.setOptions(mapOptions);
  
  field1d.startup();
  
  for(UInt i=1; i <= field1d.size(); ++i)
  {
    field1d.setDofL(i, grid1d->getNodeL(i) );
  }
  
  
  //Search 1d Init
  typedef search1dA<GEOSHAPE1D,PMAPTYPE,PMAPTYPE>  SEARCH1D;
  typedef searchData<PMAPTYPE>                   SEARCHDATA;
  
  SEARCH1D search(world);
  search.setMesh(grid1d,connectGrid1d);
  
  search.localInit(3.0);
  search.globalInit();
 
  
  //Global search point
  pVect<point3d,PMAPTYPE>    Pvect;
  pVect<SEARCHDATA,PMAPTYPE> dataVect;
  
  point3d P,V,Y;
  
  P.setX(1.6); P.setY(0.0); P.setZ(0.0);
  Pvect.push_back(P,PMAPTYPE(1,1,pid,false,false));
  Pvect.updateFinder();
  
  search.findGlobal(Pvect,dataVect);
  
  UInt elL = dataVect(1).getElMap().getLid();
         Y = dataVect(1).getLocCoord();
  
  
  //Evaluate
  if(pid == dataVect(1).getElMap().getPid())
  {
    field1d.evalL(elL,Y,V);
    cout << "pid " << pid << " eval " << V;
  }
}
