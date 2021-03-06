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
#include "search2dA.hpp"

#include "feQr2d.hpp"
#include "feStaticField2d.hpp"

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
  
  typedef linearQuad     GEOSHAPE;
  typedef pMapItemShare  PMAPTYPE;  
  
  string meshFile  = "../tests/morganaMeshes/mignonQuad2dB.unv";
  string colorFile = "../tests/morganaMeshes/mignonQuad2dB_color.unv";
  
  meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.femap_to_stdB(meshFile, colorFile, false);
  
  
  //Download data
  typedef point3d            DOFTYPE;
  typedef feQr2d<1,PMAPTYPE> FETYPE;
  
  typedef meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH2D     MESH2D;
  typedef meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT2D  CONNECT2D;
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  
  //Set the field
  dofMapStatic2d_options mapOptions;
  mapOptions.addGeoId(1);
  mapOptions.addGeoId(2);
  
  feStaticField2d<FETYPE,DOFTYPE,dms2d_vectMajor,dms2d_geoIdMode> field2d;
  field2d.setCommunicator(world);
  field2d.setGeometry(grid2d,connectGrid2d);
  field2d.setOptions(mapOptions);
  
  field2d.startup();
  
  for(UInt i=1; i <= field2d.size(); ++i)
  {
    field2d.setDofL(i, grid2d->getNodeL(i) );
  }
  
  
  //Search 3d Init
  typedef search2dA<GEOSHAPE,PMAPTYPE,PMAPTYPE>   SEARCH2D;
  typedef searchData<PMAPTYPE>                  SEARCHDATA;
  
  SEARCH2D search(world);
  search.setMesh(grid2d,connectGrid2d);
  
  search.localInit(3.0);
  search.globalInit();
  
  
  //Global search point
  pVect<point3d,PMAPTYPE>    Pvect;
  pVect<SEARCHDATA,PMAPTYPE> dataVect;
  
  point3d P,V,Y;
  
  P.setX(1.3); P.setY(0.3); P.setZ(0.0);
  Pvect.push_back(P,PMAPTYPE(1,1,pid,false,false));
  Pvect.updateFinder();
  
  search.findGlobal(Pvect,dataVect);
  
  UInt elL = dataVect(1).getElMap().getLid();
         Y = dataVect(1).getLocCoord();
  
  
  //Evaluate
  if(pid == dataVect(1).getElMap().getPid())
  {
    field2d.evalL(elL,Y,V);
    cout << "pid " << pid << " eval " << V;
  }
}
