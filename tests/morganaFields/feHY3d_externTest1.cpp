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

#include "fePr2d.hpp"
#include "fePr3d.hpp"
#include "feRt0LT3d.hpp"
#include "feRt0LT3d_extern.hpp"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"
#include "feStaticField2d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "operatorEL3d.hpp"
#include "opEL3dB.hpp"

#include "localMatrixEL3d.hpp"
#include "matrixBuilder.hpp"

#include "feHY3d_extern.hpp"
#include "elCardFeeder3d.hpp"


using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef pMapItemShare        PMAPTYPE;
  typedef feRt0LT3d<PMAPTYPE>  FIELD_FETYPE;
  typedef Real                 FIELD_DOFTYPE;
  typedef fePr3d<0,PMAPTYPE>   TEST_FETYPE;
  typedef Real                 TEST_DOFTYPE;

  typedef typename FIELD_FETYPE::GEOSHAPE          GEOSHAPE3D;
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> MESHINIT;
  typedef typename MESHINIT::MESH2D                MESH2D;
  typedef typename MESHINIT::MESH3D                MESH3D;
  typedef typename MESHINIT::CONNECT2D             CONNECT2D;
  typedef typename MESHINIT::CONNECT3D             CONNECT3D;

  assert(world.size() == 2);
  
  //Loading
  string meshFile = "../tests/morganaMeshes/mignon3dB.msh";
  
  MESHINIT init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Download data
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //Testing extern
  typedef feHY3d_extern<PMAPTYPE>  EXTERN;
  typedef typename EXTERN::FECARDS FECARDS;
  typedef elCardFeeder3d<GEOSHAPE3D,PMAPTYPE> ELCARDFEEDER;
  
  feHY3d_extern<PMAPTYPE> feeder(grid3d,connectGrid3d,grid2d,connectGrid2d);
  feeder.setCommDev(world);
  
  set<UInt> geoIds3d, geoIds2d;
  geoIds2d.insert(1);
  geoIds2d.insert(2);
  geoIds3d.insert(1);
  geoIds3d.insert(2);
  
  FECARDS feCards = feeder.buildFeCards(geoIds3d,geoIds2d);
  
  if(world.rank() == 1)
  {
    cout << grid3d->getElements() << endl;
    cout << feCards << endl;
  }
}
