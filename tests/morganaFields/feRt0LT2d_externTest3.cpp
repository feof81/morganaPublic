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
#include "meshInit3d.hpp"

#include "fePr2d.hpp"
#include "fePr3d.hpp"
#include "feRt0LT3d.hpp"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"
#include "feStaticField2d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "feRt0LT2d_extern.hpp"


using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  //Comm-------------------------------------------------------------
  environment  env(argc,argv);
  communicator world;
  
  //Typedefs---------------------------------------------------------
  typedef pMapItemShare        PMAPTYPE;
  typedef feRt0LT2d<PMAPTYPE>  FIELD_FETYPE;
  typedef Real                 FIELD_DOFTYPE;

  typedef linearTetra                      GEOSHAPE3D;
  typedef typename FIELD_FETYPE::GEOSHAPE  GEOSHAPE2D;
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> MESHINIT;
 
  //Loading----------------------------------------------------------
  string meshFile = "./tests/morganaMeshes/mignon3dM.msh";
  
  MESHINIT init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Download data----------------------------------------------------
  typedef typename MESHINIT::MESH2D     MESH2D;
  typedef typename MESHINIT::CONNECT2D  CONNECT2D;

  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  //Testing----------------------------------------------------------
  typedef feRt0LT2d_extern<PMAPTYPE>  EXTERNAL;
  typedef typename EXTERNAL::FECARDS  FECARDS;
  
  EXTERNAL external(grid2d,connectGrid2d);
  external.setCommDev(world);
  
  std::set<UInt> activeGeoIds;
  activeGeoIds.insert(3);
  activeGeoIds.insert(4);
  activeGeoIds.insert(5);
  
  std::set<UInt> majorGeoIds;
  majorGeoIds.insert(5);
  
  if(world.rank() == 0)
  { cout << grid2d->getElements() << endl; }
  
  if(world.rank() == 0)
  { cout << grid2d->getEdges() << endl; }
  
  FECARDS cards = external.buildFeCards(activeGeoIds,majorGeoIds);
  
  if(world.rank() == 0)
  { cout << cards << endl; }
}
