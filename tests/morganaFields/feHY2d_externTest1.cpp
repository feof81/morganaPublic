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

#include "fePr2d.hpp"
#include "feRt0LT2d.hpp"
#include "feRt0LT2d_extern.hpp"

#include "elCardFeeder2d.hpp"
#include "feHY2d_extern.hpp"


using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef pMapItemShare        PMAPTYPE;
  typedef feRt0LT2d<PMAPTYPE>  FIELD_FETYPE;
  typedef Real                 FIELD_DOFTYPE;
  typedef fePr2d<0,PMAPTYPE>   TEST_FETYPE;
  typedef Real                 TEST_DOFTYPE;

  typedef typename FIELD_FETYPE::GEOSHAPE          GEOSHAPE2D;
  typedef meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> MESHINIT;
  typedef typename MESHINIT::MESH2D                MESH2D;
  typedef typename MESHINIT::CONNECT2D             CONNECT2D;

  assert(world.size() <= 2);
  
  //Loading
  string meshFile = "../tests/morganaMeshes/mignon2dB.msh";
  
  MESHINIT init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Download data
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  //Testing extern
  typedef feHY2d_extern<PMAPTYPE>  EXTERN;
  typedef typename EXTERN::FECARDS FECARDS;
  typedef elCardFeeder2d<GEOSHAPE2D,PMAPTYPE> ELCARDFEEDER;
  
  feHY2d_extern<PMAPTYPE> feeder(grid2d,connectGrid2d);
  feeder.setCommDev(world);
  
  FECARDS feCards = feeder.buildFeCards();
  
  if(world.rank() == 0)
  {
    cout << grid2d->getElements() << endl;
    cout << feCards << endl;
  }
}

