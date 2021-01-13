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

#include "printMesh.hpp"
#include "elCardFeeder1d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 1);
  
  typedef linearQuad     GEOSHAPE2D;
  typedef linearLine     GEOSHAPE1D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>     MESH2D;
  typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>     MESH1D;
  typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>  CONNECT2D;
  typedef connect1d<GEOSHAPE1D,ELMAP,NODEMAP>  CONNECT1D;
  
  
  string meshFile  = "./tests/morganaMeshes/mignonQuad2dA.unv";
  string colorFile = "./tests/morganaMeshes/mignonQuad2dA_color.unv";
  
  meshInit2d<GEOSHAPE2D,ELMAP,NODEMAP> init(world);
  init.femap_to_stdB(meshFile, colorFile, false);
  
  
  //Extract data
  RCP<MESH1D>           grid1d = init.getGrid1d();
  RCP<CONNECT1D> connectGrid1d = init.getConnectGrid1d();
  
  
  //Test the feeder
  UInt el = 3;
  
  elCardFeeder1d<GEOSHAPE1D,ELMAP> elFeeder(grid1d,connectGrid1d);
  elCard1d<GEOSHAPE1D,ELMAP> elCard = elFeeder.getCardLocal(el);
  
  
  if(world.rank() == 0)
  {
    //cout << grid1d->getElements() << endl;
    cout << elCard << endl;
  }
}
