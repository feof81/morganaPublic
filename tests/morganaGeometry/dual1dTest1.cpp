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
#include "dual1d.hpp"

#include "printMesh.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef linearQuad     GEOSHAPE2D;
  typedef linearLine     GEOSHAPE1D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>  MESH2D;
  typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>  MESH1D;
  typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>  CONNECT2D;
  typedef connect1d<GEOSHAPE1D,ELMAP,NODEMAP>  CONNECT1D;
  
  
  string meshFile  = "./tests/morganaMeshes/mignonQuad2dA.unv";
  string colorFile = "./tests/morganaMeshes/mignonQuad2dA_color.unv";
  
  meshInit2d<GEOSHAPE2D,ELMAP,NODEMAP> init(world);
  init.femap_to_stdB(meshFile, colorFile);
  
  
  //Extract the meshes
  RCP<MESH1D>           grid1d = init.getGrid1d();
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT1D> connectGrid1d = init.getConnectGrid1d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  
  //Dual 1d
  dual1d<GEOSHAPE1D,ELMAP,NODEMAP> dualGrid1d(world);
  dualGrid1d.setMesh1d(grid1d);
  dualGrid1d.setConnect1d(connectGrid1d);
  
  dualGrid1d.buildDualMesh();
  dualGrid1d.buildFiniteVolumeData(grid2d,connectGrid2d);
}
