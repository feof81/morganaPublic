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
#include "ggmGridSearch.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  //Communicator-----------------------------------------------------
  environment  env(argc,argv);
  communicator world;
  
  //Typedefs---------------------------------------------------------
  typedef linearTetra    GEOSHAPE3D;
  typedef linearTriangle GEOSHAPE2D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  typedef mesh3d<GEOSHAPE3D,ELMAP,NODEMAP>     MESH3D;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>     MESH2D;
  typedef connect3d<GEOSHAPE3D,ELMAP,NODEMAP>  CONNECT3D;
  typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>  CONNECT2D;
  typedef meshInit3d<GEOSHAPE3D,ELMAP,NODEMAP> MESHINIT;
  typedef ggmGridSearch<MESH2D,MESH3D>         GRIDSEARCH;
  typedef GRIDSEARCH::OUTVECT OUTVECT;
  
  //Grids------------------------------------------------------------
  string meshFile = "./tests/morganaMeshes/mignon3dC.msh";
  
  MESHINIT init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  //Gen input vector-------------------------------------------------  
  sVect<point3d> inLocCoords;
  inLocCoords.push_back(point3d(1.0/3.0, 1.0/3.0, 0.0));
  
  //Test-------------------------------------------------------------
  GRIDSEARCH gridMatcher(world);
  gridMatcher.setMesh(grid2d,
                      connectGrid2d,
                      grid3d,
                      connectGrid3d,
                      grid2d,
                      connectGrid2d);
 
  OUTVECT outVect = gridMatcher.search(inLocCoords);
  
  //Output-----------------------------------------------------------
  if(world.rank() == 0)
  { cout << outVect << endl;}
  
  world.barrier();
  
  if(world.rank() == 1)
  { cout << outVect << endl;}
  
  world.barrier();
}