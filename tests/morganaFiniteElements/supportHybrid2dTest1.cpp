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

#include "supportHybrid2d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef linearTriangle                         GEOSHAPE;
  typedef pMapItem                               PMAPTYPE;
  typedef meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> INIT;

  //Loading
  string meshFile = "./tests/morganaMeshes/mignon2dB.msh";
  
  INIT init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Typedefs
  typedef typename INIT::MESH2D    MESH2D;
  typedef typename INIT::CONNECT2D CONNECT2D;
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  //Test support
  typedef supportHybrid2d<GEOSHAPE,PMAPTYPE> SUPPORT;
  
  SUPPORT support(grid2d,connectGrid2d);  
  support.setElement2d(2);
  support.setLocalEdge(2);
  
  point3d Yd(0.5, 0.0, 0.0);
  
  cout << "Yv : " << support.mapVolumeY(Yd);
  cout << "N  : " << support.computeNormal(Yd);
  cout << "Edge points : " << endl << support.getGlobEdgePoints() << endl;
}
