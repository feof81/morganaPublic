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


/* Grid non nested obtained through a shift x+10 */
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
  typedef GRIDSEARCH::OUTDATA OUTDATA;
  
  //Grids------------------------------------------------------------
  string meshFile = "./geometries/cubes3d/testCubeB.msh";
  
  MESHINIT init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH3D>       srcGrid3d = init.getGrid3d();
  RCP<MESH2D>       srcGrid2d = init.getGrid2d();
  RCP<CONNECT3D> srcConnect3d = init.getConnectGrid3d();
  RCP<CONNECT2D> srcConnect2d = init.getConnectGrid2d();
  
  //Map the faces----------------------------------------------------
  RCP<MESH2D> tgtGrid2d(new MESH2D(srcGrid3d->getNodes(),
                                   srcGrid3d->getIsVertex(),
                                   srcGrid3d->getFaces()));

  RCP<CONNECT2D> tgtConnect2d(new CONNECT2D(world));
  tgtConnect2d->setMesh2d(tgtGrid2d);
  
  //Map search-------------------------------------------------------
  sVect<point3d> inLocCoords;
  
  GRIDSEARCH gridMatcher(world);
  gridMatcher.setMesh(tgtGrid2d,
                      tgtConnect2d,
                      srcGrid3d,
                      srcConnect3d,
                      srcGrid2d,
                      srcConnect2d);
 
  OUTVECT outVect = gridMatcher.search(inLocCoords);
  
  //Map correspondence-----------------------------------------------
  typedef pMap<ELMAP> FACEMAP;
  FACEMAP oldMap = outVect.getMapRef();
  FACEMAP newMap(outVect.size());
  
  for(UInt i=1; i <= outVect.size(); ++i)
  {
    assert(outVect(i).getIsNested());
    newMap(i) = outVect.getMapL(i);
  }
  
  for(UInt k=0; k < world.size(); ++k)
  {
    if(world.rank() == k)
    {
      for(UInt i=1; i <= outVect.size(); ++i)
      {
        cout << "oldMap : " << oldMap(i);
        cout << "newMap : " << newMap(i) << endl;
      }
    }
    
    world.barrier();
  }
}