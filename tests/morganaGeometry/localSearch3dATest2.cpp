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

#include "connect2d.hpp"
#include "connect3d.hpp"
#include "localSearch3dA.hpp"

#include "printMesh.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 1);
  
  typedef linearTetra    GEOSHAPE3D;
  typedef linearTriangle GEOSHAPE2D;
  typedef linearLine     GEOSHAPE1D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  typedef mesh3d<GEOSHAPE3D,ELMAP,NODEMAP>         MESH3D;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>         MESH2D;
  typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>         MESH1D;
  typedef connect3d<GEOSHAPE3D,ELMAP,NODEMAP>      CONNECT3D;
  typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>      CONNECT2D;
  typedef localSearch3dA<GEOSHAPE3D,ELMAP,NODEMAP> SEARCH3D;
  
  
  //Mesh Init 
  string meshFile = "./geometries/cubes3d/testCubeE.msh";
  
  meshInit3d<GEOSHAPE3D,ELMAP,NODEMAP> init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
   
  
  //Search 3d
  time_t start, end;
  Real totDist  = 0.0;
  UInt steps = 100000;
  searchData<ELMAP> outData;
  point3d P, Pmin(0.0, 0.5, -1.0), Pmax(1.0, 0.5, -1.0);
  
  
  SEARCH3D search;
  search.setMesh(grid2d,
                 grid3d,
                 connectGrid2d,
                 connectGrid3d);
  
  cout << "Local search init" << " "; time(&start); cout << endl;
  search.localInit(3.0);
  time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;
  
  
  
  cout << "Search" << " "; time(&start); cout << endl;
  geoMapInterface<GEOSHAPE3D> globMapper;
  sVect<point3d> points;
  point3d Ps;
  
  for(UInt i=0; i <= steps; ++i)
  {
    P = Pmin * (Real(i) / steps) + Pmax * (1.0 - Real(i) / steps);
    outData = search.findLocal(P);
    
    points = grid3d->getElementNodesL(outData.getElMap().getLid());
    Ps     = globMapper.getPosition(points,outData.getLocCoord());
    
    P.setZ(0.0);
    totDist += point3d::norm2(P-Ps);
  }
  
  time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;
  
  cout << "Last point     : " << Ps;
  cout << "Total residual : " << (totDist < 1e-10) << endl;
  cout << "Mean elements  : " << search.numFeasibElements / search.localQuery << endl;
}
