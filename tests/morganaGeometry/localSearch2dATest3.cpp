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

#include "connect1d.hpp"
#include "connect2d.hpp"
#include "localSearch2dA.hpp"

#include "printMesh.hpp"


using namespace std;
using namespace boost::mpi;
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef linearTriangle                           GEOSHAPE2D;
  typedef linearLine                               GEOSHAPE1D;
  typedef pMapItemShare                            ELMAP;
  typedef pMapItemShare                            NODEMAP;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>         MESH2D;
  typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>         MESH1D;
  typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>      CONNECT2D;
  typedef connect1d<GEOSHAPE1D,ELMAP,NODEMAP>      CONNECT1D;
  typedef localSearch2dA<GEOSHAPE2D,ELMAP,NODEMAP> SEARCH2D;
  typedef searchData<ELMAP>                        SEARCHDATA;
  
  //Mesh Init
  string meshFile  = "./geometries/rectangles2d/rectangleF.msh";
  
  meshInit2d<GEOSHAPE2D,ELMAP,NODEMAP> init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH1D>           grid1d = init.getGrid1d();
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT1D> connectGrid1d = init.getConnectGrid1d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  
  //Search 2d
  time_t start, end;
  Real totDist  = 0.0;
  UInt steps = 1000;
  searchData<ELMAP> outData;
  point3d P, Pmin(0.0, 0.0, 1.0), Pmax(1.0, 1.0, 0.0);
  
  
  SEARCH2D search;
  search.setMesh(grid2d,connectGrid2d);
  
  cout << "Local search init" << " "; time(&start); cout << endl;
  search.localInit(3.0);
  time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;
  
  
  cout << "Search" << " "; time(&start); cout << endl;
  geoMapInterface<GEOSHAPE2D> globMapper;
  sVect<point3d> points;
  point3d Ps;
  
  for(UInt i=0; i <= steps; ++i)
  {
    P = Pmin * (Real(i) / steps) + Pmax * (1.0 - Real(i) / steps);
    outData = search.findLocal(P);
    
    points = grid2d->getElementNodesL(outData.getElMap().getLid());
    Ps     = globMapper.getPosition(points,outData.getLocCoord());
    
    P.setZ(0.0);
    
    totDist += point3d::norm2(P-Ps);
  }
  
  time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;
  
  cout << "Last point     : " << Ps;
  cout << "Total residual : " << totDist << endl;
  cout << "Mean elements  : " << search.numFeasibElements / search.localQuery << endl;
}

