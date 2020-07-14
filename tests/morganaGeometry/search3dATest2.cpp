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
#include "search3dA.hpp"

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
  
  typedef linearHexa     GEOSHAPE3D;
  typedef linearQuad     GEOSHAPE2D;
  typedef linearLine     GEOSHAPE1D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  typedef mesh3d<GEOSHAPE3D,ELMAP,NODEMAP>     MESH3D;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>     MESH2D;
  typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>     MESH1D;
  typedef connect3d<GEOSHAPE3D,ELMAP,NODEMAP>  CONNECT3D;
  typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>  CONNECT2D;
  typedef search3dA<GEOSHAPE3D,ELMAP,NODEMAP>  SEARCH3D;
  
  
  //Mesh Init
  UInt pid = world.rank();
  time_t start, end;
  
  string meshFile  = "./tests/morganaMeshes/mignonHexaB.unv";
  string colorFile = "./tests/morganaMeshes/mignonHexaB_color.unv";
  
  meshInit3d<GEOSHAPE3D,ELMAP,NODEMAP> init(world);
  init.femap_to_stdA(meshFile, colorFile);
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  
  //Search 3d
  SEARCH3D search(world);
  search.setMesh(grid2d,grid3d,connectGrid2d,connectGrid3d);
  
  search.localInit(3.0);
  search.globalInit();
  
  if(world.rank() == 0)
  {
    sVect<UInt> indices = search.globalMatchingPids(point3d(0.5,0.5,0.0));
  
    for(UInt i=1; i <= indices.size(); ++i)
    {
      cout << indices(i) << " ";
    }
    cout << endl << endl;
  
  
    point3d P(2.0, 2.0, 0.5);
    point3d Y;
    
    typedef SEARCH3D::SEARCHDATA SEARCHDATA;
    SEARCHDATA data = search.findLocal(P);
    cout << "Search: " << data << endl;
  
    geoMapInterface<GEOSHAPE3D> support3d;
    sVect<point3d> points = grid3d->getElementNodesL(data.getElMap().getLid());
    point3d Ps = support3d.getPosition(points,Y);
  
    cout << "Y : " << Y;
    cout << "P : " << P;
    cout << "Ps: " << Ps; 
  }
  

  //search.printCells(1);
  //search.printGlobal(1);
}

/*
Search: 
Element          : pid: 0 lid: 2 gid: 2 shared: 1 owned 1
Local coordinate : 1 1 0.5
Distance         : 1
Found            : 1

Y : 0 0 0
P : 2 2 0.5
Ps: 1 0 0
 */
