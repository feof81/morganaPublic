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
  typedef searchData<ELMAP>  SEARCHDATA;
  
  //Mesh Init
  UInt pid = world.rank();
  time_t start, end;
  string meshFile  = "./geometries/hexa3d/hexa3dA.neu";
  string colorFile = "./geometries/hexa3d/hexa3dA_color.neu";
  
  meshInit3d<GEOSHAPE3D,ELMAP,NODEMAP> init(world);
  init.neutral_to_stdA(meshFile, colorFile);
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //Search 3d Init
  SEARCH3D search(world);
  search.setMesh(grid2d,grid3d,connectGrid2d,connectGrid3d);
  
  if(pid == 0) {cout << "Local search init" << " "; time(&start); cout << endl;}
  search.localInit();
  if(pid == 0) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  if(pid == 0) {cout << "Global search init" << " "; time(&start); cout << endl;}
  search.globalInit();
  if(pid == 0) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Search 3d
  pVect<point3d,NODEMAP>    Pvect;
  pVect<SEARCHDATA,NODEMAP> dataVect;
  
  UInt gid, steps = 5000;
  point3d P, Pmin(0.0, -1.0, 0.0), Pmax(1.0, -1.0, 0.0);
  
  for(UInt i=0; i <= steps; ++i)
  {
    gid = 1 + i + pid * (steps + 1);
    P   = Pmin * (Real(i) / steps) + Pmax * (1.0 - Real(i) / steps);
    Pvect.push_back(P,pMapItemShare(i+1, gid, pid, false, false));
  }
  Pvect.updateFinder();
  
  
  //Finding
  world.barrier();
  
  if(pid == 0) {cout << "Searching time" << " "; time(&start); cout << endl;}
  search.findGlobal(Pvect,dataVect);
  if(pid == 0) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  search.printPerformances(0);
  
  
  //Checking
  UInt totCheck  = 0;
  Real locDist   = 0.0;
  geoMapInterface<GEOSHAPE3D> interfaceMap;
  sVect<point3d> elNodes;
  point3d Y;
  
  assert(Pvect.size() == dataVect.size());
  
  for(UInt i=1; i <= dataVect.size(); ++i)
  {
    if(pid == dataVect(i).getElMap().getPid())
    {
      gid     = dataVect(i).getElMap().getGid();
      elNodes = grid3d->getElementNodesG(gid);
      Y       = dataVect(i).getLocCoord();
      
      P = interfaceMap.getPosition(elNodes,Y);
      P.setY(-1.0);
      
      locDist += point3d::norm2(Pvect(i) - P);
      
      assert( !(dataVect.getMapL(i) != Pvect.getMapL(i)) );
      
      totCheck++;
    }
  }
  
  cout << endl;
  cout << "Pid        : " << pid << endl;
  cout << "Local Res  : " << locDist << endl;
  cout << "TotCheck   : " << totCheck << endl;
}
