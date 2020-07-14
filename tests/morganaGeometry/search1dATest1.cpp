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
#include "search1dA.hpp"

#include "printMesh.hpp"


using namespace std;
using namespace boost::mpi;
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef linearTriangle GEOSHAPE2D;
  typedef linearLine     GEOSHAPE1D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>     MESH2D;
  typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>     MESH1D;
  typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>  CONNECT2D;
  typedef connect1d<GEOSHAPE1D,ELMAP,NODEMAP>  CONNECT1D;
  typedef search1dA<GEOSHAPE1D,ELMAP,NODEMAP>  SEARCH1D;
  typedef searchData<ELMAP>  SEARCHDATA;
  
  //Mesh Init
  UInt pid = world.rank();
  time_t start, end;
  string meshFile  = "./geometries/rectangles2d/rectangleE.msh";
  
  meshInit2d<GEOSHAPE2D,ELMAP,NODEMAP> init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH1D>           grid1d = init.getGrid1d();
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT1D> connectGrid1d = init.getConnectGrid1d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  //Search 2d Init
  SEARCH1D search(world);
  search.setMesh(grid1d,connectGrid1d);
  
  if(pid == 0) {cout << "Local search init" << " "; time(&start); cout << endl;}
  search.localInit(3.0);
  if(pid == 0) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  if(pid == 0) {cout << "Global search init" << " "; time(&start); cout << endl;}
  search.globalInit();
  if(pid == 0) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Search 2d
  pVect<point3d,NODEMAP>    Pvect;
  pVect<SEARCHDATA,NODEMAP> dataVect;
  
  UInt gid, steps = 10000;
  point3d P, Pmin(0.0, 0.0, 0.0), Pmax(1.0,0.0,0.0);
  
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
  Real globDist  = 0.0;
  geoMapInterface<GEOSHAPE1D> interfaceMap;
  sVect<point3d> elNodes;
  point3d Y;
  
  assert(Pvect.size() == dataVect.size());
  
  for(UInt i=1; i <= dataVect.size(); ++i)
  {
    if(pid == dataVect(i).getElMap().getPid())
    {
      gid     = dataVect(i).getElMap().getGid();
      elNodes = grid1d->getElementNodesG(gid);
      Y       = dataVect(i).getLocCoord();
      
      P        = interfaceMap.getPosition(elNodes,Y);
      locDist += point3d::norm2(Pvect(i) - P);
      
      assert( !(dataVect.getMapL(i) != Pvect.getMapL(i)) );
      
      totCheck++;
    }
    
    globDist += dataVect(i).getDistance();
  }
  
  cout << endl;
  cout << "Pid        : " << pid << endl;
  cout << "Local Res  : " << locDist << endl;
  cout << "Global Res : " << globDist << endl;
  cout << "TotCheck   : " << totCheck << endl;
  cout << "NumEls     : " << grid1d->getNumElements() << endl;
}

