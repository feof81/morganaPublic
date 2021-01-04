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
  typedef GRIDSEARCH::OUTDATA OUTDATA;
  
  //Grids------------------------------------------------------------
  string meshFile = "./geometries/cubes3d/testCubeD.msh";
  
  MESHINIT init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH3D>       srcGrid3d = init.getGrid3d();
  RCP<MESH2D>       srcGrid2d = init.getGrid2d();
  RCP<CONNECT3D> srcConnect3d = init.getConnectGrid3d();
  RCP<CONNECT2D> srcConnect2d = init.getConnectGrid2d();
  
  //Mesh rebalance---------------------------------------------------
  mesh2dGlobalManip<GEOSHAPE2D,ELMAP,NODEMAP> manipMesh2d(world);
  
  RCP<MESH2D> tgtGrid2d(new MESH2D(*srcGrid2d));
  manipMesh2d.destroyOverlap(tgtGrid2d);
  manipMesh2d.meshBalancing(tgtGrid2d);
  
  RCP<CONNECT2D> tgtConnect2d(new CONNECT2D(world));
  tgtConnect2d->setMesh2d(tgtGrid2d);
  tgtConnect2d->buildConnectivity();
  
  //Gen input vector-------------------------------------------------  
  sVect<point3d> inLocCoords;
  inLocCoords.push_back(point3d(1.0/3.0, 1.0/3.0, 0.0));
  
  //Test-------------------------------------------------------------
  GRIDSEARCH gridMatcher(world);
  gridMatcher.setMesh(tgtGrid2d,
                      tgtConnect2d,
                      srcGrid3d,
                      srcConnect3d,
                      srcGrid2d,
                      srcConnect2d);
 
  OUTVECT outVect = gridMatcher.search(inLocCoords);
  
  //Output-----------------------------------------------------------
  assert(inLocCoords.size() == 1);
  
  UInt el, lid, gid;
  UInt pid = world.rank();
  point3d Pt, Ps, Y;
  Real dist = 0.0;
  
  pMapItemSendRecv      sendItem;
  pMap<pMapItemSendRecv> commMap;
  
  OUTVECT vectOther;
  pVect<point3d,ELMAP>        pCompareVect;
  geoMapSupport2d<GEOSHAPE2D> geoSupport2d;
  geoMapSupport3d<GEOSHAPE3D> geoSupport3d;
  
  for(UInt i=1; i <= outVect.size(); ++i)
  {
    sendItem.setPid(pid);
    sendItem.setLid(outVect.getMapL(i).getLid());
    sendItem.setGid(outVect.getMapL(i).getGid());
    sendItem.setSid(pid);
    
    if(outVect(i).getIsNested())
    { sendItem.setRid(outVect(i).getElMap().getPid()); }
    else
    { sendItem.setRid(outVect(i).getPData(1).getElMap().getPid()); }
    
    commMap.push_back(sendItem);
  }
  
  pVectComm<OUTDATA,ELMAP> otherComm(gridMatcher.commDev);
  otherComm.sendRecv(commMap,outVect,vectOther);  
  
  pCompareVect.resize(vectOther.size());
  
  for(UInt i=1; i <= vectOther.size(); ++i)
  {
    if(vectOther(i).getIsNested())
    {
      el = vectOther(i).getElMap().getLid();
      Y  = vectOther(i).getLocCoord(1);
    }
    else
    {
      el = vectOther(i).getPData(1).getElMap().getLid();
      Y  = vectOther(i).getPData(1).getLocCoord();
    }
    
    Pt = geoSupport3d.getPosition(srcGrid3d->getElementNodesL(el),Y);
   
    pCompareVect.getMapL(i) = vectOther.getMapL(i);
    pCompareVect(i) = Pt;
  }
  
  pVectComm<point3d,ELMAP> pCompareComm(gridMatcher.commDev);
  pCompareComm.vectorPid(pCompareVect);
  
  assert(tgtGrid2d->getNumElements() == outVect.size());
  
  for(UInt i=1; i <= tgtGrid2d->getNumElements(); ++i)
  {
    Pt = geoSupport2d.getPosition(tgtGrid2d->getElementNodesL(i),inLocCoords(1));
    
    if(outVect(i).getElMap().getPid() == pid)
    {
      if(outVect(i).getIsNested())
      {
        lid = outVect(i).getElMap().getLid();
        Y   = outVect(i).getLocCoord(1);
      }
      else
      {
        lid = vectOther(i).getPData(1).getElMap().getLid();
        Y   = vectOther(i).getPData(1).getLocCoord();
      }
      
      Ps    = geoSupport3d.getPosition(srcGrid3d->getElementNodesL(lid),Y);
      dist += point3d::norm2(Pt - Ps);
      
    }
    else
    {
      gid = outVect.getMapL(i).getGid();
      assert(pCompareVect.isG(gid));
      Ps = pCompareVect.getDataG(gid);

      dist += point3d::norm2(Pt - Ps);
    }
  }
   
  
  for(UInt i=0; i < world.size(); ++i)
  {
    if(world.rank() == i)
    { cout << "dist : " << dist << endl; }
  
    world.barrier();
  }
 
}