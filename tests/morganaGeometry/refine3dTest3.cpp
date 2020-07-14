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
#include "geoMapSupport3d.hpp"
#include "meshInit3d.hpp"

#include "connect2d.hpp"
#include "connect3d.hpp"
#include "printMesh.hpp"
#include "printMesh2dHDF5.hpp"
#include "printMesh3dHDF5.hpp"

#include "refine3d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
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
  
  //! Mesh Init------------------------------------------------------
  string meshFile = "./geometries/cubes3d/testCubeC.msh";
  
  meshInit3d<GEOSHAPE3D,ELMAP,NODEMAP> init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //! Elements tbr---------------------------------------------------
  sVect<point3d> nodes;
  sVect<UInt> elementsTbr;
  point3d P;
  Real toll = 0.05;
  
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    nodes = grid3d->getElementNodesL(i);
    P = (nodes(1) + nodes(2) + nodes(3) + nodes(4)) / 4.0;
    
    if((P.getX() <= toll) && (P.getY() <= toll) && (P.getZ() <= toll))
    { elementsTbr.push_back(i); }
  }
  
  cout << "pid : " << world.rank() << " numElements : " << elementsTbr.size() << endl;
  
  
  //! Refinement-----------------------------------------------------
  typedef refine3d<GEOSHAPE3D,ELMAP,NODEMAP>::POINT_EDGE POINT_EDGE;
  
  refine3d<GEOSHAPE3D,ELMAP,NODEMAP> raffinator(world);
  raffinator.upload(grid3d,grid2d);
  
  Real tollH      = 0.1;
  UInt rFactor    = 1;
  UInt edgesBlock = 5;
  UInt maxSteps   = 600;
  
  raffinator.setRefinementParams(tollH,rFactor,elementsTbr);
  raffinator.refineLeb(edgesBlock, maxSteps);
  raffinator.download(grid3d,grid2d);
  
  //! Printing-------------------------------------------------------  
  printMesh3dHDF5<GEOSHAPE3D,ELMAP> printer3d(world);
  printer3d.printGeoIds("griglia3d",0,grid3d);
  
  printMesh2dHDF5<GEOSHAPE2D,ELMAP> printer2d(world);
  printer2d.printGeoIds("griglia2d",0,grid2d);
  
  //! Compute volume-------------------------------------------------
  geoMapSupport3d<GEOSHAPE3D> geoSupport;
  Real locVol, vol = 0.0;
  
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    if(grid3d->getElements().getRowMapL(i).getOwned())
    {
      geoSupport.setPoints(grid3d->getElementNodesL(i));
      locVol = geoSupport.volume<STANDARD,2> ();
      vol   += locVol;
      
      assert(locVol);
    }
  }
  
  std::plus<Real> op;
  all_reduce(world, vol, vol, op);
  cout << "vol: " << vol << endl;
}
