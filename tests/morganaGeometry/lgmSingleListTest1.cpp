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
#include "lgmSingleList.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef linearTetra    GEOSHAPE3D;
  typedef linearTriangle GEOSHAPE2D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  
  typedef mesh3d<GEOSHAPE3D,ELMAP,NODEMAP>     MESH3D;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>     MESH2D;
  typedef connect3d<GEOSHAPE3D,ELMAP,NODEMAP>  CONNECT3D;
  typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>  CONNECT2D;
  typedef meshInit3d<GEOSHAPE3D,ELMAP,NODEMAP> MESHINIT;
  
  typedef lgmSingleList<MESH3D,CONNECT3D,3,3> LGMLIST33;
  typedef lgmSingleList<MESH3D,CONNECT3D,3,2> LGMLIST32;
  typedef lgmSingleList<MESH3D,CONNECT3D,3,1> LGMLIST31;
  
  typedef LGMLIST33::STDMAP STDMAP33;
  typedef LGMLIST32::STDMAP STDMAP32;
  typedef LGMLIST31::STDMAP STDMAP31;
  
  typedef STDMAP33::iterator ITER33;
  typedef STDMAP32::iterator ITER32;
  typedef STDMAP31::iterator ITER31;
  
  string meshFile = "./tests/morganaMeshes/mignon3dA.msh";
  
  MESHINIT init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH3D> grid3d = init.getGrid3d();
  RCP<MESH2D> grid2d = init.getGrid2d();
  
  RCP<CONNECT3D> connect3d = init.getConnectGrid3d();
  RCP<CONNECT2D> connect2d = init.getConnectGrid2d();
  
  
  LGMLIST33 listFeeder33(grid3d,connect3d);
  LGMLIST32 listFeeder32(grid3d,connect3d);
  LGMLIST31 listFeeder31(grid3d,connect3d);
  
  STDMAP33 stdMap33 = listFeeder33.getList();
  STDMAP32 stdMap32 = listFeeder32.getList();
  STDMAP31 stdMap31 = listFeeder31.getList();
  
  cout << "Grid------------------------------------------" << endl;
  cout << *grid3d << endl;
    
  cout << "StdMap-3-3------------------------------------" << endl;  
  for(ITER33 iter = stdMap33.begin(); iter != stdMap33.end(); ++iter)
  { cout << iter->second << " " << iter->first << endl; }
  
  cout << "StdMap-3-2------------------------------------" << endl;  
  for(ITER32 iter = stdMap32.begin(); iter != stdMap32.end(); ++iter)
  { cout << iter->second << " " << iter->first << endl; }
  
  cout << "StdMap-3-1------------------------------------" << endl;  
  for(ITER31 iter = stdMap31.begin(); iter != stdMap31.end(); ++iter)
  { cout << iter->second << " " << iter->first << endl; }
}