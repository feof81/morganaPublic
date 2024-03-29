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

#include "printMesh.hpp"
#include "printMesh2dHDF5.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef linearTetra    GEOSHAPE;
  typedef linearTriangle GEOSHAPE2D;
  typedef pMapItemShare  PMAPTYPE;
  typedef mesh3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH3D;
  typedef mesh2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>  MESH2D;
  
  
  //Mesh Init
  string meshFile = "./geometries/cadInterface/Parte1.msh";
  
  meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdA(meshFile);
  
  
  //Download data
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH2D     MESH2D;
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT2D  CONNECT2D;
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  //Printer
  printMesh2dHDF5<GEOSHAPE2D,PMAPTYPE> printer(world);
  printer.printGeoIds("testGrid",0,grid2d);
}
