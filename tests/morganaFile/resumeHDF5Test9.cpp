/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "pMapItem.h"
#include "pMapItemShare.h"
#include "geoShapes.h"
#include "meshInit3d.hpp"

#include "resumeHDF5.h"


using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef linearTetra    GEOSHAPE;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  typedef mesh3d<GEOSHAPE,ELMAP,NODEMAP>  MESH3D;
  
  string meshFile = "./tests/morganaMeshes/mignon3dB.msh";
  //string meshFile = "./geometries/cubes3d/testCubeA.msh";
  
  meshInit3d<GEOSHAPE,ELMAP,NODEMAP> init(world);
  init.gmMesh_to_stdA(meshFile);
  
  
  RCP<MESH3D> grid3d = init.getGrid3d(); 
  MESH3D rGrid;
  
  
  //Loading and downloading
  resumeHDF5 resumer(world);
  resumer.printToFile("mesh3d",*grid3d);
  resumer.loadFromFile("mesh3d", rGrid);
  
  if(world.rank() == 1)
  { cout << rGrid.getElements() << endl; }
}
