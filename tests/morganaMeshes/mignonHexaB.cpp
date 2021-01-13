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

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef linearHexa     GEOSHAPE;
  typedef linearQuad     GEOSHAPE2D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  typedef mesh3d<GEOSHAPE,ELMAP,NODEMAP>    MESH3D;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>  MESH2D;
  
  
  string meshFile  = "./tests/morganaMeshes/mignonHexaB.unv";
  string colorFile = "./tests/morganaMeshes/mignonHexaB_color.unv";
  
  meshInit3d<GEOSHAPE,ELMAP,NODEMAP> init(world);
  init.femap_to_stdA(meshFile, colorFile);
  
  //Download
  RCP<MESH2D> grid2d = init.getGrid2d();
  RCP<MESH3D> grid3d = init.getGrid3d();
  
  //Printout
  if(world.rank() == 1)
  {
    cout << "Nodes 3d" << endl;
    cout << grid3d->getNodes() << endl << endl;
    
    //cout << "Elements 3d" <<  endl;
    //cout << grid3d->getElements() << endl << endl;
    
    /*cout << "Faces 3d" <<  endl;
    cout << grid3d->getFaces() << endl << endl;
    
    cout << "Edges 3d" <<  endl;
    cout << grid3d->getEdges() << endl << endl;
    
    
    cout << "Nodes 2d" << endl;
    cout << grid2d->getNodes() << endl << endl;
    
    cout << "Elements 2d" <<  endl;
    cout << grid2d->getElements() << endl << endl;
    
    cout << "Edges 2d" <<  endl;
    cout << grid2d->getEdges() << endl << endl;*/
  }
}
