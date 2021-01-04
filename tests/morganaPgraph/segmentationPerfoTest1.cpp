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

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 1);
  
  typedef linearTetra    GEOSHAPE;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  
  //string meshFile = "./tests/morganaMeshes/mignon3dC.msh";
  string meshFile = "./geometries/cubes3d/testCubeF.msh";
  
  meshInit3d<GEOSHAPE,ELMAP,NODEMAP> init(world);
  init.gmMesh_to_stdA(meshFile);
  
  typedef typename meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::MESH2D  MESH2D;
  typedef typename meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::MESH3D  MESH3D;
  typedef typename MESH3D::GRAPH2D                             GRAPH2D;
  typedef typename MESH3D::GEOELEMENT2D                        GEOELEMENT2D;
  
  RCP<MESH2D> grid2d = init.getGrid2d();
  RCP<MESH3D> grid3d = init.getGrid3d();

  GRAPH2D faces = grid3d->getFaces();
  
  //Testing
  GEOELEMENT2D dataMax, dataMin;
  
  pVectGlobalManip<GEOELEMENT2D,ELMAP> globManip(world);
  globManip.dataMinMax(faces,dataMin,dataMax);
  
  cout << "NumNodes: " << grid3d->getNumNodes() << endl << endl;
  cout << "NumFaces: " << faces.size() << endl;
  cout << "Data Min: " << dataMin << endl;
  cout << "Data Max: " << dataMax << endl << endl;
  
  //Segmentation
  UInt numSegments = 10;
  
  typedef pVect<GEOELEMENT2D,ELMAP> PVECT;
  sVect<PVECT> targetVects;
  
  pVectManip<GEOELEMENT2D,ELMAP> segmenter;
  segmenter.segmentationData(targetVects,numSegments,dataMin,dataMax,faces);
  
  for(UInt i=1; i <= numSegments; ++i)
  {
    cout << "pid " << i << ": " << (Real(targetVects(i).size()) / Real(faces.size())) * 100.0  << "%" << endl;
  }
}