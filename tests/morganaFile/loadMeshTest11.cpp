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

#include "point3d.h"
#include "simpleFormats.hpp"
#include "pVect.hpp"
#include "pMapItem.h"

#include "loadMesh.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;

int main(int argc, char *argv[])
{
  typedef pMapItem ELMAP;
  typedef pMapItem NODEMAP;
  typedef loadMesh<linearTetra,ELMAP,NODEMAP>::NODESVECT  NODESVECT;
  typedef loadMesh<linearTetra,ELMAP,NODEMAP>::GRAPH3D    GRAPH3D;
  typedef loadMesh<linearTetra,ELMAP,NODEMAP>::GRAPH2D    GRAPH2D;
  
  environment  env(argc,argv);
  Teuchos::RCP<communicator> world(new communicator);
  
  string meshfile = "./tests/morganaMeshes/mignon3dD.msh";
  
  UInt pid      = world->rank();
  UInt printPid = 0;
  UInt numPids  = world->size();
  
  NODESVECT nodes, nodesB;
  GRAPH3D   elements3d, elements3dB;
  GRAPH2D   elements2d, elements2dB;
  
  loadMesh<linearTetra,ELMAP,NODEMAP> loader;
  loader.setParam(pid,printPid,numPids);
  loader.gmMesh(meshfile,nodes,elements3d,elements2d);
  loader.gmMeshParallel(meshfile,world,nodesB,elements3dB,elements2dB);
  
  
  world->barrier();
  if(pid == 0)
  {
    cout << "Nodes" << endl;
    cout << nodes << endl;
    
    cout << "NodesB" << endl;
    cout << nodesB << endl;
    
    cout << "Elements2d" << endl;
    cout << elements2d << endl;
    
    cout << "Elements2dB" << endl;
    cout << elements2dB << endl;
    
    cout << "Elements3d" << endl;
    cout << elements3d << endl;
    
    cout << "Elements3dB" << endl;
    cout << elements3dB << endl;
  }
  sleep(1.0);
  
  world->barrier();
  if(pid == 1)
  {
    cout << "Nodes" << endl;
    cout << nodes << endl;
    
    cout << "NodesB" << endl;
    cout << nodesB << endl;
    
    cout << "Elements2d" << endl;
    cout << elements2d << endl;
    
    cout << "Elements2dB" << endl;
    cout << elements2dB << endl;
    
    cout << "Elements3d" << endl;
    cout << elements3d << endl;
    
    cout << "Elements3dB" << endl;
    cout << elements3dB << endl;
  }
  sleep(1.0);
}

