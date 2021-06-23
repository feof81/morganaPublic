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
  typedef loadMesh<linearHexa,ELMAP,NODEMAP>::NODESVECT  NODESVECT;
  typedef loadMesh<linearHexa,ELMAP,NODEMAP>::GRAPH3D    GRAPH3D;
  typedef loadMesh<linearHexa,ELMAP,NODEMAP>::NODESCOLOR NODESCOLOR;
  
  string meshfile = "./geometries/hexa3d/hexa3dA.neu";
  
  UInt pid      = 1;
  UInt printPid = 0;
  UInt numPids  = 2;
  
  NODESVECT  nodes;
  GRAPH3D    elements3d;
  NODESCOLOR nodesColor;
  
  loadMesh<linearHexa,ELMAP,NODEMAP> loader;
  loader.setParam(pid,printPid,numPids);
  loader.neutral(meshfile,nodes,elements3d,nodesColor);
  
  cout << "Nodes-------------------------------------------" << endl;
  cout << nodes << endl << endl;
  
  cout << "Nodes color-------------------------------------" << endl;
  cout << nodesColor << endl << endl;
  
  cout << "Elements color----------------------------------" << endl;
  cout << elements3d << endl << endl;
}
