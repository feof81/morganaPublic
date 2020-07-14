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
#include "pMap.hpp"

#include "pVect.hpp"
#include "pVectComm.hpp"

#include "pGraphItem.h"
#include "pGraph.hpp"
#include "pGraphManip.hpp"
#include "pGraphComm.hpp"
#include "pGraphGlobalManip.hpp"

#include "mesh3d.hpp"
#include "mesh3dGlobalManip.hpp"

#include "traitsMpiOptimization.hpp"

using namespace Teuchos;


//! Run with two processors
int main(int argc, char *argv[])
{
  //Typedefs---------------------------------------------------------
  typedef linearTetra          GEOSHAPE;
  typedef geoElement<GEOSHAPE> GEOELEMENT;
  typedef pMapItemShare ELMAP;
  typedef pMapItem      NODEMAP;
  
  //Comm dev---------------------------------------------------------
  boost::mpi::environment  env(argc, argv);
  boost::mpi::communicator newCommDev;
  
  bool commActive = (newCommDev.rank() < 2);
  boost::mpi::communicator oldCommDev = newCommDev.split(commActive);
  
  assert(newCommDev.size() == 4);
  
  //Alloc------------------------------------------------------------
  mesh3d<GEOSHAPE,ELMAP,NODEMAP>            newGrid3d;
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> fixer;
  
  //Testing----------------------------------------------------------
  if(commActive)
  {
    mesh3d<GEOSHAPE,ELMAP,NODEMAP> oldGrid3d;
    
    if(oldCommDev.rank() == 0)
    {
      //Nodes
      pVect<point3d,NODEMAP> nodes;
      
      nodes.reserve(5);
      nodes.push_back(point3d(0.0, 0.0, 0.0),pMapItem(1,1));
      nodes.push_back(point3d(1.0, 0.0, 0.0),pMapItem(2,2));
      nodes.push_back(point3d(1.0, 1.0, 0.0),pMapItem(3,3));
      nodes.push_back(point3d(0.0, 1.0, 0.0),pMapItem(4,4));
      nodes.push_back(point3d(0.5, 0.5, 1.0),pMapItem(5,5));
      nodes.updateFinder();
      
      //Elements
      GEOELEMENT tet(true);
      pGraph<GEOELEMENT,ELMAP,NODEMAP> elList;
      
      elList.reserve(2);
      
      tet.setGeoId(1); tet(1) = 1; tet(2) = 2; tet(3) = 3; tet(4) = 5;
      elList.push_back(tet, pMapItemShare(1,1,true,true));
      tet.setGeoId(1); tet(1) = 1; tet(2) = 3; tet(3) = 4; tet(4) = 5;
      elList.push_back(tet, pMapItemShare(2,2,true,false));
      
      //The grid3d
      oldGrid3d.setNodes(nodes);
      oldGrid3d.setElements(elList);
    }
  
    if(oldCommDev.rank() == 1)
    {
      //Nodes
      pVect<point3d,NODEMAP> nodes;
      
      nodes.reserve(5);
      nodes.push_back(point3d(0.0, 0.0, 0.0),pMapItem(1,1));
      nodes.push_back(point3d(1.0, 0.0, 0.0),pMapItem(2,2));
      nodes.push_back(point3d(1.0, 1.0, 0.0),pMapItem(3,3));
      nodes.push_back(point3d(0.0, 1.0, 0.0),pMapItem(4,4));
      nodes.push_back(point3d(0.5, 0.5, 1.0),pMapItem(5,5));
      nodes.updateFinder();
      
      //Elements
      GEOELEMENT tet(true);
      pGraph<GEOELEMENT,ELMAP,NODEMAP> elList;
      
      elList.reserve(2);
      
      tet.setGeoId(1); tet(1) = 1; tet(2) = 2; tet(3) = 3; tet(4) = 5;
      elList.push_back(tet, pMapItemShare(1,1,true,false));
      tet.setGeoId(1); tet(1) = 1; tet(2) = 3; tet(3) = 4; tet(4) = 5;
      elList.push_back(tet, pMapItemShare(2,2,true,true));
      
      //The grid3d
      oldGrid3d.setNodes(nodes);
      oldGrid3d.setElements(elList);
    }
    
    fixer.expandCommunicator(commActive,
                             oldCommDev,
                             oldGrid3d,
                             newCommDev,
                             newGrid3d);
  }
  
  newGrid3d.computeNumVertices();
  
  //Printout---------------------------------------------------------
  for(UInt pid=0; pid < newCommDev.size(); pid++)
  {
    if(pid == newCommDev.rank())
    {
      cout << "Process : " << pid << " -------------------------- " << endl;
      cout << newGrid3d.getNodes()    << endl;
      cout << newGrid3d.getNumVertices() << endl;
      cout << newGrid3d.getElements() << endl;
      cout << newGrid3d.getFaces()    << endl;
      cout << newGrid3d.getEdges()    << endl;
    }
    
    newCommDev.barrier();
  }
}
