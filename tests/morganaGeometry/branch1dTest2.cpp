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
#include "meshInit1d.hpp"
#include "printMesh1dHDF5.hpp"
#include "branch1d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  //Comm startup-----------------------------------------------------
  environment  env(argc,argv);
  communicator world;
  
  //Typedefs---------------------------------------------------------
  typedef linearLine     GEOSHAPE1D;
  typedef pMapItemShare  ELMAP;
  typedef pMapItemShare  NODEMAP;
  
  typedef meshInit1d<GEOSHAPE1D,ELMAP,NODEMAP> INIT;
  typedef INIT::MESH1D    MESH1D;
  typedef INIT::CONNECT1D CONNECT1D;
  typedef branch1d<GEOSHAPE1D,ELMAP,NODEMAP> BRANCH1D;
  
  //Loading----------------------------------------------------------
  string meshFile = "./geometries/cadInterface/mignonBranch.msh";
  
  INIT init(world);
  init.gmMesh_to_stdA(meshFile);
  
  Teuchos::RCP<MESH1D>    mesh1d        = init.getGrid1d();
  Teuchos::RCP<CONNECT1D> meshConnect1d = init.getConnectGrid1d();
  
  //Printing---------------------------------------------------------
  printMesh1dHDF5<GEOSHAPE1D,ELMAP> printer(world);
  printer.printGeoIds("outGrid",0,mesh1d);
  
  //Branch test------------------------------------------------------
  Teuchos::RCP<BRANCH1D> branch(new BRANCH1D(*mesh1d));
  branch->setCommunicator(world);
  branch->setConnect1d(meshConnect1d);
  branch->computeNumVertices();
  branch->transferMap();
  branch->buildTriplePoints();
  
  UInt numEndNodes = branch->getNumEndNodes();
  
  for(UInt k=0; k <= world.size(); ++k)
  {
    if(k == world.rank())
    {
      cout << "Proc - " << k << "-------------------------------" << endl;
      //cout << "Mesh - nodes" << endl;
      //cout << mesh1d->getNodes() << endl;
      
      //cout << "Mesh - elements" << endl;
      //cout << mesh1d->getElements() << endl;
      
      cout << "numEndNodes : " << numEndNodes << endl;
      cout << "endNodesL   : "; 
      for(UInt i=1; i <= numEndNodes; ++i) { cout << branch->getEndNodeToNodeL(i) << " "; } cout << endl;
  
      cout << "endNodesG   : ";
      for(UInt i=1; i <= numEndNodes; ++i) { cout << branch->getEndNodeToNodeG(i) << " "; } cout << endl;
  
      cout << "endNodeToIntElementL : ";
      for(UInt i=1; i <= numEndNodes; ++i) { cout << branch->getEndNodeToIntElementL(i) << " "; } cout << endl;
  
      cout << "endNodeToIntElementG : ";
      for(UInt i=1; i <= numEndNodes; ++i) { cout << branch->getEndNodeToIntElementG(i) << " "; } cout << endl;
      
      cout << "endNodeToIntNormal   : ";
      for(UInt i=1; i <= numEndNodes; ++i) { cout << branch->getEndNodeToIntNormal(i) << " "; } cout << endl;
                    
      cout << "getEndNodeToExtNodes : ";
      for(UInt i=1; i <= numEndNodes; ++i) { cout << branch->getEndNodeToExtNodeG(i) << " "; } cout << endl;
                    
      cout << "getNumEndNodeToExtElements : " << endl;
      for(UInt i=1; i <= numEndNodes; ++i)
      {
        cout << branch->getNumEndNodeToExtElements(i) << " / ";
    
        for(UInt j=1; j <= branch->getNumEndNodeToExtElements(i); ++j)
        { cout << branch->getEndNodeToExtElementsG(i,j) << " "; }
    
        cout << endl;
      }
      cout << endl;

      cout << "numNormals: " << endl;
      for(UInt i=1; i <= numEndNodes; ++i)
      {
        cout << branch->getNumEndNodeToExtNormals(i) << " / ";

        for(UInt j=1; j <= branch->getNumEndNodeToExtNormals(i); ++j)
        { cout << branch->getEndNodeToExtNormals(i,j); }
    
        cout << endl;
      }
      
      cout << endl << endl;
    }
    
    world.barrier();
  }
}
