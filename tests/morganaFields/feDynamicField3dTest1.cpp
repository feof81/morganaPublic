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
#include "search3dA.hpp"

#include "feXFEM_LT3d.hpp"
#include "feDynamicField3d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  Teuchos::RCP<communicator> world(new communicator);
  
  assert(world->size() == 2);
  
  UInt pid = world->rank();
  
  typedef linearTetra    GEOSHAPE;
  typedef pMapItemShare  PMAPTYPE;
  
  string meshFile = "../tests/morganaMeshes/mignon3dC.msh";
  //string meshFile = "./geometries/cubes3d/testCubeD.msh";
  
  meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Download data
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH2D     MESH2D;
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH3D     MESH3D;
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT2D  CONNECT2D;
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT3D  CONNECT3D;
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  
  //Initialize the mapper
  typedef Real                    DOFTYPE;
  typedef feXFEM_LT3d<1,PMAPTYPE> FETYPE;
  
  dofMapDynamic3d_options mapOptions;
  feDynamicField3d<FETYPE,DOFTYPE,dmd3d_vectMajor,dmd3d_standard> field;
  
  field.setCommunicator(world);
  field.setGeometry(grid3d,connectGrid3d);
  field.setOptions(mapOptions);
  
  
  //Load the feCards
  typedef typename FETYPE::FECARD FECARD;
  FECARD FeCard;
  
  FeCard.setActive(true);
  FeCard.setXtoll(0.1);  
  
  if(world->rank() == 0)
  {
    //el 1
    FeCard.getPhi(1)  = 0.0;
    FeCard.getPhi(2)  = 0.0;
    FeCard.getPhi(3)  = 0.0;
    FeCard.getPhi(4)  = 0.0;
    
    FeCard.resizeEnrichedDof(3);
    FeCard.getEnrichedDof(1)  = 1;
    FeCard.getEnrichedDof(2)  = 2;
    FeCard.getEnrichedDof(3)  = 3;
    field.setFeCardL(1,FeCard);
    
    
    //el 2
    FeCard.getPhi(1)  = 0.0;
    FeCard.getPhi(2)  = 0.0;
    FeCard.getPhi(3)  = 0.0;
    FeCard.getPhi(4)  = 0.0;
    
    FeCard.resizeEnrichedDof(3);
    FeCard.getEnrichedDof(1)  = 1;
    FeCard.getEnrichedDof(2)  = 2;
    FeCard.getEnrichedDof(3)  = 3;
    field.setFeCardL(2,FeCard);
    
    
    //el3    
    FeCard.getPhi(1)  = 0.0;
    FeCard.getPhi(2)  = 0.0;
    FeCard.getPhi(3)  = 0.0;
    FeCard.getPhi(4)  = 1.0;
    
    FeCard.resizeEnrichedDof(4);
    FeCard.getEnrichedDof(1)  = 1;
    FeCard.getEnrichedDof(2)  = 2;
    FeCard.getEnrichedDof(3)  = 3;
    FeCard.getEnrichedDof(4)  = 4;
    field.setFeCardL(3,FeCard);
    
    
    //el4
    FeCard.getPhi(1)  = 0.0;
    FeCard.getPhi(2)  = 0.0;
    FeCard.getPhi(3)  = 0.0;
    FeCard.getPhi(4)  = 1.0;
    
    FeCard.resizeEnrichedDof(4);
    FeCard.getEnrichedDof(1)  = 1;
    FeCard.getEnrichedDof(2)  = 2;
    FeCard.getEnrichedDof(3)  = 3;
    FeCard.getEnrichedDof(4)  = 4;
    field.setFeCardL(4,FeCard);
  }
  else
  {
    //el1
    FeCard.getPhi(1)  = 0.0;
    FeCard.getPhi(2)  = 0.0;
    FeCard.getPhi(3)  = 0.0;
    FeCard.getPhi(4)  = 1.0;
    
    FeCard.resizeEnrichedDof(4);
    FeCard.getEnrichedDof(1)  = 1;
    FeCard.getEnrichedDof(2)  = 2;
    FeCard.getEnrichedDof(3)  = 3;
    FeCard.getEnrichedDof(4)  = 4;
    field.setFeCardL(1,FeCard);
    
    //el2
    FeCard.getPhi(1)  = 0.0;
    FeCard.getPhi(2)  = 0.0;
    FeCard.getPhi(3)  = 0.0;
    FeCard.getPhi(4)  = 1.0;
    
    FeCard.resizeEnrichedDof(4);
    FeCard.getEnrichedDof(1)  = 1;
    FeCard.getEnrichedDof(2)  = 2;
    FeCard.getEnrichedDof(3)  = 3;
    FeCard.getEnrichedDof(4)  = 4;
    field.setFeCardL(2,FeCard);
    
    //el3
    FeCard.getPhi(1)  = 0.0;
    FeCard.getPhi(2)  = 0.0;
    FeCard.getPhi(3)  = 0.0;
    FeCard.getPhi(4)  = 0.0;
    
    FeCard.resizeEnrichedDof(3);
    FeCard.getEnrichedDof(1)  = 1;
    FeCard.getEnrichedDof(2)  = 2;
    FeCard.getEnrichedDof(3)  = 3;
    field.setFeCardL(3,FeCard);
    
    //el4
    FeCard.getPhi(1)  = 0.0;
    FeCard.getPhi(2)  = 0.0;
    FeCard.getPhi(3)  = 0.0;
    FeCard.getPhi(4)  = 0.0;
    
    FeCard.resizeEnrichedDof(3);
    FeCard.getEnrichedDof(1)  = 1;
    FeCard.getEnrichedDof(2)  = 2;
    FeCard.getEnrichedDof(3)  = 3;
    field.setFeCardL(4,FeCard);
  }
  
  
  //Startup
  field.startup();
  
  //Dofs loading
  if(world->rank() == 0)
  {
    field.setDofL(1, 1.0);
    field.setDofL(2, 1.0);
    field.setDofL(3, 1.0);
    field.setDofL(4, 1.0);
    field.setDofL(5, 1.0);
    field.setDofL(6, 1.0);
    
    field.setDofL(7,   0.0);
    field.setDofL(8,   0.0);
    field.setDofL(9,   0.0);
    field.setDofL(10, -2.0);
    field.setDofL(11,  0.0);
  }
  else
  {
    field.setDofL(1, 1.0);
    field.setDofL(2, 1.0);
    field.setDofL(3, 1.0);
    field.setDofL(4, 1.0);
    field.setDofL(5, 1.0);
    field.setDofL(6, 1.0);
    
    field.setDofL(7,   0.0);
    field.setDofL(8,   0.0);
    field.setDofL(9,   0.0);
    field.setDofL(10, -2.0);
    field.setDofL(11,  0.0);
  }
  
  
  //Search 3d Init
  typedef search3dA<GEOSHAPE,PMAPTYPE,PMAPTYPE>  SEARCH3D;
  typedef searchData<PMAPTYPE>                   SEARCHDATA;
  
  SEARCH3D search(world);
  search.setMesh(grid2d,
                 grid3d,
                 connectGrid2d,
                 connectGrid3d);
  
  search.localInit(3.0);
  search.globalInit();
  
  
  //Global search point
  pVect<point3d,PMAPTYPE>    Pvect;
  pVect<SEARCHDATA,PMAPTYPE> dataVect;
  
  point3d P,Y;
  Real V, Vx, Vy, Vz;
  
  P.setX(0.5 + 1e-6); P.setY(0.5); P.setZ(0.5);
  Pvect.push_back(P,PMAPTYPE(1,1,pid,false,false));
  Pvect.updateFinder();
  
  search.findGlobal(Pvect,dataVect);
  
  UInt elL = dataVect(1).getElMap().getLid();
  UInt elG = dataVect(1).getElMap().getGid();
         Y = dataVect(1).getLocCoord();
  
  
  //Evaluate
  if(pid == dataVect(1).getElMap().getPid())
  {
    field.evalL(elL,Y,V);
    field.gradG(elG,Y,Vx,Vy,Vz);
    cout << "pid " << pid << endl << " eval " << V << " gradX " << Vx << " gradY " << Vy << " gradZ " << Vz << endl;
  }
}
