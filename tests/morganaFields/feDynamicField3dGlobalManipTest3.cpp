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

#include "feSpectralLH3d.hpp"
#include "feDynamicField3d.hpp"
#include "feDynamicField3dGlobalManip.hpp"

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
  
  
  //Mesh Init------------------------------------------------------------------
  typedef pMapItemShare             PMAPTYPE;
  typedef Real                      DOFTYPE;
  typedef feSpectralLH3d<PMAPTYPE>  FETYPE;
  typedef typename FETYPE::GEOSHAPE GEOSHAPE;
  
  string meshFile  = "../tests/morganaMeshes/mignonHexaB.unv";
  string colorFile = "../tests/morganaMeshes/mignonHexaB_color.unv";
  
  meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.femap_to_stdA(meshFile, colorFile, false);
  
  
  //Download data--------------------------------------------------------------
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH3D     MESH3D;
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT3D  CONNECT3D;
  
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  
  //Initialize the mapper
  dofMapDynamic3d_options mapOptions;
  feDynamicField3d<FETYPE,DOFTYPE,dmd3d_vectMajor,dmd3d_standard> field;
  
  field.setCommunicator(world);
  field.setGeometry(grid3d,connectGrid3d);
  field.setOptions(mapOptions);
  
  
  //Load the feCards
  typedef typename FETYPE::FECARD FECARD;
  FECARD FeCard;
  
  if(world->rank() == 0)
  {
    //el 1
    FeCard.setIsActive(true);
    FeCard.setR(1,1,2);
    field.setFeCardL(1,FeCard);
    
    //el 2
    FeCard.setIsActive(true);
    FeCard.setR(2,1,1);
    field.setFeCardL(2,FeCard);
    
    //el 3
    FeCard.setIsActive(true);
    FeCard.setR(1,2,1);
    field.setFeCardL(3,FeCard);
  }
  else
  {
    //el 1
    FeCard.setIsActive(true);
    FeCard.setR(1,2,1);
    field.setFeCardL(1,FeCard);
    
    //el 2
    FeCard.setIsActive(true);
    FeCard.setR(1,1,2);
    field.setFeCardL(2,FeCard);
    
    //el 3
    FeCard.setIsActive(true);
    FeCard.setR(2,1,1);
    field.setFeCardL(3,FeCard);
  }
  
  //Startup
  field.startup();
  
  
  //Dofs loading---------------------------------------------------------------
  if(world->rank() == 0)
  {
    for(UInt i=0; i<=11; ++i)
    {
      field.setDofL(3 * i + 1, 1.0);
      field.setDofL(3 * i + 2, 2.0);
      field.setDofL(3 * i + 3, 3.0);
    }
  }
  else
  {
    for(UInt i=0; i<=11; ++i)
    {
      field.setDofL(3 * i + 1, 3.0);
      field.setDofL(3 * i + 2, 4.0);
      field.setDofL(3 * i + 3, 2.0);
    }
  }
  
  
  //Polluting
  for(UInt i=1; i <= field.getDofVect().size(); ++i)
  {
    if( (field.getDofVect().getMapL(i).getShared()) && (!field.getDofVect().getMapL(i).getOwned()) )
    {
      field.setDofL(i, 10.0);
    }
  }
  
  
  //Testing
  fdf3dGlobalManip<PMAPTYPE,FETYPE,DOFTYPE,dmd3d_vectMajor,dmd3d_standard> tester;
  tester.setCommunicator(world);
  
  //Create the communication pattern
  typedef pMap<pMapItemSendRecv> SENDRECV;
  SENDRECV mapSend, mapRecv;
  tester.createSendRecvMap(field,mapSend,mapRecv);
  
  //Update data
  tester.updateData(field,mapSend);
  
  
  //Evaluation
  if(world->rank() == 1)
  {
    //cout << grid3d->getElements() << endl;
    //cout << field.getDofVect() << endl;
    
    Real val;
    point3d Y(0.5,0.5,0.5);
    UInt elG = 3;
    
    field.evalG(elG,Y,val);
    
    cout << "Val : " << val << endl;
  }
  
}
