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
#include "meshInit2d.hpp"

#include "feSpectralLH2d.hpp"
#include "feDynamicField2d.hpp"
#include "feDynamicField2dGlobalManip.hpp"

#include "printMesh.hpp"
#include <map>

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  
  UInt pid = world.rank();
  
  typedef linearQuad     GEOSHAPE;
  typedef pMapItemShare  PMAPTYPE;  
  
  string meshFile  = "../tests/morganaMeshes/mignonQuad2dB.unv";
  string colorFile = "../tests/morganaMeshes/mignonQuad2dB_color.unv";
  
  meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.femap_to_stdA(meshFile, colorFile, false);
  
  
  //Download data
  typedef meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH2D     MESH2D;
  typedef meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT2D  CONNECT2D;
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  
  //Initialize the mapper
  typedef Real                    DOFTYPE;
  typedef feSpectralLH2d<PMAPTYPE> FETYPE;
  
  dofMapDynamic2d_options mapOptions;
  feDynamicField2d<FETYPE,DOFTYPE,dmd2d_vectMajor,dmd2d_standard> field;
  
  field.setCommunicator(world);
  field.setGeometry(grid2d,connectGrid2d);
  field.setOptions(mapOptions);
  

  //Load the feCards
  typedef typename FETYPE::FECARD FECARD;
  FECARD FeCard;
  
  if(pid == 0)
  {
    FeCard.setIsActive(true);
    FeCard.setR(2,1);
    field.setFeCardL(1,FeCard);
  
    FeCard.setIsActive(true);
    FeCard.setR(1,2);
    field.setFeCardL(2,FeCard);
    
    FeCard.setIsActive(true);
    FeCard.setR(2,1);
    field.setFeCardL(3,FeCard);
  }
  else
  {
    FeCard.setIsActive(true);
    FeCard.setR(2,1);
    field.setFeCardL(1,FeCard);
  
    FeCard.setIsActive(true);
    FeCard.setR(1,2);
    field.setFeCardL(2,FeCard);
    
     FeCard.setIsActive(true);
    FeCard.setR(1,2);
    field.setFeCardL(3,FeCard);
  }
  
  //Startup
  field.startup();
  
  //Dofs loading
  if(world.rank() == 0)
  {
    for(UInt i=1; i <= 6; ++i)
    {
      field.setDofL(3 * i -2, 1.0);
      field.setDofL(3 * i -1, 2.0);
      field.setDofL(3 * i   , 3.0);
    }
  }
  else
  {
    for(UInt i=1; i <= 6; ++i)
    {
      field.setDofL(3 * i -2, 3.0);
      field.setDofL(3 * i -1, 4.0);
      field.setDofL(3 * i   , 2.0);
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
  feDynamicField2dGlobalManip<FETYPE,DOFTYPE,dmd2d_vectMajor,dmd2d_standard> tester;
  tester.setCommunicator(world);
  
  //Create the communication pattern
  typedef pMap<pMapItemSendRecv> SENDRECV;
  SENDRECV mapSend, mapRecv;
  tester.createSendRecvMap(field,mapSend,mapRecv);
  
  //Update data
  tester.updateData(field,mapSend);
  
  
  //Evaluation
  if(pid == 1)
  {    
    Real val;
    point3d Y(0.5,0.5,0.0);
    UInt elG = 2;
    
    field.evalG(elG,Y,val);
    
    cout << "Val : " << val << endl;
  }
}
