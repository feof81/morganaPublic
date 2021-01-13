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

#include "feSpectralLH1d.hpp"
#include "feDynamicField1d.hpp"
#include "feDynamicField1dGlobalManip.hpp"

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
  
  typedef linearLine     GEOSHAPE1D;
  typedef linearQuad     GEOSHAPE2D;
  typedef pMapItemShare  PMAPTYPE;  
  
  string meshFile  = "../tests/morganaMeshes/mignonQuad2dB.unv";
  string colorFile = "../tests/morganaMeshes/mignonQuad2dB_color.unv";
  
  meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> init(world);
  init.femap_to_stdA(meshFile, colorFile, false);
  
  
  //Download data
  typedef meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>::MESH1D     MESH1D;
  typedef meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>::CONNECT1D  CONNECT1D;
  
  RCP<MESH1D>           grid1d = init.getGrid1d();
  RCP<CONNECT1D> connectGrid1d = init.getConnectGrid1d();
  
  
//Initialize the mapper
  typedef Real                    DOFTYPE;
  typedef feSpectralLH1d<PMAPTYPE> FETYPE;
  
  dofMapDynamic1d_options mapOptions;
  feDynamicField1d<FETYPE,DOFTYPE,dmd1d_vectMajor,dmd1d_standard> field;
  
  field.setCommunicator(world);
  field.setGeometry(grid1d,connectGrid1d);
  field.setOptions(mapOptions); 
  
  
  //Load the feCards
  typedef typename FETYPE::FECARD FECARD;
  FECARD FeCard;
  FeCard.setR(1);
  
  for(UInt i=1; i <= grid1d->getNumElements(); ++i)
  {
    field.setFeCardL(i,FeCard);
  }
  
  //Startup
  field.startup();
  
  //Dofs loading
  for(UInt i=1; i <= 2 * grid1d->getNumElements(); ++i)
  {
    field.setDofL(i,4.0);
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
  feDynamicField1dGlobalManip<FETYPE,DOFTYPE,dmd1d_vectMajor,dmd1d_standard> tester;
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
    point3d Y(0.5,0.0,0.0);
    UInt elG = 9;
    
    field.evalG(elG,Y,val);
    
    cout << "Val : " << val << endl;
  }
}
