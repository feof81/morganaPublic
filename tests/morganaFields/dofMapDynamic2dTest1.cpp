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
#include "dofMapDynamic2d.hpp"

#include "printMesh.hpp"
#include <map>

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef linearQuad     GEOSHAPE;
  typedef pMapItemShare  PMAPTYPE;  
  
  string meshFile  = "../tests/morganaMeshes/mignonQuad2dB.unv";
  string colorFile = "../tests/morganaMeshes/mignonQuad2dB_color.unv";
  
  meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.femap_to_stdB(meshFile, colorFile, false);
  
  
  //Initialize the mapper
  typedef point3d                  DOFTYPE;
  typedef feSpectralLH2d<PMAPTYPE> FETYPE;
  
  dofMapDynamic2d_options mapOptions;
  
  dofMapDynamic2d<FETYPE,DOFTYPE,dmd2d_vectMajor,dmd2d_standard> dofMap;
  
  dofMap.setCommunicator(world);
  dofMap.setGeometry(init.getGrid2d(), init.getConnectGrid2d());
  dofMap.setOptions(mapOptions);
  
  
  //Load the feCards
  typedef typename FETYPE::FECARD FECARD;
  FECARD FeCard;
  
  FeCard.setIsActive(true);
  FeCard.setR(2,1);
  dofMap.setFeCardL(1,FeCard);
  
  FeCard.setIsActive(true);
  FeCard.setR(1,2);
  dofMap.setFeCardL(2,FeCard);
  
  //Startup
  dofMap.startup();

  
  world.barrier();
  if(world.rank() == 0)
  {
    cout << "Dof Map Data" << endl;
    cout << dofMap << endl;
    
    //The map
    feDynamicDofCard2d card;
    card.setGeoType(VOLUME);
    card.setLevel(2);
    card.setLocalId(1);
    card.setLocalElId(2);
  
    UInt I = 1, J = 1;
    
    cout << "dofL  : " << dofMap.mapDofL(card) << endl;
    cout << "dofG  : " << dofMap.mapDofG(card) << endl;
    cout << "listL : " << dofMap.mapListL(I,J,card) << endl;
    cout << "listG : " << dofMap.mapListG(I,J,card) << endl;
  }
  sleep(1);
  
  world.barrier();
  if(world.rank() == 1)
  {
    cout << "Dof Map Data" << endl;
    cout << dofMap << endl;
    
    //The map
    feDynamicDofCard2d card;
    card.setGeoType(VOLUME);
    card.setLevel(2);
    card.setLocalId(1);
    card.setLocalElId(2);
  
    UInt I = 1, J = 1;
    
    cout << "dofL  : " << dofMap.mapDofL(card) << endl;
    cout << "dofG  : " << dofMap.mapDofG(card) << endl;
    cout << "listL : " << dofMap.mapListL(I,J,card) << endl;
    cout << "listG : " << dofMap.mapListG(I,J,card) << endl;
  }
  sleep(1);
}
