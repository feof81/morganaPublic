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

#include "feQr2d.hpp"
#include "dofMapStatic2d.hpp"

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
  init.femap_to_stdB(meshFile, colorFile,false);
  
  
  //Print mesh --------------------------------------------------------------------------
  printMesh<GEOSHAPE,PMAPTYPE,PMAPTYPE> printer;
  
  string local = "localMesh" + num2str(world.rank()) + "pid.inp";
  printer.paraviewLocal(local, world.rank(), *init.getGrid2d());
  
  
  //Map test-----------------------------------------------------------------------------
  typedef point3d            DOFTYPE;
  typedef feQr2d<1,PMAPTYPE> FETYPE;
  
  dofMapStatic2d_options mapOptions;
  mapOptions.addGeoId(1);
  
  dofMapStatic2d<FETYPE,DOFTYPE,dms2d_vectMajor,dms2d_geoIdMode> dofMap;
  
  dofMap.setCommunicator(world);
  dofMap.setGeometry(init.getGrid2d(), init.getConnectGrid2d());
  dofMap.setOptions(mapOptions);
  dofMap.startup();

  
  //The map------------------------------------------------------------------------------
  feStaticDofCard2d card;
  card.setGeoType(VERTEX);
  card.setLevel(1);
  card.setLocalId(4);
  card.setLocalElId(2);
  
  UInt I = 1, J = 1;
  
  
  world.barrier();
  if(world.rank() == 0)
  {
    cout << "Dof Map Data" << endl;
    cout << dofMap << endl;
    
    cout << "dofL  : " << dofMap.mapDofL(card) << endl;
    cout << "dofG  : " << dofMap.mapDofG(card) << endl;
    cout << "listL : " << dofMap.mapListL(I,J,card) << endl;
    cout << "listG : " << dofMap.mapListG(I,J,card) << endl;
  }
}
