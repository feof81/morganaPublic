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

#include "komplex.h"
#include "pMapItem.h"
#include "pMapItemShare.h"

#include "geoShapes.h"
#include "meshInit3d.hpp"

#include "feOscHYB3d.hpp"

#include "dofMapStatic3d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "feOscHYB3d_extern.hpp"
#include "elCardFeeder3d.hpp"


using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef pMapItemShare             PMAPTYPE;
  typedef feOscHYB3d<PMAPTYPE>  FIELD_FETYPE;
  typedef Real                 FIELD_DOFTYPE;

  typedef typename FIELD_FETYPE::GEOSHAPE          GEOSHAPE3D;
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> MESHINIT;
  typedef typename MESHINIT::MESH2D                MESH2D;
  typedef typename MESHINIT::MESH3D                MESH3D;
  typedef typename MESHINIT::CONNECT2D             CONNECT2D;
  typedef typename MESHINIT::CONNECT3D             CONNECT3D;

  assert(world.size() == 2);
  
  //Loading
  string meshFile = "./tests/morganaMeshes/mignon3dC.msh";
  
  MESHINIT init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Download data
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //Testing extern
  typedef feOscHYB3d_extern<PMAPTYPE>  EXTERN;
  typedef typename EXTERN::FECARDS    FECARDS;
  typedef elCardFeeder3d<GEOSHAPE3D,PMAPTYPE> ELCARDFEEDER;
  
  feOscHYB3d_extern<PMAPTYPE> feeder(grid3d,connectGrid3d,grid2d,connectGrid2d);
  feeder.setCommDev(world);
  
  set<UInt> geoIds3d, geoIds2d;
  geoIds2d.insert(1);
  geoIds2d.insert(2);
  geoIds3d.insert(1);
  geoIds3d.insert(2);
  
  FECARDS feCards = feeder.buildFeCards(geoIds3d,geoIds2d);
  
  //Element feeder
  ELCARDFEEDER elFeeder(grid3d,connectGrid3d);
  
  //Finite element  
  if(world.rank() == 0)
  {
    cout << feCards << endl;
    
    UInt el = 1;
    point3d Y(0.25, 0.25, 0.0);
  
    feOscHYB3d<PMAPTYPE> finiteElement;
    finiteElement.setCards(feCards(el), elFeeder.getCardLocal(el));
    
    UInt num = finiteElement.getNumBasis();
    sVect<komplex> val(num);
    
    finiteElement.globalEval(Y,val);
    
    cout << "Num basis : " << num << endl;
    cout << "Dofcards  : " << endl << finiteElement.getDofCards();
    cout << "FeCard    : " << endl << feCards(el);
    cout << "Val       : " << endl << val << endl;
  }
}
