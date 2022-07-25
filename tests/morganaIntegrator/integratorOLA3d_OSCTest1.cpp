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

#include "feOscPr2d.hpp"
#include "feOscPr3d.hpp"
#include "feOscPr2d_extern.h"
#include "feOscPr3d_extern.h"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"
#include "feStaticField2d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "opLA3dA_CX.hpp"
#include "traitsIntegratorOLA3d_OSC.hpp"

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef pMapItem                  PMAPTYPE;
  typedef feOscPr3d<0,1,PMAPTYPE>   FIELD_FETYPE;
  typedef komplex                   FIELD_DOFTYPE;
  typedef feOscPr2d<0,0,PMAPTYPE>   TEST_FETYPE;
  typedef komplex                   TEST_DOFTYPE;

  typedef typename FIELD_FETYPE::GEOSHAPE GEOSHAPE3D;
  typedef typename TEST_FETYPE::GEOSHAPE  GEOSHAPE2D;

  //Loading
  string meshFile = "./tests/morganaMeshes/mignon3dC.msh";
  
  meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Typedefs
  typedef feStaticField3d<FIELD_FETYPE,FIELD_DOFTYPE,dms3d_vectMajor,dms3d_allMode> FIELD;
  typedef feStaticField2d<TEST_FETYPE,TEST_DOFTYPE,dms2d_vectMajor,dms2d_allMode>   TEST;
  
  typedef typename TEST::MESH2D     MESH2D;
  typedef typename FIELD::MESH3D    MESH3D;
  typedef typename TEST::CONNECT2D  CONNECT2D;
  typedef typename FIELD::CONNECT3D CONNECT3D;
  typedef opLA3dA_CX<FIELD,TEST>    OPERATOR;
  
  typedef OPERATOR::FIELD_OPTIONS FIELD_OPTIONS;
  typedef OPERATOR::TEST_OPTIONS  TEST_OPTIONS;
  
  //Download info
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //Build the feCards
  typedef feOscPr3d_extern<0,1,PMAPTYPE>::FECARDS FECARDS;
  
  Real K = 2.5;
  
  feOscPr2d_extern<0,0,PMAPTYPE> feExtern2d;
  feExtern2d.setGeometry(grid2d,connectGrid2d);
  feExtern2d.setCommDev(world);
  FECARDS feCards2d = feExtern2d.buildFeCards(K);
  
  feOscPr3d_extern<0,1,PMAPTYPE> feExtern3d;
  feExtern3d.setGeometry(grid3d,connectGrid3d);
  feExtern3d.setCommDev(world);
  FECARDS feCards3d = feExtern3d.buildFeCards(K);
  
  //Operator startup
  FIELD_OPTIONS fieldOptions;
  TEST_OPTIONS  testOptions;
  
  OPERATOR operatore(world,*grid3d,*connectGrid3d,*grid2d,*connectGrid2d);
  operatore.setOptions_field(fieldOptions);
  operatore.setOptions_test(testOptions);
  operatore.setFeCards_field(feCards3d);
  operatore.setFeCards_test(feCards2d);
  operatore.startup();
  
  //Eval
  UInt el = 1;
  UInt fc = 2;
  
  operatore.setElement(el);
  operatore.setLocFace(fc);
  operatore.setIJ(1,1,1,1);
  
  //Integration 
  traitsIntegratorOLA3d_OSC<OPERATOR,linearTriangle,STANDARD,1> integrator;
  
  UInt pid = world.rank();
  
  if(pid == 0)
  {
    cout << "isBoundary : " << operatore.isBoundary()        << endl;
    cout << "points 2d  : " << operatore.getGlobFacePoints() << endl;
    cout << "points 3d  : " << operatore.getIntegrationGrid()->getElementNodesL(el) << endl;
    
    sVect<Real> mat = integrator.integration(operatore);

    cout << "Matrix" << endl;
    cout << mat << endl;
  }
}
