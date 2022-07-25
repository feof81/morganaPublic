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

#include "feOscPr3d.hpp"
#include "feOscCr3d.hpp"
#include "feOscCr3d_extern.h"

#include "feOscHYB3d.hpp"
#include "feOscHYB3d_extern.hpp"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"
#include "feDynamicField3d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "opHY3dA_CX.hpp"
#include "traitsIntegratorOHY3d_OSC.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef pMapItem              PMAPTYPE;
  typedef feOscCr3d<1,PMAPTYPE> FIELD_FETYPE;
  typedef komplex               FIELD_DOFTYPE;
  typedef feOscHYB3d<PMAPTYPE>  TEST_FETYPE;
  typedef komplex               TEST_DOFTYPE;

  typedef typename FIELD_FETYPE::GEOSHAPE    GEOSHAPE;
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> INIT;

  //Loading
  string meshFile = "./tests/morganaMeshes/mignon3dC.msh";
  
  INIT init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Typedefs
  typedef feStaticField3d<FIELD_FETYPE,FIELD_DOFTYPE,dms3d_vectMajor,dms3d_allMode> FIELD;
  typedef feDynamicField3d<TEST_FETYPE,TEST_DOFTYPE,dmd3d_vectMajor,dmd3d_standard>  TEST;
  
  typedef typename INIT::MESH2D     MESH2D;
  typedef typename INIT::MESH3D     MESH3D;
  typedef typename INIT::CONNECT2D  CONNECT2D;
  typedef typename INIT::CONNECT3D  CONNECT3D;
  typedef opHY3dA_CX<FIELD,TEST>    OPERATOR;
  
  typedef OPERATOR::FIELD_OPTIONS FIELD_OPTIONS;
  typedef OPERATOR::TEST_OPTIONS  TEST_OPTIONS;
  
  //Download info
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //Build field feCards
  typedef feOscCr3d_extern<1,PMAPTYPE>::FECARDS FIELD_FECARDS;
  
  Real K = 2.5;
  
  feOscCr3d_extern<1,PMAPTYPE> fieldExtern;
  fieldExtern.setGeometry(grid3d,connectGrid3d);
  fieldExtern.setCommDev(world);
  FIELD_FECARDS fieldFeCards = fieldExtern.buildFeCards(K);
  
  //Build test feCards
  typedef feOscHYB3d_extern<PMAPTYPE>::FECARDS TEST_FECARDS;
  
  set<UInt> geoIds3d, geoIds2d;
  geoIds2d.insert(1);
  geoIds2d.insert(2);
  geoIds3d.insert(1);
  geoIds3d.insert(2);
  
  feOscHYB3d_extern<PMAPTYPE> feederTest(grid3d,connectGrid3d,grid2d,connectGrid2d);
  feederTest.setCommDev(world);
  TEST_FECARDS testFeCards = feederTest.buildFeCards(geoIds3d,geoIds2d);

  //Operator startup
  FIELD_OPTIONS fieldOptions;
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<OPERATOR> operatore(new OPERATOR(world,*grid3d,*connectGrid3d));
  operatore->setOptions_field(fieldOptions);
  operatore->setOptions_test(testOptions);
  operatore->setFeCards_field(fieldFeCards);
  operatore->setFeCards_test(testFeCards);
  operatore->startup();
  
  //Eval
  UInt el = 1;
  operatore->setElement(el);
  operatore->setIJ(1,1,1,1);
  
  //Integration 
  traitsIntegratorOHY3d_OSC<OPERATOR,linearTriangle,STANDARD,4> integrator;
  
  UInt pid = world.rank();
  
  if(pid == 0)
  {
    sVect<Real> mat = integrator.integration(operatore);

    cout << "element 3d" << endl;
    cout << grid3d->getElementNodesL(el) << endl << endl;
    
    cout << "Matrix" << endl;
    cout << mat << endl;
  }
}
