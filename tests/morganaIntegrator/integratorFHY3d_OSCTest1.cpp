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

#include "feHYB3d.hpp"
#include "feHYB3d_extern.hpp"
#include "feOscHYB3d.hpp"
#include "feOscHYB3d_extern.hpp"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"

#include "feDynamicField3d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "fnHY3dA_CX.hpp"
#include "localVectorHY3d.hpp"
#include "traitsIntegratorFHY3d_OSC.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef pMapItemShare         PMAPTYPE;
  typedef feOscHYB3d<PMAPTYPE>  TEST_FETYPE;
  typedef komplex               TEST_DOFTYPE;
  typedef feHYB3d<PMAPTYPE>     COEFF_FETYPE;
  typedef Real                  COEFF_DOFTYPE;
  
  typedef typename COEFF_FETYPE::FECARD   COEFF_FECARD;
  typedef typename TEST_FETYPE::FECARD    TEST_FECARD;
  typedef typename TEST_FETYPE::GEOSHAPE  GEOSHAPE;
  
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> INIT;

  //Loading----------------------------------------------------------------------------------------
  string meshFile = "./tests/morganaMeshes/mignon3dB.msh";
  
  INIT init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Typedefs---------------------------------------------------------------------------------------
  typedef feDynamicField3d<TEST_FETYPE,  TEST_DOFTYPE,  dmd3d_vectMajor, dmd3d_standard> TEST;
  typedef feDynamicField3d<COEFF_FETYPE, COEFF_DOFTYPE, dmd3d_vectMajor, dmd3d_standard> COEFF;
  
  typedef typename INIT::MESH2D     MESH2D;
  typedef typename INIT::MESH3D     MESH3D;
  typedef typename INIT::CONNECT2D  CONNECT2D;
  typedef typename INIT::CONNECT3D  CONNECT3D;
  typedef fnHY3dA_CX<TEST,COEFF>    FUNCTIONAL;
 
  typedef TEST::OPTIONS  TEST_OPTIONS;
  typedef TEST::OUTTYPE  TEST_OUTTYPE;
  
  typedef COEFF::OPTIONS  COEFF_OPTIONS;
  typedef COEFF::OUTTYPE  COEFF_OUTTYPE;
  
  //Download info----------------------------------------------------------------------------------
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //Test function feCards---------------------------------------------------------------------------
  typedef feOscHYB3d_extern<PMAPTYPE>   TEST_EXTERN;
  typedef typename TEST_EXTERN::FECARDS TEST_FECARDS;
  
  set<UInt> geoIds3d, geoIds2d;
  geoIds2d.insert(1);
  geoIds2d.insert(2);
  geoIds3d.insert(1);
  geoIds3d.insert(2);
  
  feOscHYB3d_extern<PMAPTYPE> testFeeder(grid3d,connectGrid3d,grid2d,connectGrid2d);
  testFeeder.setCommDev(world);
  TEST_FECARDS testFeCards = testFeeder.buildFeCards(geoIds3d,geoIds2d);
  
  //Coeff function feCards--------------------------------------------------------------------------
  typedef feHYB3d_extern<PMAPTYPE>       COEFF_EXTERN;
  typedef typename COEFF_EXTERN::FECARDS COEFF_FECARDS;
  
  feHYB3d_extern<PMAPTYPE> coeffFeeder(grid3d,connectGrid3d,grid2d,connectGrid2d);
  coeffFeeder.setCommDev(world);
  COEFF_FECARDS coeffFeCards = coeffFeeder.buildFeCards(geoIds3d,geoIds2d);
  
  //Coeff------------------------------------------------------------------------------------------
  COEFF_OPTIONS coeffOptions;
  
  COEFF coeffField;
  coeffField.setCommunicator(world);
  coeffField.setGeometry(grid3d,connectGrid3d);
  coeffField.setOptions(coeffOptions);
  coeffField.setFeCards(coeffFeCards);
  coeffField.startup();
  
  if(world.rank() == 0)
  {
    coeffField.setDofL(1, 3.0);
    coeffField.setDofL(2, 1.0);
    coeffField.setDofL(3, 4.0);
  }
  else
  {
    coeffField.setDofL(1, 3.0);
    coeffField.setDofL(2, 2.0);
    coeffField.setDofL(3, 4.0);
  }
  
  //Operator startup-------------------------------------------------------------------------------
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<FUNCTIONAL> functional(new FUNCTIONAL(world,*grid3d,*connectGrid3d));
  functional->setFeCards_test(testFeCards);
  functional->setOptions_test(testOptions);
  functional->setCoeff(coeffField);
  functional->startup();
  
  UInt el = 2;
  UInt fc = 1;
 
  functional->setElement(el);
  functional->setIJ(1,1);
  functional->setLocFace(fc);
  
  //Integration------------------------------------------------------------------------------------
  traitsIntegratorFHY3d_OSC<FUNCTIONAL,linearTriangle,STANDARD,4> integrator;
  
  UInt pid = world.rank();
  
  if(pid == 1)
  {
    sVect<Real> vec = integrator.integration(functional);

    cout << "element 3d" << endl;
    cout << grid3d->getElementNodesL(el) << endl << endl;
    
    cout << "Vector" << endl;
    cout << vec << endl;
  }
  
  //Global vector-----------------------------------------------------------------------------------
  /*typedef localVectorHY3d<TEST,intDefault,STANDARD,1> LOCVECTOR;
  typedef vectorBuilder<LOCVECTOR>                    FEBUILDER;
  typedef typename FEBUILDER::EPETRA_VECTOR           EPETRA_VECTOR;
  
  RCP<EPETRA_VECTOR> vector;
  
  FEBUILDER feVector;
  feVector.setFunctional(functional);
  feVector.setCommDev(world);
  
  feVector.buildEpetraVector(vector);
  cout << *vector << endl;*/
}
