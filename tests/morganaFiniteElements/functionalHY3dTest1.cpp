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

#include "feHY3d.hpp"
#include "feHY3d_extern.hpp"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"

#include "feDynamicField3d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "fnHY3dA.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef pMapItemShare    PMAPTYPE;
  typedef feHY3d<PMAPTYPE> TEST_FETYPE;
  typedef Real             TEST_DOFTYPE;
  
  typedef typename TEST_FETYPE::FECARD   TEST_FECARD;
  typedef typename TEST_FETYPE::GEOSHAPE GEOSHAPE;
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> INIT;

  //Loading----------------------------------------------------------------------------------------
  string meshFile = "./tests/morganaMeshes/mignon3dB.msh";
  
  INIT init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Typedefs---------------------------------------------------------------------------------------
  typedef feDynamicField3d<TEST_FETYPE,TEST_DOFTYPE,dmd3d_vectMajor,dmd3d_standard> TEST;
  typedef feDynamicField3d<TEST_FETYPE,TEST_DOFTYPE,dmd3d_vectMajor,dmd3d_standard> COEFF;
  
  typedef typename INIT::MESH2D     MESH2D;
  typedef typename INIT::MESH3D     MESH3D;
  typedef typename INIT::CONNECT2D  CONNECT2D;
  typedef typename INIT::CONNECT3D  CONNECT3D;
  typedef fnHY3dA<TEST,COEFF>       FUNCTIONAL;
  
  typedef FUNCTIONAL::TEST_OPTIONS  TEST_OPTIONS;
  typedef FUNCTIONAL::TEST_OUTTYPE  TEST_OUTTYPE;
  
  typedef TEST_OPTIONS  COEFF_OPTIONS;
  typedef TEST_OUTTYPE  COEFF_OUTTYPE;
  
  //Download info----------------------------------------------------------------------------------
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //Test function feCards--------------------------------------------------------------------------
  typedef feHY3d_extern<PMAPTYPE>  EXTERN;
  typedef typename EXTERN::FECARDS FECARDS;
  
  feHY3d_extern<PMAPTYPE> feeder(grid3d,connectGrid3d,grid2d,connectGrid2d);
  feeder.setCommDev(world);
  
  set<UInt> geoIds3d, geoIds2d;
  geoIds2d.insert(1);
  geoIds2d.insert(2);
  geoIds3d.insert(1);
  geoIds3d.insert(2);
  
  FECARDS feCards = feeder.buildFeCards(geoIds3d, geoIds2d);
  
  //Coeff------------------------------------------------------------------------------------------
  COEFF_OPTIONS options;
  
  COEFF coeffField;
  coeffField.setCommunicator(world);
  coeffField.setGeometry(grid3d,connectGrid3d);
  coeffField.setOptions(options);
  coeffField.setFeCards(feCards);
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
  functional->setFeCards_test(feCards);
  functional->setOptions_test(testOptions);
  functional->setCoeff(coeffField);
  functional->startup();
  
  if(world.rank() == 1)
  {
    //Eval
    UInt el = 2;
    UInt fc = 4;
    point3d Yf(0.25, 0.25, 0.0);
  
    functional->setElement(el);
    functional->setIJ(1,1);
    functional->setLocFace(fc);
    
    UInt numTest = functional->numIndex_test();
    
    sVect<UInt>          indicesLid(numTest);
    sVect<UInt>          indicesGid(numTest);
    sVect<TEST_OUTTYPE>  testVal(numTest);
    sVect<Real>          vect(numTest);
    
    functional->indexL_test(indicesLid); 
    functional->indexG_test(indicesGid);
    functional->eval_test(Yf,testVal);
    functional->eval(Yf,vect);
    
    cout << "NODES - - - - - - " << endl;
    cout << grid3d->getNodes() << endl;
    
    cout << "GEOMETRY - - - - - - " << endl;
    cout << grid3d->getElements() << endl;
    
    cout << "DATA - - - - - - - - " << endl;
    cout << "Num dofsL test  : " << functional->getNumDofsL_test() << endl;
    cout << "Num dofsG test  : " << functional->getNumDofsG_test() << endl << endl;
    
    cout << "MAPS - - - - - - - - " << endl;
    cout << "Map test  : " << endl << functional->getListMap_test() << endl;
    
    cout << "MAPPINGS - - - - - - " << endl;
    cout << "Yv : " << functional->mapVolumeY(Yf);
    cout << "N  : " << functional->computeNormal(Yf) << endl;
    cout << "Face nodes : " << endl << functional->getGlobFacePoints() << endl;
    
    cout << "INDICES - - - - - - -" << endl;
    cout << "Indices Lid : " << endl << indicesLid << endl;
    cout << "Indices Gid : " << endl << indicesGid << endl;
    
    cout << "EVAL - - - - - - - - " << endl;
    cout << "Test  : " << endl << testVal  << endl;
    
    cout << "VECTOR - - - - - - - " << endl;
    cout << vect << endl;
  }
}
