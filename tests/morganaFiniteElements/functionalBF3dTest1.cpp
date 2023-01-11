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

#include "fePr2d.hpp"
#include "feRt0LT3d.hpp"
#include "feRt0LT3d_extern.hpp"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"

#include "feStaticField2d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField2dGlobalManip.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "functionalBF3d.hpp"
#include "fnBF3dA.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef pMapItemShare        PMAPTYPE;
  typedef feRt0LT3d<PMAPTYPE>  TEST_FETYPE;
  typedef Real                 FIELD_DOFTYPE;
  typedef fePr2d<1,PMAPTYPE>   COEFF_FETYPE;
  typedef Real                 COEFF_DOFTYPE;
  
  typedef typename TEST_FETYPE::FECARD   TEST_FECARD;
  typedef typename TEST_FETYPE::GEOSHAPE GEOSHAPE;
  
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> INIT;

  
  //Loading
  string meshFile = "./tests/morganaMeshes/mignon3dB.msh";
  
  meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  
  
  //Typedefs
  typedef feStaticField3d<TEST_FETYPE,FIELD_DOFTYPE,dms3d_vectMajor,dms3d_allMode>  TEST;
  
  typedef typename INIT::MESH2D     MESH2D;
  typedef typename INIT::MESH3D     MESH3D;
  typedef typename INIT::CONNECT2D  CONNECT2D;
  typedef typename INIT::CONNECT3D  CONNECT3D;  
  
  //Download info
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  

  //Set the coeff field
  typedef feStaticField2d<COEFF_FETYPE,COEFF_DOFTYPE,dms2d_vectMajor,dms2d_allMode> COEFF;
  
  dofMapStatic2d_options mapOptions;
  
  COEFF field2d;
  field2d.setCommunicator(world);
  field2d.setGeometry(grid2d,connectGrid2d);
  field2d.setOptions(mapOptions);
  field2d.startup();
  
  //Field initialization
  feStaticField2dGlobalManip<COEFF_FETYPE,COEFF_DOFTYPE,dms2d_vectMajor,dms2d_allMode> manipulator(world);
  
  sArray<string> functions(1,1);
  functions(1,1) = "2";
  
  manipulator.initilize(functions,field2d);
  
  
  
  //Extern FE
  typedef feRt0LT3d_extern<PMAPTYPE>  EXTERN;
  typedef typename EXTERN::FECARDS    FECARDS;
  
  EXTERN feeder(grid3d,connectGrid3d);
  feeder.setCommDev(world);
  FECARDS feCards = feeder.buildFeCards();
  
  
  //Operator startup
  typedef fnBF3dA<TEST,COEFF>       FUNCTIONAL;
  typedef FUNCTIONAL::TEST_OPTIONS  TEST_OPTIONS;
  typedef FUNCTIONAL::TEST_OUTTYPE  TEST_OUTTYPE;
  
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<FUNCTIONAL> funzionale(new FUNCTIONAL(world,*grid3d,*connectGrid3d,*grid2d,*connectGrid2d));
  funzionale->setOptions_test(testOptions);
  funzionale->setFeCards_test(feCards);
  funzionale->setCoeff(field2d);
  funzionale->startup();
  
  
   
  if(world.rank() == 1)
  {
    //Eval
    UInt el = 1;
    UInt fc = 3;
    point3d Ys(0.0, 1.0, 0.0);
  
    funzionale->setElement(el);
    funzionale->setIJ(1,1);
    funzionale->setLocFace(fc);
    
    if(funzionale->isBoundary())
    {
      cout << "IsBoundary" << endl;
      cout << "Num boundary : " << funzionale->getNumBoundary() << endl;
      cout << "El2d         : " << funzionale->getEl2d() << endl;
      cout << "ElNodes      : " << endl << grid3d->getElementNodesL(el) << endl;
      
      UInt numTest = funzionale->numIndex_test();
      
      sVect<UInt> indicesTest(numTest);
      sVect<TEST_OUTTYPE> testVal(numTest);
      sVect<Real> vect(numTest);
      
      point3d N  = funzionale->mapVolumeY(Ys);
      point3d Yv = funzionale->computeNormal(Ys);
      
      funzionale->indexG_test(indicesTest);
      funzionale->eval_test(Ys,testVal);
      funzionale->eval(Ys,vect);
      
      
      cout << "GEOMETRY - - - - - - " << endl;
      cout << grid3d->getElements() << endl;
      
      cout << "DATA - - - - - - - - " << endl;
      cout << "Num dofsL test  : " << funzionale->getNumDofsL_test() << endl;
      cout << "Num dofsG test  : " << funzionale->getNumDofsG_test() << endl << endl;
  
      cout << "MAPS - - - - - - - - " << endl;
      cout << "Map test  : " << endl << funzionale->getListMap_test() << endl;
    
      cout << "MAPPINGS - - - - - - " << endl;
      cout << "Yv : " << funzionale->mapVolumeY(Ys);
      cout << "N  : " << funzionale->computeNormal(Ys) << endl;
      cout << "Face nodes : " << endl << funzionale->getGlobFacePoints() << endl;
      
      cout << "INDICES - - - - - - -" << endl;
      cout << "Test  : " << endl << indicesTest  << endl;
    
      cout << "EVAL - - - - - - - - " << endl;
      cout << "Test  : " << endl << testVal  << endl;
      cout << "Vect  : " << endl << vect     << endl;
    }
    else
    {
      cout << "NotBoundary" << endl;
    }
  }


}
