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

#include "fePr1d.hpp"
#include "fePr2d.hpp"
#include "feRt0LT2d.hpp"

#include "dofMapStatic2d.hpp"
#include "elCardFeeder2d.hpp"

#include "feStaticField1d.hpp"
#include "feStaticField2d.hpp"
#include "feStaticField1dGlobalManip.hpp"
#include "feStaticField2dGlobalManip.hpp"

#include "operatorLA2d.hpp"
#include "fnLA2dA.hpp"


using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef pMapItemShare        PMAPTYPE;
  typedef fePr1d<1,PMAPTYPE>   TEST_FETYPE;
  typedef Real                 TEST_DOFTYPE;
  typedef fePr1d<1,PMAPTYPE>   COEFF_FETYPE;
  typedef point3d              COEFF_DOFTYPE;

  typedef                 linearTriangle  GEOSHAPE2D;
  typedef typename TEST_FETYPE::GEOSHAPE  GEOSHAPE1D;

  typedef meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> INIT;

  
  //Loading
  string meshFile = "./tests/morganaMeshes/mignon2dB.msh";
  
  INIT init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  
  //Typedefs
  typedef feStaticField1d<TEST_FETYPE,TEST_DOFTYPE,dms1d_vectMajor,dms1d_allMode>  TEST;
  
  typedef typename INIT::MESH1D     MESH1D;
  typedef typename INIT::MESH2D     MESH2D;
  typedef typename INIT::CONNECT1D  CONNECT1D;
  typedef typename INIT::CONNECT2D  CONNECT2D;  
  
  //Download info
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<MESH1D>           grid1d = init.getGrid1d();
  RCP<CONNECT1D> connectGrid1d = init.getConnectGrid1d();
  
 
  //Set the coeff field
  typedef feStaticField1d<COEFF_FETYPE,COEFF_DOFTYPE,dms1d_vectMajor,dms1d_geoIdMode> COEFF;
  
  dofMapStatic1d_options mapOptions;
  mapOptions.addGeoId(1);
  mapOptions.addGeoId(2);
  
  COEFF field1d;
  field1d.setCommunicator(world);
  field1d.setGeometry(grid1d,connectGrid1d);
  field1d.setOptions(mapOptions);
  field1d.startup();
  
  //Field initialization
  feStaticField1dGlobalManip<COEFF_FETYPE,COEFF_DOFTYPE,dms1d_vectMajor,dms1d_geoIdMode> manipulator(world);
  
  sArray<string> functions(3,1);
  functions(1,1) = "2 * x";
  functions(2,1) = "0";
  functions(3,1) = "0";
  
  manipulator.initilize(functions,field1d);
  
  
  //Operator startup
  typedef fnLA2dA<TEST,GEOSHAPE2D,COEFF>  FUNCTIONAL;
  typedef FUNCTIONAL::TEST_OPTIONS        TEST_OPTIONS;
  typedef FUNCTIONAL::TEST_OUTTYPE        TEST_OUTTYPE;
  
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<FUNCTIONAL> funzionale(new FUNCTIONAL(world,*grid2d,*connectGrid2d,*grid1d,*connectGrid1d));
  funzionale->setOptions_test(testOptions);
  funzionale->setCoeff(field1d);
  funzionale->startup();
  
  if(world.rank() == 1)
  {
    //Grid
    cout << "Grid" << endl;
    cout << grid1d->getElements() << endl;
    
    //Eval
    UInt el = 2;
    UInt ed = 2;
    point3d Ys(0.5, 0.0, 0.0);
  
    funzionale->setElement(el);
    funzionale->setIJ(1,1);
    funzionale->setLocEdge(ed);
    
    if(funzionale->isBoundary())
    {
      cout << "IsBoundary" << endl;
      cout << "Num boundary : " << funzionale->getNumBoundary() << endl;
      cout << "El2d         : " << funzionale->getEl1d() << endl;
      cout << "ElNodes      : " << endl << grid2d->getElementNodesL(el) << endl;
      
      UInt numTest = funzionale->numIndex_test();
      
      sVect<UInt> indicesTest(numTest);
      sVect<TEST_OUTTYPE> testVal(numTest);
      sVect<Real> vect(numTest);
      
      point3d N  = funzionale->mapVolumeY(Ys);
      point3d Yv = funzionale->computeNormal(Ys);
      
      funzionale->indexG_test(indicesTest);
      funzionale->eval_test(Ys,testVal);
      funzionale->eval(Ys,vect);
      
      cout << "DATA - - - - - - - - " << endl;
      cout << "Num dofsL test  : " << funzionale->getNumDofsL_test() << endl;
      cout << "Num dofsG test  : " << funzionale->getNumDofsG_test() << endl << endl;
  
      cout << "MAPS - - - - - - - - " << endl;
      cout << "Map test  : " << endl << funzionale->getListMap_test() << endl;
    
      cout << "MAPPINGS - - - - - - " << endl;
      cout << "Yv : " << funzionale->mapVolumeY(Ys);
      cout << "N  : " << funzionale->computeNormal(Ys) << endl;
      cout << "Face nodes : " << endl << funzionale->getGlobEdgePoints() << endl;
      
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
