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

#include "feHY2d.hpp"
#include "feHY2d_extern.hpp"
#include "feRt0Loc2d.hpp"

#include "dofMapStatic2d.hpp"
#include "elCardFeeder2d.hpp"

#include "feDynamicField2d.hpp"
#include "feStaticField2d.hpp"
#include "feStaticField2dGlobalManip.hpp"

#include "localMatrixHY2d.hpp"
#include "opHY2dA.hpp"
//#include "integratorHYB2d.hpp"

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
  typedef feRt0Loc2d<PMAPTYPE>  FIELD_FETYPE;
  typedef Real                  FIELD_DOFTYPE;
  typedef feHY2d<PMAPTYPE>      TEST_FETYPE;
  typedef Real                  TEST_DOFTYPE;
  
  typedef typename FIELD_FETYPE::FECARD FIELD_FECARD;
  typedef typename TEST_FETYPE::FECARD  TEST_FECARD;

  typedef typename FIELD_FETYPE::GEOSHAPE GEOSHAPE;

  
  //Loading
  string meshFile = "./tests/morganaMeshes/mignon2dB.msh";
  
  meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  
  //Typedefs
  typedef feStaticField2d<FIELD_FETYPE,FIELD_DOFTYPE,dms2d_vectMajor,dms2d_allMode>  FIELD;
  typedef feDynamicField2d<TEST_FETYPE,TEST_DOFTYPE,dmd2d_vectMajor,dmd2d_standard>  TEST;
  
  typedef typename FIELD::MESH2D     MESH2D;
  typedef typename FIELD::CONNECT2D  CONNECT2D;
  typedef opHY2dA<FIELD,TEST>        OPERATOR;
  
  typedef OPERATOR::FIELD_OPTIONS FIELD_OPTIONS;
  typedef OPERATOR::TEST_OPTIONS  TEST_OPTIONS;
  typedef OPERATOR::FIELD_OUTTYPE FIELD_OUTTYPE;
  typedef OPERATOR::TEST_OUTTYPE  TEST_OUTTYPE;
  
  
  //Download info
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  
  //Test function feCards
  typedef feHY2d_extern<PMAPTYPE>  EXTERN;
  typedef typename EXTERN::FECARDS FECARDS;
  
  feHY2d_extern<PMAPTYPE> feeder(grid2d,connectGrid2d);
  feeder.setCommDev(world);
  
  FECARDS feCards = feeder.buildFeCards();
  
  
  //Operator startup
  FIELD_OPTIONS fieldOptions;
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<OPERATOR> operatore(new OPERATOR(world,*grid2d,*connectGrid2d));
  
  operatore->setFeCards_test(feCards);
  operatore->setOptions_field(fieldOptions);
  operatore->setOptions_test(testOptions);
  
  operatore->startup();
  
   
  if(world.rank() == 1)
  {
    //Eval
    UInt el = 1;
    UInt ed = 1;
    point3d Ys(0.5, 0.0, 0.0);
  
    operatore->setElement(el);
    operatore->setIJ(1,1,1,1);
    operatore->setLocEdge(ed);
    
    UInt numField = operatore->numIndex_field();
    UInt numTest  = operatore->numIndex_test();
    
    sVect<UInt> indicesField(numField), indicesTest(numTest);
    sVect<FIELD_OUTTYPE> fieldVal(numField);
    sVect<TEST_OUTTYPE>  testVal(numTest);
    sVect<Real> mat(numField * numTest);
   
    operatore->indexG_field(indicesField);
    operatore->indexG_test(indicesTest);
    
    operatore->eval_field(Ys,fieldVal);
    operatore->eval_test(Ys,testVal);
    operatore->eval(Ys,mat);
    
    
    cout << "GEOMETRY - - - - - - " << endl;
    cout << grid2d->getEdges() << endl;
    
    cout << "DATA - - - - - - - - " << endl;
    cout << "Num dofsL field : " << operatore->getNumDofsL_field() << endl;
    cout << "Num dofsG field : " << operatore->getNumDofsG_field() << endl;
    cout << "Num dofsL test  : " << operatore->getNumDofsL_test() << endl;
    cout << "Num dofsG test  : " << operatore->getNumDofsG_test() << endl << endl;
  
    cout << "MAPS - - - - - - - - " << endl;
    cout << "Map field : " << endl << operatore->getListMap_field() << endl;
    cout << "Map test  : " << endl << operatore->getListMap_test() << endl;
    
    cout << "MAPPINGS - - - - - - " << endl;
    cout << "Yv : " << operatore->mapVolumeY(Ys);
    cout << "N  : " << operatore->computeNormal(Ys) << endl;
    cout << "Face nodes : " << endl << operatore->getGlobEdgePoints() << endl;
    
    cout << "INDICES - - - - - - -" << endl;
    cout << "Field : " << endl << indicesField << endl;
    cout << "Test  : " << endl << indicesTest  << endl;
    
    cout << "EVAL - - - - - - - - " << endl;
    cout << "Field : " << endl << fieldVal << endl;
    cout << "Test  : " << endl << testVal  << endl;
    
    cout << "MATRIX - - - - - - - " << endl;
    cout << mat << endl;
  }
}
