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
#include "fePr3d.hpp"
#include "feHY3d.hpp"
#include "feHY3d_extern.hpp"
#include "feRt0Loc3d.hpp"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"

#include "feDynamicField3d.hpp"
#include "feStaticField2d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField2dGlobalManip.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "opHY3dB.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  typedef pMapItemShare         PMAPTYPE;
  typedef feRt0Loc3d<PMAPTYPE>  FIELD_FETYPE;
  typedef Real                  FIELD_DOFTYPE;
  typedef feHY3d<PMAPTYPE>      TEST_FETYPE;
  typedef Real                  TEST_DOFTYPE;
  
  typedef typename FIELD_FETYPE::FECARD FIELD_FECARD;
  typedef typename TEST_FETYPE::FECARD  TEST_FECARD;

  typedef typename FIELD_FETYPE::GEOSHAPE GEOSHAPE;
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> INIT;
  
  //Loading----------------------------------------------------------------------------------------
  string meshFile = "./tests/morganaMeshes/mignon3dD.msh";
  //string meshFile = "./geometries/cubes3d/testCubeF.msh";
  
  INIT init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Typedefs---------------------------------------------------------------------------------------
  typedef feStaticField3d<FIELD_FETYPE,FIELD_DOFTYPE,dms3d_vectMajor,dms3d_allMode>  FIELD;
  typedef feDynamicField3d<TEST_FETYPE,TEST_DOFTYPE,dmd3d_vectMajor,dmd3d_standard>  TEST;
  
  typedef typename INIT::MESH2D     MESH2D;
  typedef typename INIT::MESH3D     MESH3D;
  typedef typename INIT::CONNECT2D  CONNECT2D;
  typedef typename INIT::CONNECT3D  CONNECT3D;
  
  //Coeff typedefs---------------------------------------------------------------------------------
  typedef feHY3d<PMAPTYPE>    COEFFH_FETYPE;
  typedef Real                COEFFH_DOFTYPE;
  
  typedef fePr3d<0,PMAPTYPE>  COEFFV_FETYPE;
  typedef Real                COEFFV_DOFTYPE;
  
  typedef fePr2d<0,PMAPTYPE>  COEFFS_FETYPE;
  typedef Real                COEFFS_DOFTYPE;
  
  typedef feDynamicField3d<COEFFH_FETYPE, COEFFH_DOFTYPE, dmd3d_vectMajor, dmd3d_standard>  COEFFH;
  typedef feStaticField3d<COEFFV_FETYPE, COEFFV_DOFTYPE, dms3d_vectMajor, dms3d_geoIdMode>  COEFFV;
  typedef feStaticField2d<COEFFS_FETYPE, COEFFS_DOFTYPE, dms2d_vectMajor, dms2d_geoIdMode>  COEFFS;
  
  typedef COEFFH::OPTIONS  COEFFH_OPTIONS;
  typedef COEFFV::OPTIONS  COEFFV_OPTIONS;
  typedef COEFFS::OPTIONS  COEFFS_OPTIONS;
  
  typedef typename COEFFH::DOFCARD DOFCARD_H;
  typedef typename COEFFV::DOFCARD DOFCARD_V;
  typedef typename COEFFS::DOFCARD DOFCARD_S;
  
  //Operator typedefs------------------------------------------------------------------------------
  typedef opHY3dB<FIELD,TEST,COEFFH,COEFFV>  OPERATOR;
  
  typedef OPERATOR::FIELD_OPTIONS FIELD_OPTIONS;
  typedef OPERATOR::TEST_OPTIONS  TEST_OPTIONS;
  typedef OPERATOR::FIELD_OUTTYPE FIELD_OUTTYPE;
  typedef OPERATOR::TEST_OUTTYPE  TEST_OUTTYPE;
  
  
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
  
  //Startup the coeffH-----------------------------------------------------------------------------
  COEFFH_OPTIONS coeffH_options;
  
  RCP<COEFFH> coeffH(new COEFFH);
  coeffH->setCommunicator(world);
  coeffH->setGeometry(grid3d,connectGrid3d);
  coeffH->setOptions(coeffH_options);
  coeffH->setFeCards(feCards);
  coeffH->startup();
  
  //Startup the coeffV-----------------------------------------------------------------------------
  COEFFV_OPTIONS coeffV_options;
  coeffV_options.addGeoId(1);
  coeffV_options.addGeoId(2);
  
  RCP<COEFFV> coeffV(new COEFFV);
  coeffV->setCommunicator(world);
  coeffV->setGeometry(grid3d,connectGrid3d);
  coeffV->setOptions(coeffV_options);
  coeffV->startup();
  
  feStaticField3dGlobalManip<COEFFV_FETYPE,COEFFV_DOFTYPE,dms3d_vectMajor,dms3d_geoIdMode> manipulatorV(world);
  sArray<string> functionsV(1,1);
  functionsV(1,1) = "3";
  manipulatorV.initilize(functionsV,coeffV);
  
  //Startup surface conductivity field-------------------------------------------------------------
  COEFFS_OPTIONS coeffS_options;
  coeffS_options.addGeoId(3);
  
  RCP<COEFFS> coeffS(new COEFFS);
  coeffS->setCommunicator(world);
  coeffS->setGeometry(grid2d,connectGrid2d);
  coeffS->setOptions(coeffS_options);
  coeffS->startup();
  
  feStaticField2dGlobalManip<COEFFS_FETYPE,COEFFS_DOFTYPE,dms2d_vectMajor,dms2d_geoIdMode> manipulatorS(world);
  sArray<string> functionsS(1,1);
  functionsS(1,1) = "2";
  manipulatorS.initilize(functionsS,coeffS);
  
  //Interpolation----------------------------------------------------------------------------------
  DOFCARD_H dofCardH;
  dofCardH.setGeoType(FACE);
  dofCardH.setLevel(1);
  
  Real dof;
  UInt fc, el2d, lid3d;
  point3d B(1.0/3.0, 1.0/3.0, 0.0); 
  
  for(UInt el3d=1; el3d <= grid3d->getNumElements(); el3d++)
  {
    for(UInt f=1; f <= GEOSHAPE::numFaces; ++f)
    {
      fc = connectGrid3d->getElementToFace(el3d,f);
      
      dofCardH.setLocalElId(el3d);
      dofCardH.setLocalId(f);

      if(connectGrid3d->getFaceIsBoundary(fc) && coeffH->getDofMapper().isActive(dofCardH))
      {
	el2d = connectGrid3d->getFaceBFace(fc);
	coeffS->evalL(el2d,B,dof);
	
	lid3d = coeffH->getDofMapper().mapDofL(dofCardH);
	coeffH->setDofL(lid3d,dof);
      }
    }
  }

  //Operator startup
  FIELD_OPTIONS fieldOptions;
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<OPERATOR> operatore(new OPERATOR(world,*grid3d,*connectGrid3d));
  
  operatore->setFeCards_test(feCards);
  operatore->setOptions_field(fieldOptions);
  operatore->setOptions_test(testOptions);
  operatore->setCoeff(coeffH,coeffV);
  operatore->startup();
   
  if(world.rank() == 1)
  {
    //Eval
    UInt el = 1;
    UInt fc = 1;
    point3d Ys(0.25, 0.25, 0.0);
  
    operatore->setElement(el);
    operatore->setIJ(1,1,1,1);
    operatore->setLocFace(fc);
    
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
    cout << grid3d->getElements() << endl;
    
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
    cout << "Face nodes : " << endl << operatore->getGlobFacePoints() << endl;
    
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
