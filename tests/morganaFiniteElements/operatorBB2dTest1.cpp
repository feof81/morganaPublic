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

#include "dofMapStatic2d.hpp"
#include "elCardFeeder2d.hpp"
#include "feStaticField1d.hpp"
#include "feStaticField2d.hpp"
#include "feStaticField2dGlobalManip.hpp"

#include "opBB2dA.hpp"


using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  UInt pid = world.rank();

  assert(world.size() == 2);
  
  typedef pMapItemShare        PMAPTYPE;
  typedef fePr2d<1,PMAPTYPE>   FIELD_FETYPE;
  typedef Real                 FIELD_DOFTYPE;
  typedef fePr2d<1,PMAPTYPE>   TEST_FETYPE;
  typedef Real                 TEST_DOFTYPE;

  typedef typename FIELD_FETYPE::GEOSHAPE  GEOSHAPE2D;

  
  //Loading
  typedef meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> INIT;
  
  string meshFile = "./tests/morganaMeshes/mignon2dD.msh";
  
  INIT init(world);
  init.gmMesh_to_stdA(meshFile, false);
  
  
  //Testing
  typedef feStaticField2d<FIELD_FETYPE,FIELD_DOFTYPE,dms2d_vectMajor,dms2d_allMode>  FIELD;
  typedef feStaticField2d<TEST_FETYPE,TEST_DOFTYPE,dms2d_vectMajor,dms2d_allMode>    TEST;
  
  typedef typename INIT::MESH2D     MESH2D;
  typedef typename INIT::CONNECT2D  CONNECT2D;
  
  typedef typename INIT::MESH1D      MESH1D;
  typedef typename INIT::CONNECT1D   CONNECT1D;
  
  typedef opBB2dA<FIELD,TEST>        OPERATOR;
  
  typedef OPERATOR::FIELD_OPTIONS    FIELD_OPTIONS;
  typedef OPERATOR::TEST_OPTIONS     TEST_OPTIONS;
  
  
  //Download data
  RCP<MESH1D>           grid1d = init.getGrid1d();
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT1D> connectGrid1d = init.getConnectGrid1d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  //Operator
  FIELD_OPTIONS fieldOptions;
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<OPERATOR> operatore(new OPERATOR(world,*grid2d,*connectGrid2d,*grid1d,*connectGrid1d));
  
  operatore->setOptions_field(fieldOptions);
  operatore->setOptions_test(testOptions);
  operatore->startup();
  
  if(pid == 0)
  {
    //cout << "Elements" << endl;
    //cout << grid1d->getElements() << endl;
    
    UInt el    = 2;
    UInt locEd = 1;
    UInt  I=1, J=1;
    
    point3d Yf(0.5, 0.0, 0.0);
    
    set<UInt> geoIds2d;
    geoIds2d.insert(1);
    
    operatore->setGeoIds(geoIds2d);
    operatore->setElement(el);
    operatore->setLocEdge(locEd);
    operatore->setIJ(I,J,I,J);
    operatore->setCoeff(2.0);
    
    UInt numIndexField = operatore->numIndex_field();
    UInt numIndexTest  = operatore->numIndex_test();
    
    sVect<UInt> indicesField(numIndexField);
    sVect<UInt> indicesTest(numIndexTest);
    sVect<Real> mat(numIndexField * numIndexTest), valTest(numIndexTest), valField(numIndexField);
    
    operatore->indexG_field(indicesField);
    operatore->indexG_test(indicesTest);
    operatore->eval(Yf,mat);
    point3d N  = operatore->computeNormal(Yf);
    point3d Ys = operatore->mapSurfaceY(Yf);
    point3d Yv = operatore->mapVolumeY(Yf);
    
    operatore->eval_field(Yf,valField);
    operatore->eval_test(Yf,valTest);
    
    cout << "Is Boundary : " << operatore->isBoundary() << endl << endl;
    cout << "List Field  : " << endl << operatore->getListMap_field() << endl;
    cout << "List Test   : " << endl << operatore->getListMap_test() << endl;
    
    cout << "Indices field : " << endl << indicesField << endl;
    cout << "Indices test  : " << endl << indicesTest  << endl;
    
    cout << "valTest  : " << endl << valTest << endl;
    cout << "valField : " << endl << valField << endl;
    
    cout << "Matrix        : " << endl << mat << endl;
    cout << "N  : " << N;
    cout << "Ys : " << Ys;
    cout << "Yv : " << Yv;
  }
}
