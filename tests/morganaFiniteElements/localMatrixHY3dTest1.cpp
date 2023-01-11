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
#include "feRt0Loc3d.hpp"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"

#include "feDynamicField3d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "localMatrixHY3d.hpp"
#include "integratorOHY3d_STD.hpp"
#include "opHY3dA.hpp"

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
  typedef feRt0Loc3d<PMAPTYPE>  FIELD_FETYPE;
  typedef Real                  FIELD_DOFTYPE;
  typedef feHY3d<PMAPTYPE>      TEST_FETYPE;
  typedef Real                  TEST_DOFTYPE;
  
  typedef typename FIELD_FETYPE::FECARD FIELD_FECARD;
  typedef typename TEST_FETYPE::FECARD  TEST_FECARD;

  typedef typename FIELD_FETYPE::GEOSHAPE GEOSHAPE;
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> INIT;

  
  //Loading
  string meshFile = "./tests/morganaMeshes/mignon3dE.msh";
  
  meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile);
  
  
  //Typedefs
  typedef feStaticField3d<FIELD_FETYPE,FIELD_DOFTYPE,dms3d_vectMajor,dms3d_allMode>  FIELD;
  typedef feDynamicField3d<TEST_FETYPE,TEST_DOFTYPE,dmd3d_vectMajor,dmd3d_standard>  TEST;
  
  typedef typename INIT::MESH2D     MESH2D;
  typedef typename INIT::MESH3D     MESH3D;
  typedef typename INIT::CONNECT2D  CONNECT2D;
  typedef typename INIT::CONNECT3D  CONNECT3D;
  typedef opHY3dA<FIELD,TEST>       OPERATOR;
  
  typedef OPERATOR::FIELD_OPTIONS FIELD_OPTIONS;
  typedef OPERATOR::TEST_OPTIONS  TEST_OPTIONS;
  typedef OPERATOR::FIELD_OUTTYPE FIELD_OUTTYPE;
  typedef OPERATOR::TEST_OUTTYPE  TEST_OUTTYPE;
  
  
  //Download info
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  
  //Test function feCards
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
  
  //Operator startup
  FIELD_OPTIONS fieldOptions;
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<OPERATOR> operatore(new OPERATOR(world,*grid3d,*connectGrid3d));
  
  operatore->setFeCards_test(feCards);
  operatore->setOptions_field(fieldOptions);
  operatore->setOptions_test(testOptions);
  
  operatore->startup();
  
   
  if(world.rank() == 1)
  {
    //Eval
    UInt el = 2;
    operatore->setElement(el);
    operatore->setIJ(1,1,1,1);
    
    //Integrator 
    typedef localMatrixHY3d<FIELD,TEST> LOCALMATRIX;
    LOCALMATRIX localMatrix(operatore);
    
    UInt numRow = localMatrix.numIndex_row();
    UInt numCol = localMatrix.numIndex_col();
    
    sVect<UInt> indicesRow(numRow), indicesCol(numCol);
    sVect<Real> mat;
    
    localMatrix.indexG_row(indicesRow);
    localMatrix.indexG_col(indicesCol);
    mat = localMatrix.matrix();
    
    cout << "Row : " << endl << indicesRow << endl;
    cout << "Col : " << endl << indicesCol << endl;
    cout << "Mat : " << endl << mat << endl;
  }
}
