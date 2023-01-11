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

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"
#include "feStaticField2d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "opBB3dA.hpp"
#include "localMatrixBB3d.hpp"
#include "intPolicyFeeder3d_STC.hpp"


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
  typedef fePr3d<1,PMAPTYPE>   FIELD_FETYPE;
  typedef Real                 FIELD_DOFTYPE;
  typedef fePr3d<1,PMAPTYPE>   TEST_FETYPE;
  typedef Real                 TEST_DOFTYPE;

  typedef typename FIELD_FETYPE::GEOSHAPE  GEOSHAPE3D;
  typedef typename TEST_FETYPE::GEOSHAPE   GEOSHAPE2D;

  
  //Loading
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> INIT;
  
  string meshFile = "./tests/morganaMeshes/mignon3dD.msh";
  
  INIT init(world);
  init.gmMesh_to_stdA(meshFile, false);
  
  
  //Testing
  typedef feStaticField3d<FIELD_FETYPE,FIELD_DOFTYPE,dms3d_vectMajor,dms3d_allMode>  FIELD;
  typedef feStaticField3d<TEST_FETYPE,TEST_DOFTYPE,dms3d_vectMajor,dms3d_allMode>    TEST;
  
  typedef typename INIT::MESH3D     MESH3D;
  typedef typename INIT::CONNECT3D  CONNECT3D;
  
  typedef typename INIT::MESH2D      MESH2D;
  typedef typename INIT::CONNECT2D   CONNECT2D;
  
  typedef opBB3dA<FIELD,TEST>        OPERATOR;
  
  typedef OPERATOR::FIELD_OPTIONS    FIELD_OPTIONS;
  typedef OPERATOR::TEST_OPTIONS     TEST_OPTIONS;
  
  
  //Download data
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //Operator
  FIELD_OPTIONS fieldOptions;
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<OPERATOR> operatore(new OPERATOR(world,*grid3d,*connectGrid3d,*grid2d,*connectGrid2d));
  
  operatore->setOptions_field(fieldOptions);
  operatore->setOptions_test(testOptions);
  operatore->startup();
  
  
  //Int cards--------------------------------------------------------
  sVect<UInt> geoIds3d;
  geoIds3d.push_back(1);
  geoIds3d.push_back(2);
  
  typedef intPolicyFeeder3d_STC<GEOSHAPE3D,PMAPTYPE> POLICYFEEDER;
  typedef POLICYFEEDER::INTCARDS                     INTCARDS;
  
  POLICYFEEDER policyFeeder(grid3d);
  INTCARDS intCards = policyFeeder.getIntCards(geoIds3d);
  
  
  //Local matrix
  typedef localMatrixBB3d<FIELD,TEST,OBB3d_STC,STANDARD,2>  LOCALMATRIX;
  
  if(pid == 1)
  {
    //cout << "Grid3d" << endl;
    //cout << grid3d->getElements() << endl;
    
    //Geo data
    UInt el = 1;
    UInt I_test=1, J_test=1, I_field=1, J_field=1;
    
    //Active ids2d
    set<UInt> geoIds2d;
    geoIds2d.insert(2);
    
    operatore->setElement(el);
    operatore->setIJ(I_test, J_test, I_field, J_field);
    operatore->setGeoIds(geoIds2d);
    
    //Local matrix
    LOCALMATRIX locMatrix(operatore);
    locMatrix.setIntCards(intCards);
    
    sVect<UInt> indRow(locMatrix.numIndex_row());
    sVect<UInt> indCol(locMatrix.numIndex_col());
    sVect<Real> mat = locMatrix.matrix();
    
    locMatrix.indexG_col(indCol);
    locMatrix.indexG_row(indRow);
    
    cout << "numRows : " << locMatrix.numIndex_row() << endl;
    cout << "numCols : " << locMatrix.numIndex_col() << endl;
    
    cout << "rows   : " << endl << indRow << endl;
    cout << "cols   : " << endl << indCol << endl;
    cout << "matrix : " << endl << mat << endl;
  }
}

