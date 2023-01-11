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
#include "feStaticField2dGlobalManip.hpp"

#include "operatorLA2d.hpp"
#include "opLA2dA.hpp"
#include "localMatrixLA2d.hpp"
#include "matrixBuilder.hpp"


using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  Teuchos::RCP<communicator> world(new communicator);
  
  UInt pid = world->rank();

  assert(world->size() == 2);
  
  typedef pMapItem             PMAPTYPE;
  typedef fePr2d<1,PMAPTYPE>   FIELD_FETYPE;
  typedef Real                 FIELD_DOFTYPE;
  typedef fePr1d<1,PMAPTYPE>   TEST_FETYPE;
  typedef Real                 TEST_DOFTYPE;

  typedef typename FIELD_FETYPE::GEOSHAPE  GEOSHAPE2D;
  typedef typename TEST_FETYPE::GEOSHAPE   GEOSHAPE1D;

  
  //Loading
  string meshFile = "./tests/morganaMeshes/mignon2dB.msh";
  
  meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  
  //Testing
  typedef feStaticField2d<FIELD_FETYPE,FIELD_DOFTYPE,dms2d_vectMajor,dms2d_allMode>  FIELD;
  typedef feStaticField1d<TEST_FETYPE,TEST_DOFTYPE,dms1d_vectMajor,dms1d_geoIdMode>  TEST;
  
  typedef typename FIELD::MESH2D     MESH2D;
  typedef typename FIELD::CONNECT2D  CONNECT2D;
  
  typedef typename TEST::MESH1D      MESH1D;
  typedef typename TEST::CONNECT1D   CONNECT1D;
  
  typedef opLA2dA<FIELD,TEST>        OPERATOR;
  
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
  
  testOptions.addGeoId(1);
  testOptions.addGeoId(2);
  
  Teuchos::RCP<OPERATOR> operatore(new OPERATOR(*world,*grid2d,*connectGrid2d,*grid1d,*connectGrid1d));
  
  operatore->setOptions_field(fieldOptions);
  operatore->setOptions_test(testOptions);
  operatore->startup();
  
  
  //Local Matrix
  typedef localMatrixLA2d<FIELD,TEST,intDefault,STANDARD,2>  LOCALMATRIX;
  typedef matrixBuilder<LOCALMATRIX>                         ASSEMBLER;
  typedef typename ASSEMBLER::EPETRA_CRS EPETRA_CRS;
  
  RCP<EPETRA_CRS> matrix;
  
  ASSEMBLER assembler(world,operatore);
  assembler.buildEpetraCrs(matrix);
 
  cout << *matrix << endl;
}
