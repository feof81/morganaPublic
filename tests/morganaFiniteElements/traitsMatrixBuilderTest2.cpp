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
#include "feRt0LT3d.hpp"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"
#include "feStaticField2d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "operatorLA3d.hpp"
#include "opLA3dA.hpp"

#include "localMatrixLA3d.hpp"
#include "traitsMatrixBuilder.hpp"


using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  Teuchos::RCP<communicator> world(new communicator);
  
  UInt pid = world->rank();
  
  typedef pMapItem             PMAPTYPE;
  typedef fePr3d<1,PMAPTYPE>   FIELD_FETYPE;
  typedef Real                 FIELD_DOFTYPE;
  typedef fePr2d<0,PMAPTYPE>   TEST_FETYPE;
  typedef Real                 TEST_DOFTYPE;

  typedef typename FIELD_FETYPE::GEOSHAPE  GEOSHAPE3D;
  typedef typename TEST_FETYPE::GEOSHAPE   GEOSHAPE2D;

  
  //Loading
  string meshFile = "./tests/morganaMeshes/mignon3dC.msh";
  
  meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile);
  
  
  //Testing
  typedef feStaticField3d<FIELD_FETYPE,FIELD_DOFTYPE,dms3d_vectMajor,dms3d_allMode>  FIELD;
  typedef feStaticField2d<TEST_FETYPE,TEST_DOFTYPE,dms2d_vectMajor,dms2d_geoIdMode>  TEST;
  
  typedef typename FIELD::MESH3D     MESH3D;
  typedef typename FIELD::CONNECT3D  CONNECT3D;
  
  typedef typename TEST::MESH2D      MESH2D;
  typedef typename TEST::CONNECT2D   CONNECT2D;
  
  typedef opLA3dA<FIELD,TEST>        OPERATOR;
  
  typedef OPERATOR::FIELD_OPTIONS    FIELD_OPTIONS;
  typedef OPERATOR::TEST_OPTIONS     TEST_OPTIONS;
  
  typedef localMatrixLA3d<FIELD,TEST,intDefault,STANDARD,1> LOCALMATRIX;
  typedef traitsMatrixBuilder<LOCALMATRIX,pMapItem>         ASSEMBLER;
  
  
  //Download data
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //Operator
  FIELD_OPTIONS fieldOptions;
  TEST_OPTIONS  testOptions;
  
  testOptions.addGeoId(2);
  
  Teuchos::RCP<OPERATOR> operatore(new OPERATOR(*world,*grid3d,*connectGrid3d,*grid2d,*connectGrid2d));
  
  operatore->setOptions_field(fieldOptions);
  operatore->setOptions_test(testOptions);
  operatore->startup();
  
  
  //Global assembler
  typedef typename ASSEMBLER::TPETRA_CRS TPETRA_CRS;
  
  RCP<TPETRA_CRS> matrix;
  
  ASSEMBLER assembler(world,operatore);
  assembler.buildTpetraCrs(matrix);
  
  
  //Printout
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  matrix->describe(*fos, Teuchos::VERB_EXTREME);
  cout << matrix->description() << endl;
}
