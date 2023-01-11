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

#include "fePr3d.hpp"
#include "feRt0LT3d.hpp"
#include "feRt0LT3d_extern.hpp"

#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "opEL3dE.hpp"
#include "localMatrixEL3d.hpp"
#include "matrixBuilder.hpp"


// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  Teuchos::RCP<communicator> world(new communicator);
  
  typedef pMapItem             PMAPTYPE;
  typedef feRt0LT3d<PMAPTYPE>  FIELD_FETYPE;
  typedef Real                 FIELD_DOFTYPE;
  typedef feRt0LT3d<PMAPTYPE>  TEST_FETYPE;
  typedef Real                 TEST_DOFTYPE;
  typedef fePr3d<1,PMAPTYPE>   COEFF_FETYPE;
  typedef Real                 COEFF_DOFTYPE;

  typedef typename FIELD_FETYPE::GEOSHAPE  GEOSHAPE;

  
  //Loading----------------------------------------------------------
  string meshFile = "./tests/morganaMeshes/mignon3dB.msh";
  
  meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  
  //Typedefs---------------------------------------------------------
  typedef feStaticField3d<FIELD_FETYPE, FIELD_DOFTYPE, dms3d_vectMajor, dms3d_allMode>  FIELD;
  typedef feStaticField3d<TEST_FETYPE,  TEST_DOFTYPE,  dms3d_vectMajor, dms3d_allMode>  TEST;
  typedef feStaticField3d<COEFF_FETYPE, COEFF_DOFTYPE, dms3d_vectMajor, dms3d_allMode>  COEFF;
  
  typedef typename FIELD::MESH3D     MESH3D;
  typedef typename FIELD::CONNECT3D  CONNECT3D;
  typedef opEL3dE<FIELD,TEST,COEFF>  OPERATOR;
  
  typedef OPERATOR::FIELD_OPTIONS FIELD_OPTIONS;
  typedef OPERATOR::TEST_OPTIONS  TEST_OPTIONS;
  typedef COEFF::OPTIONS          COEFF_OPTIONS;
  
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  
  //Coeff startup----------------------------------------------------
  COEFF_OPTIONS coeffOptions;
  
  COEFF coeff;
  coeff.setCommunicator(world);
  coeff.setGeometry(grid3d,connectGrid3d);
  coeff.setOptions(coeffOptions);
  coeff.startup();
  
  
  //Coeff initialization---------------------------------------------
  sArray<string> functions(1,1);
  functions(1,1) = "0.5";
  
  feStaticField3dGlobalManip<COEFF_FETYPE,COEFF_DOFTYPE,dms3d_vectMajor,dms3d_allMode> manipulator(world);
  manipulator.initilize(functions,coeff);
  
  
  //Testing---------------------------------------------------------
  typedef feRt0LT3d_extern<PMAPTYPE>  EXTERNAL;
  typedef typename EXTERNAL::FECARDS  FECARDS;
  
  EXTERNAL external(grid3d,connectGrid3d);
  external.setCommDev(world);
  
  FECARDS feCards = external.buildFeCards();
  
  
  //Operator---------------------------------------------------------
  FIELD_OPTIONS fieldOptions;
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<OPERATOR> operatore(new OPERATOR(*world,*grid3d,*connectGrid3d));
  operatore->setFeCards_field(feCards);
  operatore->setFeCards_test(feCards);
  operatore->setOptions_field(fieldOptions);
  operatore->setOptions_test(testOptions);
  operatore->setCoeff(coeff);
  operatore->startup();
  
  
  //Local Matrix-----------------------------------------------------
  typedef localMatrixEL3d<FIELD,TEST,intDefault,STANDARD,2>  LOCALMATRIX;
  typedef matrixBuilder<LOCALMATRIX>                         ASSEMBLER;
  typedef typename ASSEMBLER::EPETRA_CRS EPETRA_CRS;
  
  RCP<EPETRA_CRS> matrix;
  
  ASSEMBLER assembler(world,operatore);
  assembler.buildEpetraCrs(matrix);
 
  cout << *matrix << endl;
  
  if(world->rank() == 0)
  { cout << grid3d->getFaces() << endl; }
}
