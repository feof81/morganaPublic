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
#include "feRt0LT2d.hpp"
#include "feRt0LT2d_extern.hpp"

#include "dofMapStatic2d.hpp"
#include "elCardFeeder2d.hpp"

#include "feStaticField1d.hpp"
#include "feStaticField2d.hpp"
#include "feStaticField1dGlobalManip.hpp"
#include "feStaticField2dGlobalManip.hpp"

#include "functionalBF2d.hpp"
#include "fnBF2dA.hpp"
#include "localVectorBF2d.hpp"
#include "vectorBuilder.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  Teuchos::RCP<communicator> world(new communicator);

  assert(world->size() == 2);
  
  typedef pMapItemShare        PMAPTYPE;
  typedef feRt0LT2d<PMAPTYPE>  TEST_FETYPE;
  typedef Real                 FIELD_DOFTYPE;
  typedef fePr1d<1,PMAPTYPE>   COEFF_FETYPE;
  typedef Real                 COEFF_DOFTYPE;
  
  typedef typename TEST_FETYPE::FECARD   TEST_FECARD;
  typedef typename TEST_FETYPE::GEOSHAPE GEOSHAPE;
  
  typedef meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> INIT;

  
  //Loading
  string meshFile = "./tests/morganaMeshes/mignon2dB.msh";
  
  INIT init(world);
  init.gmMesh_to_stdB(meshFile);
  
  
  //Typedefs
  typedef feStaticField2d<TEST_FETYPE,FIELD_DOFTYPE,dms2d_vectMajor,dms2d_allMode>  TEST;
  
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
  typedef feStaticField1d<COEFF_FETYPE,COEFF_DOFTYPE,dms1d_vectMajor,dms1d_allMode> COEFF;
  
  dofMapStatic1d_options mapOptions;
  
  COEFF field1d;
  field1d.setCommunicator(world);
  field1d.setGeometry(grid1d,connectGrid1d);
  field1d.setOptions(mapOptions);
  field1d.startup();
  
  //Field initialization
  feStaticField1dGlobalManip<COEFF_FETYPE,COEFF_DOFTYPE,dms1d_vectMajor,dms1d_allMode> manipulator(world);
  
  sArray<string> functions(1,1);
  functions(1,1) = "2";
  
  manipulator.initilize(functions,field1d);
  
  
  //Extern FE
  typedef feRt0LT2d_extern<PMAPTYPE>  EXTERN;
  typedef typename EXTERN::FECARDS    FECARDS;
  
  EXTERN feeder(grid2d,connectGrid2d);
  feeder.setCommDev(world);
  FECARDS feCards = feeder.buildFeCards();
  
  
  //Operator startup
  typedef fnBF2dA<TEST,COEFF>       FUNCTIONAL;
  typedef FUNCTIONAL::TEST_OPTIONS  TEST_OPTIONS;
  typedef FUNCTIONAL::TEST_OUTTYPE  TEST_OUTTYPE;
  
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<FUNCTIONAL> funzionale(new FUNCTIONAL(*world,*grid2d,*connectGrid2d,*grid1d,*connectGrid1d));
  funzionale->setOptions_test(testOptions);
  funzionale->setFeCards_test(feCards);
  funzionale->setCoeff(field1d);
  funzionale->startup();
  
  
  //Vector
  typedef localVectorBF2d<TEST,intDefault,STANDARD,2> LOCVECTOR;
  typedef vectorBuilder<LOCVECTOR>                    FEBUILDER;
  typedef typename FEBUILDER::EPETRA_VECTOR           EPETRA_VECTOR;
   
  RCP<EPETRA_VECTOR> vector;
  
  FEBUILDER feVector;
  feVector.setFunctional(funzionale);
  feVector.setCommDev(world);
  
  feVector.buildEpetraVector(vector);
  cout << *vector << endl;
}
