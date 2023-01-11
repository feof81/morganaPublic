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

#include "fePr2d.hpp"

#include "feDynamicField2d.hpp"
#include "feStaticField2d.hpp"
#include "feStaticField2dGlobalManip.hpp"

#include "fnEL2dA.hpp"

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
  typedef fePr2d<1,PMAPTYPE>    TEST_FETYPE;
  typedef Real                  TEST_DOFTYPE;
  typedef fePr2d<1,PMAPTYPE>    COEFF_FETYPE;
  typedef Real                  COEFF_DOFTYPE;
  
  typedef typename TEST_FETYPE::FECARD  TEST_FECARD;
  typedef typename TEST_FETYPE::GEOSHAPE GEOSHAPE;

  
  //Loading
  string meshFile = "./tests/morganaMeshes/mignon2dB.msh";
  
  meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  
  //Typedefs
  typedef feStaticField2d<TEST_FETYPE,TEST_DOFTYPE,dms2d_vectMajor,dms2d_allMode>    TEST;
  typedef feStaticField2d<COEFF_FETYPE,COEFF_DOFTYPE,dms2d_vectMajor,dms2d_allMode>  COEFF;
  
  typedef typename TEST::MESH2D     MESH2D;
  typedef typename TEST::CONNECT2D  CONNECT2D;
  typedef fnEL2dA<TEST,COEFF>       FUNCTIONAL;
  
  typedef FUNCTIONAL::TEST_OPTIONS  TEST_OPTIONS;
  typedef COEFF::OPTIONS            COEFF_OPTIONS;
  
  //Download info
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
  //Coeff startup
  COEFF_OPTIONS coeffOptions;
  
  COEFF coeff;
  coeff.setCommunicator(world);
  coeff.setGeometry(grid2d,connectGrid2d);
  coeff.setOptions(coeffOptions);
  coeff.startup();
  
  //Coeff initialization
  feStaticField2dGlobalManip<COEFF_FETYPE,COEFF_DOFTYPE,dms2d_vectMajor,dms2d_allMode> manipulator(world);
  
  sArray<string> functions(1,1);
  functions(1,1) = "y";
  
  manipulator.initilize(functions,coeff);
  
  //Operator startup
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<FUNCTIONAL> operatore(new FUNCTIONAL(world,*grid2d,*connectGrid2d));
  operatore->setOptions_test(testOptions);  
  operatore->setCoeff(coeff);
  operatore->startup();
  
  //Eval
  UInt el = 1;
  
  operatore->setElement(el);
  operatore->setIJ(1,1);
  
  if(world.rank() == 1)
  {
    typedef FUNCTIONAL::TEST_OUTTYPE TEST_OUTTYPE;
    
    point3d Y(0.0, 0.5, 0.0);
    UInt numTest = operatore->numIndex_test();
    
    sVect<UInt>         indices(numTest);
    sVect<TEST_OUTTYPE> valTest(numTest);
    sVect<Real>         mat(numTest);
    
    operatore->indexG_test(indices);
    operatore->eval_test(Y,valTest);
    operatore->eval(Y,mat);
    
    cout << "Indices : " << endl << indices << endl;
    cout << "Test eval : " << endl << valTest << endl;
    cout << "Matrix eval :" << endl << mat << endl;
  }
}
