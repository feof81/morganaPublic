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

#include "feDynamicField3d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "localVectorEL3d.hpp"
#include "fnEL3dA.hpp"

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
  typedef fePr3d<1,PMAPTYPE>    TEST_FETYPE;
  typedef Real                  TEST_DOFTYPE;
  typedef fePr3d<0,PMAPTYPE>    COEFF_FETYPE;
  typedef Real                  COEFF_DOFTYPE;
  
  typedef typename TEST_FETYPE::FECARD  TEST_FECARD;
  typedef typename TEST_FETYPE::GEOSHAPE GEOSHAPE;

  
  //Loading
  string meshFile = "./tests/morganaMeshes/mignon3dB.msh";
  
  meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  
  //Typedefs
  typedef feStaticField3d<TEST_FETYPE,TEST_DOFTYPE,dms3d_vectMajor,dms3d_allMode>    TEST;
  typedef feStaticField3d<COEFF_FETYPE,COEFF_DOFTYPE,dms3d_vectMajor,dms3d_allMode>  COEFF;
  
  typedef typename TEST::MESH3D     MESH3D;
  typedef typename TEST::CONNECT3D  CONNECT3D;
  typedef fnEL3dA<TEST,COEFF>       FUNCTIONAL;
  
  typedef FUNCTIONAL::TEST_OPTIONS  TEST_OPTIONS;
  typedef COEFF::OPTIONS           COEFF_OPTIONS;
  
  //Download info
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //Coeff startup
  COEFF_OPTIONS coeffOptions;
  
  COEFF coeff;
  coeff.setCommunicator(world);
  coeff.setGeometry(grid3d,connectGrid3d);
  coeff.setOptions(coeffOptions);
  coeff.startup();
  
  //Coeff initialization
  feStaticField3dGlobalManip<COEFF_FETYPE,COEFF_DOFTYPE,dms3d_vectMajor,dms3d_allMode> manipulator(world);
  
  sArray<string> functions(1,1);
  functions(1,1) = "2";
  
  manipulator.initilize(functions,coeff);
  
  //Operator startup
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<FUNCTIONAL> operatore(new FUNCTIONAL(world,*grid3d,*connectGrid3d));
  operatore->setOptions_test(testOptions);  
  operatore->setCoeff(coeff);
  operatore->startup();
  
  //Local matrix test
  UInt el = 1;
  
  operatore->setElement(el);
  operatore->setIJ(1,1);
  
  
  typedef localVectorEL3d<TEST> LOCALVECTOR;
  LOCALVECTOR locVect(operatore);
  
  
  if(world.rank() == 0)
  {    
    UInt numTest = locVect.numIndex_row();
    
    sVect<UInt> indices(numTest);
    
    locVect.indexG_row(indices);
    sVect<Real> mat = locVect.vector();
    
    cout << "Indices : " << endl << indices << endl;
    cout << "Matrix eval :" << endl << mat << endl;
  }
}
