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
#include "feOscPr3d.hpp"
#include "feOscPr3d_extern.h"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "opEL3dA_CX.hpp"
#include "integratorOEL3d_OSC.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef pMapItem                  PMAPTYPE;
  typedef feOscPr3d<1,1,PMAPTYPE>   FIELD_FETYPE;
  typedef komplex                   FIELD_DOFTYPE;
  typedef feOscPr3d<1,1,PMAPTYPE>   TEST_FETYPE;
  typedef komplex                   TEST_DOFTYPE;
  typedef fePr3d<0,PMAPTYPE>        COEFF_FETYPE;
  typedef Real                      COEFF_DOFTYPE;

  typedef typename FIELD_FETYPE::GEOSHAPE GEOSHAPE;

  //Loading
  string meshFile = "./tests/morganaMeshes/mignon3dC.msh";
  
  meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Typedefs
  typedef feStaticField3d<FIELD_FETYPE,FIELD_DOFTYPE,dms3d_vectMajor,dms3d_allMode> FIELD;
  typedef feStaticField3d<TEST_FETYPE,TEST_DOFTYPE,dms3d_vectMajor,dms3d_allMode>   TEST;
  typedef feStaticField3d<COEFF_FETYPE,COEFF_DOFTYPE,dms3d_vectMajor,dms3d_allMode> COEFF;
  
  typedef typename FIELD::MESH3D        MESH3D;
  typedef typename FIELD::CONNECT3D     CONNECT3D;
  typedef opEL3dA_CX<FIELD,TEST,COEFF>  OPERATOR;
  
  typedef OPERATOR::FIELD_OPTIONS FIELD_OPTIONS;
  typedef OPERATOR::TEST_OPTIONS  TEST_OPTIONS;
  typedef COEFF::OPTIONS          COEFF_OPTIONS;
  
  //Download info
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //Field startup
  COEFF_OPTIONS coeffOptions;
  
  COEFF coeff;
  coeff.setCommunicator(world);
  coeff.setGeometry(grid3d,connectGrid3d);
  coeff.setOptions(coeffOptions);
  coeff.startup();
  
  //Field initialization
  feStaticField3dGlobalManip<COEFF_FETYPE,COEFF_DOFTYPE,dms3d_vectMajor,dms3d_allMode> manipulator(world);
  
  sArray<string> functions(1,1); functions(1,1) = "1.0";
  manipulator.initilize(functions,coeff);
  
  //Build the feCards
  typedef feOscPr3d_extern<1,1,PMAPTYPE>::FECARDS FECARDS;
  
  Real K = 2.5;
  
  feOscPr3d_extern<1,1,PMAPTYPE> feExtern;
  feExtern.setGeometry(grid3d,connectGrid3d);
  feExtern.setCommDev(world);
  FECARDS feCards = feExtern.buildFeCards(K);
  
  //Operator startup
  FIELD_OPTIONS fieldOptions;
  TEST_OPTIONS  testOptions;
  
  OPERATOR operatore(world,*grid3d,*connectGrid3d);
  operatore.setOptions_field(fieldOptions);
  operatore.setOptions_test(testOptions);
  operatore.setFeCards_field(feCards);
  operatore.setFeCards_test(feCards);
  operatore.setCoeff(coeff);
  operatore.startup();
  
  //Eval
  UInt el = 1;
  operatore.setElement(el);
  operatore.setIJ(1,1,1,1);
  
  //Integration 
  integratorOEL3d_OSC<OPERATOR,STANDARD,1> integrator;
  
  UInt pid = world.rank();
  
  if(pid == 0)
  {
    sVect<Real> mat = integrator.integration(operatore);
    
    cout << "element" << endl;
    cout << grid3d->getElementNodesL(el) << endl << endl;
    
    cout << "Matrix" << endl;
    cout << mat << endl;
  }
}
