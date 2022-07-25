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
#include "feOscPr2d.hpp"
#include "feOscPr2d_extern.h"
#include "feOscPr3d.hpp"

#include "dofMapStatic2d.hpp"
#include "elCardFeeder2d.hpp"
#include "feStaticField2d.hpp"
#include "feStaticField2dGlobalManip.hpp"

#include "fnLA3dA_CX.hpp"
#include "traitsIntegratorFLA3d_OSC.hpp"

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
  typedef feOscPr2d<0,0,PMAPTYPE>   TEST_FETYPE;
  typedef komplex                   TEST_DOFTYPE;
  typedef fePr2d<0,PMAPTYPE>        COEFF_FETYPE;
  typedef komplex                   COEFF_DOFTYPE;

  typedef linearTetra GEOSHAPE3D;

  //Loading
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> INIT;
  
  string meshFile = "./tests/morganaMeshes/mignon3dC.msh";
  
  INIT init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Typedefs
  typedef feStaticField2d<TEST_FETYPE,  TEST_DOFTYPE,  dms2d_vectMajor, dms2d_allMode>  TEST;
  typedef feStaticField2d<COEFF_FETYPE, COEFF_DOFTYPE, dms2d_vectMajor, dms2d_allMode>  COEFF;
  
  typedef typename INIT::MESH2D       MESH2D;
  typedef typename INIT::MESH3D       MESH3D;
  typedef typename INIT::CONNECT2D CONNECT2D;
  typedef typename INIT::CONNECT3D CONNECT3D;
  typedef fnLA3dA_CX<TEST,GEOSHAPE3D,COEFF>  FUNCTIONAL;
  
  typedef FUNCTIONAL::TEST_OPTIONS TEST_OPTIONS;
  typedef COEFF::OPTIONS          COEFF_OPTIONS;
  
  //Download info
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //Field startup
  COEFF_OPTIONS coeffOptions;
  
  COEFF coeff;
  coeff.setCommunicator(world);
  coeff.setGeometry(grid2d,connectGrid2d);
  coeff.setOptions(coeffOptions);
  coeff.startup();
  
  //Field initialization
  feStaticField2dGlobalManip<COEFF_FETYPE,COEFF_DOFTYPE,dms2d_vectMajor,dms2d_allMode> manipulator(world);
  
  sArray<string> functions(2,1);
  functions(1,1) = "1.0";
  functions(2,1) = "1.0";
  manipulator.initilize(functions,coeff);
  
  //Build the feCards
  typedef feOscPr2d_extern<0,0,PMAPTYPE>::FECARDS FECARDS;
  
  Real K = 2.5;
  
  feOscPr2d_extern<0,0,PMAPTYPE> feExtern;
  feExtern.setGeometry(grid2d,connectGrid2d);
  feExtern.setCommDev(world);
  FECARDS feCards = feExtern.buildFeCards(K);
  
  //Operator startup
  TEST_OPTIONS  testOptions;
  
  FUNCTIONAL funzionale(world,*grid3d,*connectGrid3d,*grid2d,*connectGrid2d);
  funzionale.setOptions_test(testOptions);
  funzionale.setFeCards_test(feCards);
  funzionale.setCoeff(coeff);
  funzionale.startup();
  
  //Eval
  UInt el = 1;
  UInt fc = 2;
  funzionale.setElement(el);
  funzionale.setLocFace(fc);
  funzionale.setIJ(1,1);
  
  //Integration 
  traitsIntegratorFLA3d_OSC<FUNCTIONAL,linearTriangle,STANDARD,1> integrator;
  
  UInt pid = world.rank();
  
  if(pid == 0)
  {
    cout << "isBoundary : " << funzionale.isBoundary() << endl;
    
    cout << "element" << endl;
    cout << grid2d->getElementNodesL(el) << endl;
    
    sVect<Real> vect = integrator.integration(funzionale);
    
    cout << "Matrix" << endl;
    cout << vect << endl;
  }
}

