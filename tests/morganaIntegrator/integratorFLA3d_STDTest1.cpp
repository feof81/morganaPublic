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
#include "feStaticField2dGlobalManip.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "operatorLA3d.hpp"
#include "fnLA3dA.hpp"
#include "integratorFLA3d_STD.hpp"
#include "morganaIntegrator.hpp"


using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef pMapItemShare        PMAPTYPE;
  typedef fePr2d<1,PMAPTYPE>   TEST_FETYPE;
  typedef Real                 TEST_DOFTYPE;
  typedef fePr2d<1,PMAPTYPE>   COEFF_FETYPE;
  typedef point3d              COEFF_DOFTYPE;

  typedef          linearTetra            GEOSHAPE3D;
  typedef typename TEST_FETYPE::GEOSHAPE  GEOSHAPE2D;

  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> INIT;

  
  //Loading
  string meshFile = "../tests/morganaMeshes/mignon3dB.msh";
  
  INIT init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  
  //Typedefs
  typedef feStaticField2d<TEST_FETYPE,TEST_DOFTYPE,dms2d_vectMajor,dms2d_allMode>  TEST;
  
  typedef typename INIT::MESH2D     MESH2D;
  typedef typename INIT::MESH3D     MESH3D;
  typedef typename INIT::CONNECT2D  CONNECT2D;
  typedef typename INIT::CONNECT3D  CONNECT3D;  
  
  //Download info
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  
 
  //Set the coeff field
  typedef feStaticField2d<COEFF_FETYPE,COEFF_DOFTYPE,dms2d_vectMajor,dms2d_geoIdMode> COEFF;
  
  dofMapStatic2d_options mapOptions;
  mapOptions.addGeoId(1);
  mapOptions.addGeoId(2);
  
  COEFF field2d;
  field2d.setCommunicator(world);
  field2d.setGeometry(grid2d,connectGrid2d);
  field2d.setOptions(mapOptions);
  field2d.startup();
  
  //Field initialization
  feStaticField2dGlobalManip<COEFF_FETYPE,COEFF_DOFTYPE,dms2d_vectMajor,dms2d_geoIdMode> manipulator(world);
  
  sArray<string> functions(3,1);
  functions(1,1) = "2";
  functions(2,1) = "0";
  functions(3,1) = "0";
  
  manipulator.initilize(functions,field2d);
  
  
  //Operator startup
  typedef fnLA3dA<TEST,GEOSHAPE3D,COEFF>  FUNCTIONAL;
  typedef FUNCTIONAL::TEST_OPTIONS        TEST_OPTIONS;
  typedef FUNCTIONAL::TEST_OUTTYPE        TEST_OUTTYPE;
  
  TEST_OPTIONS  testOptions;
  
  FUNCTIONAL funzionale(world,*grid3d,*connectGrid3d,*grid2d,*connectGrid2d);
  funzionale.setOptions_test(testOptions);
  funzionale.setCoeff(field2d);
  funzionale.startup();
  
  //Integrator
  typedef integratorFLA3d_STD<FUNCTIONAL,STANDARD,1> INTEGRATOR;
  INTEGRATOR integrator;
  
   
  if(world.rank() == 1)
  {
    //Eval
    UInt el = 1;
    UInt fc = 3;
    point3d Ys(0.0, 1.0, 0.0);
  
    funzionale.setElement(el);
    funzionale.setIJ(1,1);
    funzionale.setLocFace(fc);
    
    if(funzionale.isBoundary())
    {
      cout << "IsBoundary" << endl;
      
      //Alloc
      UInt numTest = funzionale.numIndex_test();
      sVect<UInt> indicesTest(numTest);
      
      funzionale.indexG_test(indicesTest);
      sVect<Real> vect = integrator.integration(funzionale);
      
      cout << "Indices : " << endl << indicesTest << endl;
      cout << "Eval : " << endl << vect << endl;
    }
    else
    {
      cout << "NotBoundary" << endl;
    }
  }
}
