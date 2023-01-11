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
#include "feRt0LT3d.hpp"
#include "feRt0LT3d_extern.hpp"

#include "feStaticField2d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField2dGlobalManip.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "fnBF3dC.hpp"
#include "localVectorBF3d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef pMapItemShare        PMAPTYPE;
  typedef feRt0LT3d<PMAPTYPE>  TEST_FETYPE;
  typedef Real                 FIELD_DOFTYPE;
  typedef fePr2d<1,PMAPTYPE>   COEFF_FETYPE;
  typedef Real                 COEFF_DOFTYPE;
  
  typedef typename TEST_FETYPE::FECARD   TEST_FECARD;
  typedef typename TEST_FETYPE::GEOSHAPE GEOSHAPE;
  
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> INIT;

  
  //Loading
  string meshFile = "./tests/morganaMeshes/mignon3dB.msh";
  
  INIT init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  
  //Typedefs
  typedef feStaticField3d<TEST_FETYPE,FIELD_DOFTYPE,dms3d_vectMajor,dms3d_allMode>  TEST;
  
  typedef typename INIT::MESH2D     MESH2D;
  typedef typename INIT::MESH3D     MESH3D;
  typedef typename INIT::CONNECT2D  CONNECT2D;
  typedef typename INIT::CONNECT3D  CONNECT3D;  
  
  //Download info
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();  
  
  //Extern FE
  typedef feRt0LT3d_extern<PMAPTYPE>  EXTERN;
  typedef typename EXTERN::FECARDS    FECARDS;
  
  EXTERN feeder(grid3d,connectGrid3d);
  feeder.setCommDev(world);
  FECARDS feCards = feeder.buildFeCards();
  
  
  //Operator startup 
  typedef fnBF3dC<TEST>             FUNCTIONAL;
  typedef FUNCTIONAL::TEST_OPTIONS  TEST_OPTIONS;
  typedef FUNCTIONAL::TEST_OUTTYPE  TEST_OUTTYPE;
  
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<FUNCTIONAL> funzionale(new FUNCTIONAL(world,*grid3d,*connectGrid3d,*grid2d,*connectGrid2d));
  funzionale->setOptions_test(testOptions);
  funzionale->setFeCards_test(feCards);
  funzionale->startup();
  
  
  //Integrator
  typedef localVectorBF3d<TEST,intDefault,STANDARD,2> LOCALVECTOR;
  LOCALVECTOR locvector(funzionale);
  
  if(world.rank() == 1)
  {
    //Eval
    UInt el = 2;
    
    funzionale->setElement(el);
    funzionale->setIJ(1,1);
    
    UInt numRow = locvector.numIndex_row();
    sVect<UInt> indices(numRow);
    
    locvector.indexG_row(indices);
    sVect<Real> vect = locvector.vector();
    
    cout << "Indices" << endl << indices << endl;
    cout << "Vect" << endl << vect << endl;
  }
}
