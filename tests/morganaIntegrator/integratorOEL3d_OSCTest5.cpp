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
#include "feOscCr3d.hpp"
#include "feOscCr3d_extern.h"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "opEL3dB_CX.hpp"
#include "integratorOEL3d_OSC.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  typedef pMapItem                  PMAPTYPE;
  typedef feOscCr3d<1,PMAPTYPE>     FIELD_FETYPE;
  typedef komplex                   FIELD_DOFTYPE;
  typedef feOscCr3d<1,PMAPTYPE>     TEST_FETYPE;
  typedef komplex                   TEST_DOFTYPE;

  typedef typename FIELD_FETYPE::GEOSHAPE GEOSHAPE;

  //Loading
  string meshFile = "./geometries/cubes3d/testCubeC.msh";
  
  meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  //Typedefs
  typedef feStaticField3d<FIELD_FETYPE,FIELD_DOFTYPE,dms3d_vectMajor,dms3d_allMode> FIELD;
  typedef feStaticField3d<TEST_FETYPE,TEST_DOFTYPE,dms3d_vectMajor,dms3d_allMode>   TEST;
  
  typedef typename FIELD::MESH3D    MESH3D;
  typedef typename FIELD::CONNECT3D CONNECT3D;
  typedef opEL3dB_CX<FIELD,TEST>    OPERATOR;
  
  typedef OPERATOR::FIELD_OPTIONS FIELD_OPTIONS;
  typedef OPERATOR::TEST_OPTIONS  TEST_OPTIONS;
  
  //Download info
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  //Build the feCards
  typedef feOscCr3d_extern<1,PMAPTYPE>::FECARDS FECARDS;
  
  Real K = 10.0;
  
  feOscCr3d_extern<1,PMAPTYPE> feExtern;
  feExtern.setGeometry(grid3d,connectGrid3d);
  feExtern.setCommDev(world);
  FECARDS feCardsField = feExtern.buildFeCards(K);
  
  FECARDS feCardsTest(feCardsField);
  
  for(UInt i=1; i <= feCardsTest.size(); ++i)
  { feCardsTest.getL(i).conjugate(); }
  
  //Operator startup
  FIELD_OPTIONS fieldOptions;
  TEST_OPTIONS  testOptions;
  
  OPERATOR operatore(world,*grid3d,*connectGrid3d);
  operatore.setOptions_field(fieldOptions);
  operatore.setOptions_test(testOptions);
  operatore.setFeCards_field(feCardsField);
  operatore.setFeCards_test(feCardsTest);
  operatore.startup();
  
  integratorOEL3d_OSC<OPERATOR,STANDARD,4> integrator1;
  integratorOEL3d_OSC<OPERATOR,STANDARD,5> integrator3;
  
  Real diffLoc = 0.0, diff = 0.0;
  sVect<Real> mat1, mat3;
  UInt elRef = 0.0;
  
  operatore.setIJ(1,1,1,1);
  
  /*for(UInt el=1; el <= grid3d->getNumElements(); el++)
  {
    operatore.setElement(el);
    mat1 = integrator1.integration(operatore);
    mat3 = integrator3.integration(operatore);
    
    diffLoc = 0.0;
    
    for(UInt i=1; i <= mat1.size(); ++i)
    { diffLoc = std::max(diffLoc, mat1(i) - mat3(i)); }
    
    if(diffLoc >= diff)
    {
      elRef = el;
      diff  = std::max(diff,diffLoc);
    }
  }
  
  cout << "elRef : " << elRef << endl;
  cout << "diff  : " << diff  << endl;*/
  
  elRef = 417;
  
  operatore.setElement(elRef);
  cout << "Integrator 1" << endl;
  mat1 = integrator1.integration(operatore);
  cout << "Integrator 3" << endl;
  mat3 = integrator3.integration(operatore);
  
  cout << "mat1: " << mat1(7) << endl;
  cout << "mat3: " << mat3(7) << endl;
}
