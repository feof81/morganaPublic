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
#include "feRt0LT2d.hpp"
#include "feRt0LT2d_extern.hpp"

#include "dofMapStatic2d.hpp"
#include "elCardFeeder2d.hpp"
#include "feStaticField2d.hpp"
#include "feStaticField2dGlobalManip.hpp"

#include "operatorEL2d.hpp"
#include "opEL2dA.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef pMapItem             PMAPTYPE;
  typedef fePr2d<1,PMAPTYPE>   FIELD_FETYPE;
  typedef Real                 FIELD_DOFTYPE;
  typedef fePr2d<1,PMAPTYPE>   TEST_FETYPE;
  typedef Real                 TEST_DOFTYPE;
  typedef fePr2d<0,PMAPTYPE>   COEFF_FETYPE;
  typedef Real                 COEFF_DOFTYPE;

  typedef typename FIELD_FETYPE::GEOSHAPE   GEOSHAPE;

  
  //Loading
  string meshFile = "./tests/morganaMeshes/mignon2dB.msh";
  
  meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  
  //Typedefs
  typedef feStaticField2d<FIELD_FETYPE,FIELD_DOFTYPE,dms2d_vectMajor,dms2d_allMode> FIELD;
  typedef feStaticField2d<TEST_FETYPE,TEST_DOFTYPE,dms2d_vectMajor,dms2d_allMode>   TEST;
  typedef feStaticField2d<COEFF_FETYPE,COEFF_DOFTYPE,dms2d_vectMajor,dms2d_allMode> COEFF;
  
  typedef typename FIELD::MESH2D     MESH2D;
  typedef typename FIELD::CONNECT2D  CONNECT2D;
  typedef opEL2dA<FIELD,TEST,COEFF>  OPERATOR;
  
  typedef OPERATOR::FIELD_OPTIONS FIELD_OPTIONS;
  typedef OPERATOR::TEST_OPTIONS  TEST_OPTIONS;
  typedef COEFF::OPTIONS          COEFF_OPTIONS;
  
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
  functions(1,1) = "x";
  
  manipulator.initilize(functions,coeff);
  
  
  //Operator
  FIELD_OPTIONS fieldOptions;
  TEST_OPTIONS  testOptions;
  
  Teuchos::RCP<OPERATOR> operatore(new OPERATOR(world,*grid2d,*connectGrid2d));
  operatore->setOptions_field(fieldOptions);
  operatore->setOptions_test(testOptions);
  operatore->setCoeff(coeff);
  operatore->startup();
  
  UInt pid = world.rank();
  

  if(pid == 0)
  {
    cout << "Grid 2d - nodes" << endl;
    cout << grid2d->getNodes() << endl << endl;
    
    cout << "Grid 2d - elements" << endl;
    cout << grid2d->getElements() << endl << endl;
    
    sVect<UInt> indices;
    
    //Indices map
    cout << "field dofsL " << operatore->getNumDofsL_field() << endl;
    cout << "field dofsG " << operatore->getNumDofsG_field() << endl;
    cout << "field listL " << operatore->getSizeListL_field() << endl;
    cout << "field listG " << operatore->getSizeListG_field() << endl << endl;
  
    cout << "test dofsL " << operatore->getNumDofsL_test() << endl;
    cout << "test dofsG " << operatore->getNumDofsG_test() << endl;
    cout << "test listL " << operatore->getSizeListL_test() << endl;
    cout << "test listG " << operatore->getSizeListG_test() << endl << endl;
    
    cout << "field numIndices " << operatore->numIndex_field() << endl;
    
    //Field data
    operatore->setElement(1);
    operatore->setIJ(1,1,1,1);
   
    indices.resize(operatore->numIndex_field());
    operatore->indexL_field(indices);
    
    cout << "field indices " << endl << indices << endl;
    
    //Test data
    operatore->setElement(1);
    operatore->setIJ(1,1,1,1);
   
    indices.resize(operatore->numIndex_test());
    operatore->indexL_test(indices);
    
    cout << "test indices " << endl << indices << endl;
    
    //Test eval
    sVect<Real> mat(operatore->numIndex_test() * operatore->numIndex_field());
    operatore->eval(point3d(1.0/3.0, 1.0/3.0, 0.0), mat);
    
    cout << "eval " << endl << mat << endl;
  }  
}
