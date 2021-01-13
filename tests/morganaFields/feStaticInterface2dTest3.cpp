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

#include "feD0LT2d.hpp"
#include "dofMapStatic2d.hpp"
#include "elCardFeeder2d.hpp"
#include "feStaticInterface2d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef pMapItemShare              PMAPTYPE;
  typedef feD0LT2d<PMAPTYPE>         FETYPE;
  typedef Real                       DOFTYPE;
  typedef typename FETYPE::GEOSHAPE  GEOSHAPE;
  typedef typename FETYPE::FECARD    FECARD;
  typedef typename FETYPE::ELCARD    ELCARD;
  typedef typename FETYPE::BASETYPE  BASETYPE;
  typedef typename feStaticInterface2d<FETYPE,DOFTYPE>::OUTTYPE  OUTTYPE;
  
  string meshFile = "../tests/morganaMeshes/mignon2dB.msh";
  //string meshFile = "./geometries/cubes3d/testCubeD.msh";
  
  meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  
  //Data
  UInt el = 2;
  point3d Y(0.0, 1.0, 0.0);
  
  //Cards
  elCardFeeder2d<GEOSHAPE,PMAPTYPE> cardsFeeder(init.getGrid2d(),init.getConnectGrid2d());
  
  FECARD FeCard;
  ELCARD ElCard = cardsFeeder.getCardLocal(el);
  
  
  //Finite element
  FETYPE finiteElement(FeCard,ElCard);
  
  if(world.rank() == 0)
  {
    cout << "El Card --------------------" << endl;
    cout << ElCard << endl;
  }
  

  //Eval 
  feStaticInterface2d<FETYPE,DOFTYPE> interface;
  bool flag = interface.feTest();
  assert(flag);
  
  UInt numBasis = FETYPE::numBasis;
  sVect<OUTTYPE> globalEval(numBasis), globalGradX(numBasis), globalGradY(numBasis), globalGradZ(numBasis);
  
  //Dof indices
  UInt I=1;
  UInt J=1;
  
  //Evaluation
  if(world.rank() == 0)
  {
    interface.eval(el, I,J, FeCard,ElCard, Y, globalEval);
    interface.grad(el, I,J, FeCard,ElCard, Y, globalGradX,globalGradY,globalGradZ);
    
    cout << "GlobalEval" << endl;
    cout << globalEval << endl;
    
    cout << "GlobalGradX" << endl;
    cout << globalGradX << endl;
    
    cout << "GlobalGradY" << endl;
    cout << globalGradY << endl;
    
    cout << "GlobalGradZ" << endl;
    cout << globalGradZ << endl;
  }
}
