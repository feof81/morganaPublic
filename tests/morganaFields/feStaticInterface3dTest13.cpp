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

#include "feRt0LT3d.hpp"
#include "feRt0LT3d_extern.hpp"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"
#include "feStaticInterface3d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef pMapItem                   PMAPTYPE;
  typedef feRt0LT3d<PMAPTYPE>        FETYPE;
  typedef Real                       DOFTYPE;
  typedef typename FETYPE::GEOSHAPE  GEOSHAPE;
  typedef typename FETYPE::FECARD    FECARD;
  typedef typename FETYPE::ELCARD    ELCARD;
  typedef typename FETYPE::BASETYPE  BASETYPE;
  typedef typename feStaticInterface3d<FETYPE,DOFTYPE>::OUTTYPE  OUTTYPE;
  
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> MESHINIT;
  typedef typename MESHINIT::MESH3D     MESH3D;
  typedef typename MESHINIT::CONNECT3D  CONNECT3D;
  
  
  //Init mesh
  MESHINIT init(world);
  string meshFile = "../tests/morganaMeshes/mignon3dC.msh";
  init.gmMesh_to_stdB(meshFile, false);
  
  
  //Download info
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  
  //Testing
  typedef feRt0LT3d_extern<PMAPTYPE>  EXTERNAL;
  typedef typename EXTERNAL::FECARDS  FECARDS;
  
  EXTERNAL external(grid3d,connectGrid3d);
  external.setCommDev(world);
  FECARDS cards = external.buildFeCards();
  
  
  //Data
  UInt I  = 1;
  UInt J  = 1;
  UInt el = 3;
  point3d Y(0.0, 1.0, 0.0);
  

  //Eval   
  dofMapStatic3d_options mapOptions;
  mapOptions.addGeoId(1);
  
  feStaticInterface3d<FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_geoIdMode> interface;
  interface.setOptions(mapOptions);
  interface.setCommunicator(world);
  interface.setGeometry(grid3d,connectGrid3d);
  interface.setFeCards(cards);
  interface.startup();
  
  if(world.rank() == 0)
  {
    cout << "Elements" << endl;
    cout << grid3d->getElements() << endl;
  }
  
  UInt numBasis = interface.feNumBasisL(el);
  sVect<OUTTYPE> globalEval(numBasis), globalGradX(numBasis), globalGradY(numBasis), globalGradZ(numBasis);
  sVect<UInt>    indices(numBasis);  
  
  //Evaluation
  if(world.rank() == 0)
  {
    interface.feIndicesLL(el,I,J,indices);
    interface.feEvalL(el,I,J,Y,globalEval);
    interface.feGradL(el,I,J,Y,globalGradX,globalGradY,globalGradZ);

    cout << "NumBasis : " << numBasis << endl << endl;
    
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

