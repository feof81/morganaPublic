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
  
  typedef linearTriangle  GEOSHAPE;
  typedef pMapItemShare   PMAPTYPE;
  
  string meshFile = "../tests/morganaMeshes/mignon2dB.msh";
  //string meshFile = "./geometries/cubes3d/testCubeD.msh";
  
  meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  
  //Interface test 
  typedef point3d            DOFTYPE;
  typedef fePr2d<1,PMAPTYPE> FETYPE;
  
  dofMapStatic2d_options mapOptions;  
  feStaticInterface2d<FETYPE,DOFTYPE,dms2d_componentMajor,dms2d_allMode> interface;
  
  interface.setCommunicator(world);
  interface.setGeometry(init.getGrid2d(), init.getConnectGrid2d());
  interface.setOptions(mapOptions);
  interface.startup();
  
  //Interace test 
  if(world.rank() == 0)
  {
    typedef feStaticInterface2d<FETYPE,DOFTYPE,dms2d_componentMajor,dms2d_allMode>::OUTTYPE  OUTTYPE;
    
    UInt el = 1;
    UInt I = 3, J = 1;
    point3d Y(1.0/3.0, 1.0/3.0, 0.0);
    
    sVect<UInt>    indices(interface.feNumBasisL(el));
    sVect<OUTTYPE> eval(interface.feNumBasisL(el));
    sVect<OUTTYPE> gradX(interface.feNumBasisL(el)), gradY(interface.feNumBasisL(el)), gradZ(interface.feNumBasisL(el));
  
    interface.feIndicesLG(el,I,J,indices);
    interface.feEvalL(el,I,J, Y, eval);
    interface.feGradL(el,I,J, Y, gradX, gradY, gradZ);
    
    cout << "Indices" << endl;
    cout << indices << endl;
    
    cout << "Eval" << endl;
    cout << eval << endl;
    
    cout << "GradX" << endl;
    cout << gradX << endl;
    
    cout << "GradY" << endl;
    cout << gradY << endl;
    
    cout << "GradZ" << endl;
    cout << gradZ << endl;
  }
}
