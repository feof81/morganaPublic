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

#include "feSpectralLH3d.hpp"
#include "feDynamicField3d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  UInt pid = world.rank();
  
  
  //Mesh Init------------------------------------------------------------------
  typedef pMapItemShare             PMAPTYPE;
  typedef point3d                   DOFTYPE;
  typedef feSpectralLH3d<PMAPTYPE>  FETYPE;
  typedef typename FETYPE::GEOSHAPE GEOSHAPE;
  
  string meshFile  = "../tests/morganaMeshes/mignonHexaB.unv";
  string colorFile = "../tests/morganaMeshes/mignonHexaB_color.unv";
  
  meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.femap_to_stdB(meshFile, colorFile, false);
  
  
  //Download data--------------------------------------------------------------
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH3D     MESH3D;
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT3D  CONNECT3D;
  
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  
  //Initialize the mapper
  dofMapDynamic3d_options mapOptions;
  feDynamicInterface3d<FETYPE,DOFTYPE,dmd3d_componentMajor,dmd3d_standard> interface;
  
  interface.setCommunicator(world);
  interface.setGeometry(grid3d,connectGrid3d);
  interface.setOptions(mapOptions);
  
  
  //Load the feCards
  typedef typename FETYPE::FECARD FECARD;
  FECARD FeCard;
  
  if(world.rank() == 0)
  {
    //el 1
    FeCard.setIsActive(true);
    FeCard.setR(2,1,1);
    interface.setFeCardL(1,FeCard);
    
    //el 2
    FeCard.setIsActive(true);
    FeCard.setR(1,2,1);
    interface.setFeCardL(2,FeCard);
  }
  else
  {
    //el 1
    FeCard.setIsActive(true);
    FeCard.setR(1,1,2);
    interface.setFeCardL(1,FeCard);
    
    //el 2
    FeCard.setIsActive(true);
    FeCard.setR(2,1,1);
    interface.setFeCardL(2,FeCard);
  }
  
  //Startup
  interface.startup();
  
  
  //Interace test 
  if(world.rank() == 0)
  {
    typedef feDynamicInterface3d<FETYPE,DOFTYPE,dmd3d_componentMajor,dmd3d_standard>::OUTTYPE  OUTTYPE;
    
    UInt el = 1;
    UInt I = 3, J = 1;
    point3d Y(0.23, 0.27, 0.55);
    
    cout << "Num basis" << endl;
    cout << interface.feNumBasisL(el) << endl;
  
    sVect<UInt>    indices(interface.feNumBasisL(el));
    sVect<OUTTYPE> eval(interface.feNumBasisL(el));
    sVect<OUTTYPE> gradX(interface.feNumBasisL(el)), gradY(interface.feNumBasisL(el)), gradZ(interface.feNumBasisL(el));
  
    interface.feIndicesLG(el,I,J,indices);
    interface.feEvalL(el,I,J, Y, eval);
    interface.feGradL(el,I,J, Y, gradX, gradY, gradZ);
    
    cout << "Indices" << endl;
    cout << indices << endl;
    
    
    point3d totEval, totGradX, totGradY, totGradZ;
    
    for(UInt i=1; i <= eval.size(); ++i)
    {
      totEval  += eval(i);
      totGradX += gradX(i);
      totGradY += gradY(i);
      totGradZ += gradZ(i);
    }
    
    cout << "Eval" << endl;
    cout << totEval << endl;
    
    cout << "GradX" << endl;
    cout << totGradX << endl;
    
    cout << "GradY" << endl;
    cout << totGradY << endl;
    
    cout << "GradZ" << endl;
    cout << totGradZ << endl;
  }
}
