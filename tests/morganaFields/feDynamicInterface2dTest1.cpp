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

#include "feSpectralLH2d.hpp"
#include "feDynamicField2d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef linearQuad     GEOSHAPE;
  typedef pMapItemShare  PMAPTYPE;  
  
  string meshFile  = "../tests/morganaMeshes/mignonQuad2dB.unv";
  string colorFile = "../tests/morganaMeshes/mignonQuad2dB_color.unv";
  
  meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> init(world);
  init.femap_to_stdB(meshFile, colorFile, false);
  
  
  //Download data--------------------------------------------------------------
  typedef meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::MESH2D     MESH2D;
  typedef meshInit2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>::CONNECT2D  CONNECT2D;
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT2D> connectGrid3d = init.getConnectGrid2d();
  
  
  //Types
  typedef point3d                  DOFTYPE;
  typedef feSpectralLH2d<PMAPTYPE> FETYPE;
  
  
  //Initialize the mapper
  dofMapDynamic2d_options mapOptions;
  feDynamicInterface2d<FETYPE,DOFTYPE,dmd2d_componentMajor,dmd2d_standard> interface;
  
  interface.setCommunicator(world);
  interface.setGeometry(grid2d,connectGrid3d);
  interface.setOptions(mapOptions);
  
  
  //Load the feCards
  typedef typename FETYPE::FECARD FECARD;
  FECARD FeCard;
  
  FeCard.setIsActive(true);
  FeCard.setR(2,1);
  interface.setFeCardL(1,FeCard);
  
  FeCard.setIsActive(true);
  FeCard.setR(1,2);
  interface.setFeCardL(2,FeCard);
  
  //Startup
  interface.startup();
  
 
  //Interace test 
  if(world.rank() == 0)
  {
    typedef feDynamicInterface2d<FETYPE,DOFTYPE,dmd2d_componentMajor,dmd2d_standard>::OUTTYPE  OUTTYPE;
    
    UInt el = 1;
    UInt I = 3, J = 1;
    point3d Y(0.23, 0.17, 0.0);
    
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
