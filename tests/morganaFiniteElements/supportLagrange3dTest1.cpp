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
#include "feRt0LT3d.hpp"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "supportLagrange3d.hpp"


using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef pMapItem             PMAPTYPE;
  typedef fePr3d<1,PMAPTYPE>   FIELD_FETYPE;
  typedef Real                 FIELD_DOFTYPE;
  typedef fePr3d<1,PMAPTYPE>   TEST_FETYPE;
  typedef Real                 TEST_DOFTYPE;
  typedef fePr3d<0,PMAPTYPE>   COEFF_FETYPE;
  typedef Real                 COEFF_DOFTYPE;

  typedef typename FIELD_FETYPE::GEOSHAPE   GEOSHAPE;
  
  typedef meshInit3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> MESHINIT;
  typedef typename MESHINIT::MESH3D     MESH3D;
  typedef typename MESHINIT::CONNECT3D  CONNECT3D;
  typedef typename MESHINIT::MESH2D     MESH2D;
  typedef typename MESHINIT::CONNECT2D  CONNECT2D;

  
  //Loading
  string meshFile = "./tests/morganaMeshes/mignon3dC.msh";
  
  MESHINIT init(world);
  init.gmMesh_to_stdB(meshFile, false);
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  typedef supportLagrange3d<GEOSHAPE,PMAPTYPE> SUPPORT;
  SUPPORT support(grid3d, connectGrid3d, grid2d, connectGrid2d);
  
  if(world.rank() == 0)
  {
    UInt el    = 2;
    UInt locFc = 2;
    point3d Ys(0.1,0.33,0.0);
    
    support.setElement3d(el);
    support.setLocalFace(locFc);
    
    cout << grid3d->getElementL(el) << endl;
    cout << grid3d->getNodes() << endl;
    
    cout << "isBoundary : " << support.isBoundary() << endl;
    cout << "el2d       : " << support.getElement2d() << endl;
    cout << "Yv         : " << support.mapVolumeY(Ys);
    cout << "N          : " << support.computeNormal(Ys);
    cout << "Y2d        : " << support.mapSurfaceY(Ys);
  }
}
