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

#include "meshDoctor2d.hpp"
#include "meshInit3d.hpp"
#include "meshInit2d.hpp"

#include "feStaticField1dGlobalManip.hpp"
#include "feStaticField2dGlobalManip.hpp"

#include "feStaticFieldPrinter2d.hpp"
#include "printMesh1dHDF5.hpp"
#include "printMesh2dHDF5.hpp"

#include "feStaticFieldPrinter1d.hpp"
#include "feStaticFieldPrinter2d.hpp"

#include "fvEL_linearADV2dA.hpp"

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  //Startup environment----------------------------------------------
  environment  env(argc,argv);
  communicator world;
  
  //Active geoIds----------------------------------------------------
  sVect<UInt> geoIds2d, geoIds1d;
  
  geoIds2d.push_back(1);
  geoIds2d.push_back(2);
  geoIds2d.push_back(3);
  geoIds2d.push_back(4);
  geoIds2d.push_back(5);
  geoIds2d.push_back(6);
  
  
  //Typedefs choose the fields---------------------------------------
  typedef linearTriangle         GEOSHAPE2D;
  typedef GEOSHAPE2D::GEOBSHAPE  GEOSHAPE1D;
  typedef pMapItemShare          PMAPTYPE;
  
  typedef fvEL_linearADV2dA<GEOSHAPE2D> SOLVER;
  
  typedef typename SOLVER::BOUNDARY_DOFTYPE  BOUNDARY_DOFTYPE;
  typedef typename SOLVER::SOLUTION_DOFTYPE  SOLUTION_DOFTYPE;
  typedef typename SOLVER::FLUXPARM_DOFTYPE  FLUXPARM_DOFTYPE;
  
  typedef typename SOLVER::BOUNDARY_FETYPE  BOUNDARY_FETYPE;
  typedef typename SOLVER::SOLUTION_FETYPE  SOLUTION_FETYPE;
  typedef typename SOLVER::FLUXPARM_FETYPE  FLUXPARM_FETYPE;
  
  typedef typename SOLVER::BOUNDARY_FIELD  BOUNDARY_FIELD;
  typedef typename SOLVER::SOLUTION_FIELD  SOLUTION_FIELD;
  typedef typename SOLVER::FLUXPARM_FIELD  FLUXPARM_FIELD;
  
  typedef BOUNDARY_FIELD::OPTIONS  BOUNDARY_OPTIONS;
  typedef SOLUTION_FIELD::OPTIONS  SOLUTION_OPTIONS;
  typedef FLUXPARM_FIELD::OPTIONS  FLUXPARM_OPTIONS;
  
 
  //Typedefs Loading-------------------------------------------------  
  typedef meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> INIT2D;
  typedef INIT2D::MESH1D    MESH1D;
  typedef INIT2D::CONNECT1D CONNECT1D;
  
  
  //Typedefs Loading-------------------------------------------------  
  typedef meshInit3d<linearTetra,PMAPTYPE,PMAPTYPE> INIT3D;
  typedef INIT3D::MESH3D    MESH3D;
  typedef INIT3D::CONNECT3D CONNECT3D;
  
  typedef INIT3D::MESH2D    MESH2D;
  typedef INIT3D::CONNECT2D CONNECT2D;
  
  
  //Loading and downloading------------------------------------------
  //string meshFile = "./tests/morganaMeshes/mignon3dH.msh";
  string meshFile = "./geometries/cubes3d/testCubeC.msh";

  INIT3D  init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH3D>           grid3d = init.getGrid3d();
  RCP<CONNECT3D> connectGrid3d = init.getConnectGrid3d();
  
  RCP<MESH2D>           grid2d = init.getGrid2d();
  RCP<CONNECT2D> connectGrid2d = init.getConnectGrid2d();

  //Testing----------------------------------------------------------
  bool internal, owned;
  UInt ed1, ed2, ed3;
  
  if(world.rank() == 0)
  {
    //cout << grid2d->getElements() << endl;
    
    
    for(UInt el=1; el <= grid2d->getNumElements(); ++el)
    {
      ed1 = connectGrid2d->getElementToEdge(el,1);
      ed2 = connectGrid2d->getElementToEdge(el,2);
      ed3 = connectGrid2d->getElementToEdge(el,3);
      
      owned = grid2d->getElements().getRowMapL(el).getOwned();
      
      internal = (connectGrid2d->getNumEdgeToElement(ed1) == 2) &&
                 (connectGrid2d->getNumEdgeToElement(ed2) == 2) &&
                 (connectGrid2d->getNumEdgeToElement(ed3) == 2);
 
      if(owned && (!internal))
      {
	cout << "error " << endl;
	cout << "el : " << el << endl;
      }
    }
  }
}
