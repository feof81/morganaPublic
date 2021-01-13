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
#include "feRt0LT3d_extern.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "mesh2dGlobalManip.hpp"
#include "mesh3dGlobalManip.hpp"

#include "interpolatorNodal.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  //Communicator-----------------------------------------------------
  environment  env(argc,argv);
  communicator world;
  
  //Typedefs---------------------------------------------------------
  typedef linearTetra    GEOSHAPE3D;
  typedef linearTriangle GEOSHAPE2D;
  typedef pMapItemShare  PMAPTYPE;
  
  typedef mesh3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>     MESH3D;
  typedef mesh2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>     MESH2D;
  typedef connect3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>  CONNECT3D;
  typedef connect2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>  CONNECT2D;
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> MESHINIT;
  
  typedef feRt0LT3d<PMAPTYPE>  FETYPE;
  typedef Real                 DOFTYPE;
  typedef feStaticField3d<FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_allMode>            FIELD;
  typedef feStaticField3dGlobalManip<FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_allMode> FIELDMANIP;
  
  typedef FIELD::OPTIONS FIELD_OPTIONS;
  
  typedef feRt0LT3d_extern<PMAPTYPE>  EXTERNAL;
  typedef typename EXTERNAL::FECARDS  FECARDS;
  
  typedef mesh2dGlobalManip<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> MESHMANIP2D;
  typedef mesh3dGlobalManip<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> MESHMANIP3D;
  
  
  //Grids------------------------------------------------------------
  string meshFile = "./geometries/cubes3d/testCubeE.msh";
  //string meshFile = "./tests/morganaMeshes/mignon3dD.msh";
  
  MESHINIT init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH3D>       srcGrid3d = init.getGrid3d();
  RCP<MESH2D>       srcGrid2d = init.getGrid2d();
  RCP<CONNECT3D> srcConnect3d = init.getConnectGrid3d();
  RCP<CONNECT2D> srcConnect2d = init.getConnectGrid2d();
  
  //Create new grid--------------------------------------------------
  RCP<MESH3D> tgtGrid3d(new MESH3D(*srcGrid3d));
  RCP<MESH2D> tgtGrid2d(new MESH2D(*srcGrid2d));
  
  MESHMANIP3D meshManip3d(world);
  meshManip3d.gather(0,tgtGrid3d,srcGrid3d);
  
  MESHMANIP2D meshManip2d(world);
  meshManip2d.gather(0,tgtGrid2d,srcGrid2d);
  
  RCP<CONNECT3D> tgtConnect3d(new CONNECT3D(world));
  tgtConnect3d->setMesh3d(tgtGrid3d);
  tgtConnect3d->setMesh2d(tgtGrid2d);
  tgtConnect3d->buildConnectivity();
  tgtConnect3d->buildBoundaryConnectivity();
  
  RCP<CONNECT2D> tgtConnect2d(new CONNECT2D(world));
  tgtConnect2d->setMesh2d(tgtGrid2d);
  tgtConnect2d->buildConnectivity();
  
  //RT0Cards SRC-----------------------------------------------------
  EXTERNAL srcExternal(srcGrid3d,srcConnect3d);
  srcExternal.setCommDev(world);
  FECARDS  srcFeCards = srcExternal.buildFeCards();

  //RT0Cards TGT-----------------------------------------------------
  EXTERNAL tgtExternal(tgtGrid3d,tgtConnect3d);
  tgtExternal.setCommDev(world);
  FECARDS  tgtFeCards = tgtExternal.buildFeCards();
  
  //Source field-----------------------------------------------------
  RCP<FIELD_OPTIONS> srcFieldOptions(new FIELD_OPTIONS);
  
  RCP<FIELD> srcField(new FIELD);
  srcField->setCommunicator(world);
  srcField->setGeometry(srcGrid3d,srcConnect3d);
  srcField->setFeCards(srcFeCards);
  srcField->setOptions(srcFieldOptions);
  srcField->startup();
  
  assert(srcField->getDofVect().size() == srcGrid3d->getNumFaces());
  sVect<point3d> faceNodes;
  Real dof;
  
  for(UInt f=1; f <= srcGrid3d->getNumFaces(); ++f)
  {
    faceNodes = srcGrid3d->getFaceNodesL(f);
    dof       = (faceNodes(1).getX() + faceNodes(2).getX() + faceNodes(3).getX()) / 3.0;
    srcField->setDofL(f,dof);
  }
  
  //Target field-----------------------------------------------------
  RCP<FIELD_OPTIONS> tgtFieldOptions(new FIELD_OPTIONS);
  RCP<FIELD> tgtField(new FIELD(*srcField));
  
  //Test-------------------------------------------------------------
  FIELDMANIP fieldManip(world);
  fieldManip.matchGrid(tgtGrid3d,
                       tgtConnect3d,
                       tgtFieldOptions,
                       tgtField);
  
  tgtField->setFeCards(tgtFeCards);
  
  //Evaluate fields--------------------------------------------------
  typedef point3d            INTERP_DOFTYPE;
  typedef fePr3d<0,PMAPTYPE> INTERP_FETYPE;
  typedef feStaticField3d<INTERP_FETYPE,INTERP_DOFTYPE,dms3d_vectMajor,dms3d_allMode> INTERP_FIELD;
  
  dofMapStatic3d_options interpOptions;
  
  INTERP_FIELD interpSrcField;
  interpSrcField.setCommunicator(world);
  interpSrcField.setGeometry(srcGrid3d,srcConnect3d);
  interpSrcField.setOptions(interpOptions);
  interpSrcField.startup();

  INTERP_FIELD interpTgtField;
  interpTgtField.setCommunicator(world);
  interpTgtField.setGeometry(srcGrid3d,srcConnect3d);
  interpTgtField.setOptions(interpOptions);
  interpTgtField.startup();
  
  interpolatorNodal<FIELD,INTERP_FIELD> interpolatorSrc(world);
  interpolatorSrc.setMesh(srcGrid2d,
                          srcGrid3d,
                          srcConnect2d,
                          srcConnect3d);
  
  interpolatorSrc.localInit(2.0);
  interpolatorSrc.globalInit();
  interpolatorSrc.findDofs(interpSrcField);
  interpolatorSrc.exchangeData(*srcField,interpSrcField);
  
  interpolatorNodal<FIELD,INTERP_FIELD> interpolatorTgt(world);
  interpolatorTgt.setMesh(tgtGrid2d,
                          tgtGrid3d,
                          tgtConnect2d,
                          tgtConnect3d);
  
  interpolatorTgt.localInit(2.0);
  interpolatorTgt.globalInit();
  interpolatorTgt.findDofs(interpTgtField);
  interpolatorTgt.exchangeData(*tgtField,interpTgtField);
  
  assert(interpSrcField.getDofVect().size() == interpTgtField.getDofVect().size());
  
  //Compute residual-------------------------------------------------
  point3d V;
  Real res = 0.0;
  
  for(UInt i=1; i <= interpSrcField.getDofVect().size(); ++i)
  {
    V    = interpSrcField.getDofL(i) - interpTgtField.getDofL(i);
    res += point3d::norm2(V); 
  }
  
  cout << "residual : " << res << endl;
}