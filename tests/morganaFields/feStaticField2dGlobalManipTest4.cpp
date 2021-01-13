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
#include "feStaticField2dGlobalManip.hpp"

#include "mesh1dGlobalManip.hpp"
#include "mesh2dGlobalManip.hpp"

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
  typedef linearTriangle GEOSHAPE2D;
  typedef linearLine     GEOSHAPE1D;
  typedef pMapItemShare  PMAPTYPE;
  
  typedef mesh2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>     MESH2D;
  typedef mesh1d<GEOSHAPE1D,PMAPTYPE,PMAPTYPE>     MESH1D;
  typedef connect2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>  CONNECT2D;
  typedef connect1d<GEOSHAPE1D,PMAPTYPE,PMAPTYPE>  CONNECT1D;
  typedef meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> MESHINIT;
  
  typedef feRt0LT2d<PMAPTYPE>  FETYPE;
  typedef Real                 DOFTYPE;
  typedef feStaticField2d<FETYPE,DOFTYPE,dms2d_vectMajor,dms2d_allMode>            FIELD;
  typedef feStaticField2dGlobalManip<FETYPE,DOFTYPE,dms2d_vectMajor,dms2d_allMode> FIELDMANIP;
  
  typedef FIELD::OPTIONS FIELD_OPTIONS;
  
  typedef feRt0LT2d_extern<PMAPTYPE>  EXTERNAL;
  typedef typename EXTERNAL::FECARDS  FECARDS;
  
  typedef mesh1dGlobalManip<GEOSHAPE1D,PMAPTYPE,PMAPTYPE> MESHMANIP1D;
  typedef mesh2dGlobalManip<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> MESHMANIP2D;
  
  
  //Grids------------------------------------------------------------
  string meshFile = "./geometries/rectangles2d/rectangleE.msh";
  
  MESHINIT init(world);
  init.gmMesh_to_stdA(meshFile);
  
  RCP<MESH2D>       srcGrid2d = init.getGrid2d();
  RCP<MESH1D>       srcGrid1d = init.getGrid1d();
  RCP<CONNECT2D> srcConnect2d = init.getConnectGrid2d();
  RCP<CONNECT1D> srcConnect1d = init.getConnectGrid1d();
  
  //Create new grid--------------------------------------------------
  RCP<MESH2D> tgtGrid2d(new MESH2D(*srcGrid2d));
  RCP<MESH1D> tgtGrid1d(new MESH1D(*srcGrid1d));
  
  MESHMANIP2D meshManip2d(world);
  meshManip2d.gather(0,tgtGrid2d,srcGrid2d);
  
  MESHMANIP1D meshManip1d(world);
  meshManip1d.gather(0,tgtGrid1d,srcGrid1d);
  
  RCP<CONNECT2D> tgtConnect2d(new CONNECT2D(world));
  tgtConnect2d->setMesh2d(tgtGrid2d);
  tgtConnect2d->setMesh1d(tgtGrid1d);
  tgtConnect2d->buildConnectivity();
  tgtConnect2d->buildBoundaryConnectivity();
  
  RCP<CONNECT1D> tgtConnect1d(new CONNECT1D(world));
  tgtConnect1d->setMesh1d(tgtGrid1d);
  tgtConnect1d->buildConnectivity();
  
  //RT0Cards SRC-----------------------------------------------------
  EXTERNAL srcExternal(srcGrid2d,srcConnect2d);
  srcExternal.setCommDev(world);
  FECARDS  srcFeCards = srcExternal.buildFeCards();

  //RT0Cards TGT-----------------------------------------------------
  EXTERNAL tgtExternal(tgtGrid2d,tgtConnect2d);
  tgtExternal.setCommDev(world);
  FECARDS  tgtFeCards = tgtExternal.buildFeCards();
  
  //Source field-----------------------------------------------------
  RCP<FIELD_OPTIONS> srcFieldOptions(new FIELD_OPTIONS);
  
  RCP<FIELD> srcField(new FIELD);
  srcField->setCommunicator(world);
  srcField->setGeometry(srcGrid2d,srcConnect2d);
  srcField->setFeCards(srcFeCards);
  srcField->setOptions(srcFieldOptions);
  srcField->startup();
  
  assert(srcField->getDofVect().size() == srcGrid2d->getNumEdges());
  sVect<point3d> edgeNodes;
  Real dof;
  
  for(UInt d=1; d <= srcGrid2d->getNumEdges(); ++d)
  {
    edgeNodes = srcGrid2d->getEdgeNodesL(d);
    dof       = (edgeNodes(1).getX() + edgeNodes(2).getX()) / 2.0;
    srcField->setDofL(d,dof);
  }
  
  //Target field-----------------------------------------------------
  RCP<FIELD_OPTIONS> tgtFieldOptions(new FIELD_OPTIONS);
  RCP<FIELD> tgtField(new FIELD(*srcField));
  
  //Test-------------------------------------------------------------
  FIELDMANIP fieldManip(world);
  fieldManip.matchGrid(tgtGrid2d,
                       tgtConnect2d,
                       tgtFieldOptions,
                       tgtField);
  
  tgtField->setFeCards(tgtFeCards);
  
  //Evaluate fields--------------------------------------------------
  typedef point3d            INTERP_DOFTYPE;
  typedef fePr2d<0,PMAPTYPE> INTERP_FETYPE;
  typedef feStaticField2d<INTERP_FETYPE,INTERP_DOFTYPE,dms2d_vectMajor,dms2d_allMode> INTERP_FIELD;
  
  dofMapStatic2d_options interpOptions;
  
  INTERP_FIELD interpSrcField;
  interpSrcField.setCommunicator(world);
  interpSrcField.setGeometry(srcGrid2d,srcConnect2d);
  interpSrcField.setOptions(interpOptions);
  interpSrcField.startup();

  INTERP_FIELD interpTgtField;
  interpTgtField.setCommunicator(world);
  interpTgtField.setGeometry(srcGrid2d,srcConnect2d);
  interpTgtField.setOptions(interpOptions);
  interpTgtField.startup();
  
  interpolatorNodal<FIELD,INTERP_FIELD> interpolatorSrc(world);
  interpolatorSrc.setMesh(srcGrid2d,
                          srcConnect2d);
  
  interpolatorSrc.localInit(2.0);
  interpolatorSrc.globalInit();
  interpolatorSrc.findDofs(interpSrcField);
  interpolatorSrc.exchangeData(*srcField,interpSrcField);
  
  interpolatorNodal<FIELD,INTERP_FIELD> interpolatorTgt(world);
  interpolatorTgt.setMesh(tgtGrid2d,
                          tgtConnect2d);
  
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