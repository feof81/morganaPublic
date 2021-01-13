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
#include "pMap.hpp"

#include "pVect.hpp"
#include "pVectComm.hpp"

#include "fePr3d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"
#include "feStaticFieldPrinter3d.hpp"
#include "meshInit3d.hpp"

#include "traitsMpiOptimization.hpp"

using namespace Teuchos;


//! Run with two processors
int main(int argc, char *argv[])
{
  //Typedefs---------------------------------------------------------
  typedef linearTriangle         GEOSHAPE2D;
  typedef linearTetra            GEOSHAPE3D;
  typedef geoElement<GEOSHAPE3D> GEOELEMENT;
  typedef pMapItemShare          PMAPTYPE;
  
  typedef Real               DOFTYPE;
  typedef fePr3d<1,PMAPTYPE> FETYPE;
  
  typedef meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>                                 INIT;
  typedef feStaticField3d<FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_allMode>            FIELD;
  typedef mesh3dGlobalManip<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>                          MANIP_MESH;
  typedef feStaticField3dGlobalManip<FETYPE,DOFTYPE,dms3d_vectMajor,dms3d_allMode> MANIP_FIELD;
  
  typedef INIT::MESH2D     MESH2D;
  typedef INIT::MESH3D     MESH3D;
  typedef INIT::CONNECT2D  CONNECT2D;
  typedef INIT::CONNECT3D  CONNECT3D;
  
  //Comm dev---------------------------------------------------------
  boost::mpi::environment  env(argc, argv);
  boost::mpi::communicator oldCommDev;
  
  bool commActive = (oldCommDev.rank() < 1);
  boost::mpi::communicator newCommDev = oldCommDev.split(commActive);
  
  assert(oldCommDev.size() == 2);
  
  //Old mesh---------------------------------------------------------
  string meshFile = "./tests/morganaMeshes/mignon3dC.msh";
  
  INIT oldInit(oldCommDev);
  oldInit.gmMesh_to_stdA(meshFile);
  
  RCP<MESH3D>       oldGrid3d = oldInit.getGrid3d();
  RCP<CONNECT3D> oldConnect3d = oldInit.getConnectGrid3d();
  
  //Old field--------------------------------------------------------
  dofMapStatic3d_options oldMapOptions;
  
  FIELD oldField;
  oldField.setCommunicator(oldCommDev);
  oldField.setGeometry(oldGrid3d,oldConnect3d);
  oldField.setOptions(oldMapOptions);
  oldField.startup();
   
  for(UInt i=1; i <= oldField.size(); ++i)
  { oldField.setDofL(i, oldGrid3d->getNodeL(i).getZ() ); }
  
  feStaticFieldPrinter3d<FIELD> oldPrinter(oldCommDev);
  oldPrinter.printHDF5_nodes("oldField",0,oldField);
  
  //New grid---------------------------------------------------------
  meshInit3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> newInit(newCommDev);
  RCP<MESH3D>    newGrid3d;
  RCP<CONNECT3D> newConnect3d;
    
  if(commActive)
  {
    newInit.gmMesh_to_stdA(meshFile);
    
    newGrid3d    = newInit.getGrid3d();
    newConnect3d = newInit.getConnectGrid3d();
  }
  
  //Shifted grid-----------------------------------------------------
  RCP<MESH3D> shiftGrid3d(new MESH3D);
  
  MANIP_MESH shiftMeshManip;
  shiftMeshManip.expandCommunicator(commActive,
                                    newCommDev,
                                   *newGrid3d,
                                    oldCommDev,
                                   *shiftGrid3d);
  
  RCP<CONNECT3D> shiftConnect3d(new CONNECT3D(oldCommDev));
  shiftConnect3d->setMesh3d(shiftGrid3d);
  shiftConnect3d->buildConnectivity();
  
  for(UInt pid=0; pid < oldCommDev.size(); pid++)
  {
    if(pid == oldCommDev.rank())
    {
      cout << "Process : " << pid << " -------------------------- " << endl;
      //cout << shiftGrid3d->getNodes()       << endl;
      cout << shiftGrid3d->getElements()    << endl;
      //cout << shiftGrid3d->getFaces()       << endl;
      //cout << shiftGrid3d->getEdges()       << endl;
    }
    
    oldCommDev.barrier();
  }
  
  //Shifted field----------------------------------------------------
  dofMapStatic3d_options shiftOptions;
  FIELD shiftField(oldField);

  MANIP_FIELD shiftManip(oldCommDev);
  shiftManip.matchGrid(*shiftGrid3d,
                       *shiftConnect3d,
                        shiftOptions,
                        shiftField);
  
  for(UInt pid=0; pid < oldCommDev.size(); pid++)
  {
    if(pid == oldCommDev.rank())
    {
      cout << "Process : " << pid << " -------------------------- " << endl;
      cout << shiftField.getDofVect() << endl;
    }
    
    oldCommDev.barrier();
  }
    
  //New comm---------------------------------------------------------
  if(commActive)
  {
    FIELD newField;
    
    MANIP_FIELD newManip;
    newManip.reduceCommunicator(commActive,
                                oldCommDev,
                                shiftField,
                                newCommDev,
                                newGrid3d,
                                newConnect3d,
                                newField);
    
    feStaticFieldPrinter3d<FIELD> newPrinter(newCommDev);
    newPrinter.printHDF5_nodes("newField",0,newField);
  }
  
}
