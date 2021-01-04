/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MESHINIT1D_HPP
#define MESHINIT1D_HPP

#include <map>
#include <set>

#include "staticAssert.hpp"
#include "time.h"

#include "pMapItem.h"
#include "pMapGlobalManip.h"

#include "pGraph.hpp"
#include "pGraphItem.h"
#include "pGraphItemOriented.h"
#include "pGraphItemSubLoc.h"
#include "pGraphGlobalManip.hpp"
#include "pGraphManip.hpp"

#include "loadMesh.hpp"
#include "connect1d.hpp"
#include "mesh1dGlobalManip.hpp"
#include "meshDoctor1d.hpp"


/*! Mesh initialization 1d */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class meshInit1d
{
   /*! @name Typedefs */ //@{
  public:
    typedef GEOSHAPE                            GEOSHAPE1D;
    typedef geoElement<GEOSHAPE1D>              GEOELEMENT1D;
    typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>    MESH1D;
    typedef connect1d<GEOSHAPE1D,ELMAP,NODEMAP> CONNECT1D;
    
    typedef typename MESH1D::NODESVECT  NODESVECT;
    typedef typename MESH1D::GRAPH1D    GRAPH1D;
    //@}
    
    /*! @name Control flags */ //@{
  public:
    bool commDevLoaded;
    bool mesh1dCreated;
    bool connect1dCreated;
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<communicator>  commDev;
    Teuchos::RCP<MESH1D>        grid1d;
    Teuchos::RCP<CONNECT1D>     connectGrid1d;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    meshInit1d();
    meshInit1d(const Teuchos::RCP<communicator> & CommDev);
    meshInit1d(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Get functions */ //@{
  public:
    Teuchos::RCP<MESH1D>    getGrid1d();
    Teuchos::RCP<CONNECT1D> getConnectGrid1d();
    
    Teuchos::RCP<MESH1D>    getCopyGrid1d();
    Teuchos::RCP<CONNECT1D> getCopyConnectGrid1d();
    //@}
    
    /*! @name Loading functions */ //@{
  public:
    void gmMesh_to_stdA(const string & meshfile, const bool & verbose = true);
    void gmMesh_to_stdB(const string & meshfile, const bool & verbose = true);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::
meshInit1d()
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  commDevLoaded    = false;
  mesh1dCreated    = false;
  connect1dCreated = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::
meshInit1d(const Teuchos::RCP<communicator> & CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  commDev = CommDev;
  
  commDevLoaded    = true;
  mesh1dCreated    = false;
  connect1dCreated = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::
meshInit1d(communicator & CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  commDev = Teuchos::rcpFromRef(CommDev);
  
  commDevLoaded    = true;
  mesh1dCreated    = false;
  connect1dCreated = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}


//_______________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::MESH1D>
meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::
getGrid1d()
{
  assert(mesh1dCreated);
  return(grid1d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::CONNECT1D>
meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::
getConnectGrid1d()
{
  assert(connect1dCreated);
  return(connectGrid1d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::MESH1D>
meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::
getCopyGrid1d()
{
  assert(mesh1dCreated);
  Teuchos::RCP<MESH1D> newGrid1d(new MESH1D(*grid1d));
  return(newGrid1d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::CONNECT1D>
meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::
getCopyConnectGrid1d()
{
  assert(connect1dCreated);
  Teuchos::RCP<CONNECT1D>  newConnectGrid1d(new CONNECT1D(*connectGrid1d));
  return(newConnectGrid1d);
}


//_______________________________________________________________________________________________________
// LOADING FUNCTIONS
//-------------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::
gmMesh_to_stdA(const string & meshfile, const bool & verbose)
{
  assert(staticAssert< GEOSHAPE::geoName == LINEARLINE >::returnValue);
  assert(commDevLoaded);
  
  //Allocations------------------------------------------------------------------------------------
  UInt pid      = commDev->rank();
  UInt numPids  = commDev->size();
  UInt printPid = 0;
  time_t start, end;
  
  //Allocations
  NODESVECT nodes;
  GRAPH1D   elements1d;
  
  //Mesh loading-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Start mesh loading" << " "; time(&start); cout << endl;}
  loadMesh<GEOSHAPE,ELMAP,NODEMAP> meshLoader;
  meshLoader.setParam(pid,printPid,numPids);
  meshLoader.gmMesh(meshfile,nodes,elements1d);
  
  elements1d.colIsLocal() = false;
  
  commDev->barrier();
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //The grid1d-------------------------------------------------------------------------------------
  grid1d = rcp(new mesh1d<GEOSHAPE,ELMAP,NODEMAP>);
  grid1d->setNodes(nodes);
  grid1d->setElements(elements1d);
  grid1d->transferMap();
  grid1d->setMeshStandard(STDL);
  
  //Mesh partition---------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh partition" << " "; time(&start); cout << endl;}
  
  mesh1dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip1d(commDev);
  gridManip1d.meshPartition(grid1d);
  assert(gridManip1d.check(grid1d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Mesh overlapping-------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh overlap" << " "; time(&start); cout << endl;}
  
  meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP> gridDoctor(commDev);
  gridManip1d.meshOverlapping(grid1d);
  gridDoctor.fixGeoIds(grid1d);
  assert(gridManip1d.check(grid1d));
  assert(gridDoctor.checkGeoIds(grid1d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Surface 1d" << " "; time(&start); cout << endl;}
  
  connectGrid1d = rcp(new connect1d<GEOSHAPE1D,ELMAP,NODEMAP>(commDev));
  connectGrid1d->setMesh1d(grid1d);
  connectGrid1d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Flags------------------------------------------------------------------------------------------
  grid1d->setMeshStandard(STDA);
  
  mesh1dCreated    = true;
  connect1dCreated = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit1d<GEOSHAPE,ELMAP,NODEMAP>::
gmMesh_to_stdB(const string & meshfile, const bool & verbose)
{
  assert(staticAssert< GEOSHAPE::geoName == LINEARLINE >::returnValue);
  assert(commDevLoaded);
  
  //Allocations------------------------------------------------------------------------------------
  UInt pid      = commDev->rank();
  UInt numPids  = commDev->size();
  UInt printPid = 0;
  time_t start, end;
  
  //Allocations
  NODESVECT nodes;
  GRAPH1D   elements1d;
  
  //Mesh loading-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Start mesh loading" << " "; time(&start); cout << endl;}
  loadMesh<GEOSHAPE,ELMAP,NODEMAP> meshLoader;
  meshLoader.setParam(pid,printPid,numPids);
  meshLoader.gmMesh(meshfile,nodes,elements1d);
  
  elements1d.colIsLocal() = false;
  
  commDev->barrier();
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //The grid1d-------------------------------------------------------------------------------------
  grid1d = rcp(new mesh1d<GEOSHAPE,ELMAP,NODEMAP>);
  grid1d->setNodes(nodes);
  grid1d->setElements(elements1d);
  grid1d->transferMap();
  grid1d->setMeshStandard(STDL);
  
  //Mesh partition---------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh partition" << " "; time(&start); cout << endl;}
  
  mesh1dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip1d(commDev);
  gridManip1d.meshPartition(grid1d);
  assert(gridManip1d.check(grid1d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Mesh cheking-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh overlap" << " "; time(&start); cout << endl;}
  
  meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP> gridDoctor(commDev);
  gridDoctor.fixGeoIds(grid1d);
  assert(gridManip1d.check(grid1d));
  assert(gridDoctor.checkGeoIds(grid1d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Surface 1d" << " "; time(&start); cout << endl;}
  
  connectGrid1d = rcp(new connect1d<GEOSHAPE1D,ELMAP,NODEMAP>(commDev));
  connectGrid1d->setMesh1d(grid1d);
  connectGrid1d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Flags------------------------------------------------------------------------------------------
  grid1d->setMeshStandard(STDB);
  
  mesh1dCreated    = true;
  connect1dCreated = true;
}


#endif
