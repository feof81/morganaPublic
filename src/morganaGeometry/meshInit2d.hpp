/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MESHINIT2D_HPP
#define MESHINIT2D_HPP

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
#include "connect2d.hpp"

#include "mesh1dGlobalManip.hpp"
#include "mesh2dGlobalManip.hpp"
#include "meshDoctor1d.hpp"
#include "meshDoctor2d.hpp"


/*! Mesh initialization 2d */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class meshInit2d
{
   /*! @name Typedefs */ //@{
  public:
    typedef GEOSHAPE                             GEOSHAPE2D;
    typedef typename   GEOSHAPE::GEOBSHAPE       GEOSHAPE1D;
    
    typedef geoElement<GEOSHAPE2D>               GEOELEMENT2D;
    typedef geoElement<GEOSHAPE1D>               GEOELEMENT1D;
    
    typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>     MESH2D;
    typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>     MESH1D;
    
    typedef connect1d<GEOSHAPE1D,ELMAP,NODEMAP>  CONNECT1D;
    typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>  CONNECT2D;
    
    typedef typename MESH2D::NODESVECT           NODESVECT;
    typedef typename MESH2D::GRAPH2D             GRAPH2D;
    typedef typename MESH2D::GRAPH1D             GRAPH1D;
    //@}
    
    /*! @name Control flags */ //@{
  public:
    bool commDevLoaded;
    bool mesh2dCreated, mesh1dCreated;
    bool connect2dCreated, connect1dCreated;
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<communicator>  commDev;
    Teuchos::RCP<MESH2D>        grid2d;
    Teuchos::RCP<MESH1D>        grid1d;
    Teuchos::RCP<CONNECT2D>     connectGrid2d;
    Teuchos::RCP<CONNECT1D>     connectGrid1d;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    meshInit2d();
    meshInit2d(const Teuchos::RCP<communicator> & CommDev);
    meshInit2d(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Get functions */ //@{
  public:
    Teuchos::RCP<MESH2D>    getGrid2d();
    Teuchos::RCP<MESH1D>    getGrid1d();
    Teuchos::RCP<CONNECT2D> getConnectGrid2d();
    Teuchos::RCP<CONNECT1D> getConnectGrid1d();
    
    Teuchos::RCP<MESH2D>    getCopyGrid2d();
    Teuchos::RCP<MESH1D>    getCopyGrid1d();
    Teuchos::RCP<CONNECT2D> getCopyConnectGrid2d();
    Teuchos::RCP<CONNECT1D> getCopyConnectGrid1d();
    //@}
    
    /*! @name Loading functions */ //@{
  public:
    void gmMesh_to_stdA(const string & meshfile, const bool & verbose = true);
    void gmMesh_to_stdB(const string & meshfile, const bool & verbose = true);
    void femap_to_stdA(const string & meshfile, const string & colorfile, const bool & verbose = true);
    void femap_to_stdB(const string & meshfile, const string & colorfile, const bool & verbose = true);
    void neutral_to_stdA(const string & meshfile, const string & colorfile, const bool & verbose = true);
    void neutral_to_stdB(const string & meshfile, const string & colorfile, const bool & verbose = true);
    //@}
};


template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
meshInit2d()
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  
  commDevLoaded     = false;
  mesh2dCreated     = false;
  mesh1dCreated     = false;
  connect2dCreated  = false;
  connect1dCreated  = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
meshInit2d(const Teuchos::RCP<communicator> & CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  
  commDev = CommDev;
  
  commDevLoaded     = true;
  mesh2dCreated     = false;
  mesh1dCreated     = false;
  connect2dCreated  = false;
  connect1dCreated  = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
meshInit2d(communicator & CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  
  commDev = Teuchos::rcpFromRef(CommDev);
  
  commDevLoaded     = true;
  mesh2dCreated     = false;
  mesh1dCreated     = false;
  connect2dCreated  = false;
  connect1dCreated  = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}



//_______________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::MESH2D>
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
getGrid2d()
{
  assert(mesh2dCreated);
  return(grid2d);
}
 
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::MESH1D>
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
getGrid1d()
{
  assert(mesh1dCreated);
  return(grid1d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::CONNECT2D>
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
getConnectGrid2d()
{
  assert(connect2dCreated);
  return(connectGrid2d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::CONNECT1D>
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
getConnectGrid1d()
{
  assert(connect1dCreated);
  return(connectGrid1d);
}



template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::MESH2D>
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
getCopyGrid2d()
{
  assert(mesh2dCreated);
  Teuchos::RCP<MESH2D>  newGrid2d(new MESH2D(*grid2d));
  return(newGrid2d);
}
 
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::MESH1D>
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
getCopyGrid1d()
{
  assert(mesh1dCreated);
  Teuchos::RCP<MESH1D>  newGrid1d(new MESH1D(*grid1d));
  return(newGrid1d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::CONNECT2D>
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
getCopyConnectGrid2d()
{
  assert(connect2dCreated);
  Teuchos::RCP<CONNECT2D>  newConnectGrid2d(new CONNECT2D(*connectGrid2d));
  return(newConnectGrid2d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::CONNECT1D>
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
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
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
gmMesh_to_stdA(const string & meshfile, const bool & verbose)
{
  assert(staticAssert< GEOSHAPE::geoName == LINEARTRIANGLE >::returnValue);
  assert(commDevLoaded);
  
  //Allocations------------------------------------------------------------------------------------
  UInt pid      = commDev->rank();
  UInt numPids  = commDev->size();
  UInt printPid = 0;
  time_t start, end;
  
  //Allocations
  NODESVECT nodes;
  GRAPH2D   elements2d;
  GRAPH1D   elements1d;
  
  
  //Mesh loading-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Start mesh loading" << " "; time(&start); cout << endl;}
  loadMesh<GEOSHAPE,ELMAP,NODEMAP> meshLoader;
  meshLoader.setParam(pid,printPid,numPids);
  meshLoader.gmMesh(meshfile,nodes,elements2d,elements1d);
  
  elements1d.colIsLocal() = false;
  elements2d.colIsLocal() = false;
  
  commDev->barrier();
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //The grid3d-------------------------------------------------------------------------------------
  grid2d = rcp(new mesh2d<GEOSHAPE,ELMAP,NODEMAP>);
  grid2d->setNodes(nodes);
  grid2d->setElements(elements2d);
  grid2d->transferMap();
  grid2d->setMeshStandard(STDL);
  
  
  //Mesh partition---------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh partition" << " "; time(&start); cout << endl;}
  
  mesh2dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip2d(commDev);
  gridManip2d.meshPartition(grid2d);
  assert(gridManip2d.check(grid2d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh overlapping-------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh overlap" << " "; time(&start); cout << endl;}
  
  meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP> gridDoctor(commDev);
  gridManip2d.meshOverlapping(grid2d);
  gridDoctor.fixGeoIds(grid2d);
  assert(gridManip2d.check(grid2d));
  assert(gridDoctor.checkGeoIds(grid2d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Surface 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d = rcp(new connect2d<GEOSHAPE2D,ELMAP,NODEMAP>(commDev));
  connectGrid2d->setMesh2d(grid2d);
  connectGrid2d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Edges cheking---------------------------------------------------------------------------------- 
  if((pid == printPid) && verbose) {cout << "Check Edges" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(commDev);
  assert(checkEdges.check(grid2d->getEdges()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Distribute the elements2d----------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Distribute elements 1d" << " "; time(&start); cout << endl;}
    
  //Find local cadidates - the ones that has only one element connected
  GEOELEMENT1D  edgeGeoItem(true);
  ELMAP         edgeMapItem;
  
  GRAPH1D tempEdges2d = grid2d->getEdges();
  tempEdges2d.pushToGlobal();
  
  GRAPH1D edges2d = grid2d->getEdges();
  edges2d.pushToGlobal();  
  
  //Normal communication
  GEOELEMENT1D maxEdge(true), minEdge(true);
  GEOELEMENT1D maxEl1d(true), minEl1d(true);
  GEOELEMENT1D maxItem(true), minItem(true);
  
  pVectGlobalManip<GEOELEMENT1D,ELMAP> facesManipulator(commDev);
  facesManipulator.dataMinMax(tempEdges2d,minEdge,maxEdge);
  facesManipulator.dataMinMax(elements1d,minEl1d,maxEl1d);
  
  minItem = std::min(minEdge,minEl1d);
  maxItem = std::max(maxEdge,maxEl1d);  
  
  pVectComm<GEOELEMENT1D,ELMAP> edgesCommunicator(commDev);
  edgesCommunicator.vectorData(tempEdges2d,minItem,maxItem);
  edgesCommunicator.vectorData(elements1d,minItem,maxItem);
  
  //Ids matching -> elements2d creation
  typedef pair<GEOELEMENT1D,ELMAP>                    VALUETYPE;
  typedef typename map<GEOELEMENT1D,ELMAP>::iterator  ITERATOR1D;
  typedef pair<ITERATOR1D,bool>                       PAIR1D;
  
  map<GEOELEMENT1D,ELMAP>  elements1dSet;
  VALUETYPE                value1d;
  ITERATOR1D               iter;
  PAIR1D                   ret;
  
  for(UInt i=1; i <= elements1d.size(); ++i)
  {
    value1d.first  = elements1d.getItemL(i);
    value1d.second = elements1d.getRowMapL(i);
    
    ret = elements1dSet.insert(value1d);
    assert(ret.second == true); //The elements1d should be unique
  }
  
  elements1d.clear();
  UInt fLid = 1;
  
  for(UInt i=1; i <= tempEdges2d.size(); ++i)
  {
    edgeGeoItem = tempEdges2d.getItemL(i);
    edgeMapItem = tempEdges2d.getRowMapL(i);
    
    if(elements1dSet.count(edgeGeoItem) == 1)
    {     
      iter = elements1dSet.find(edgeGeoItem);
      
      edgeGeoItem.setGeoId(iter->first.getGeoId());
      edgeMapItem.setGid(iter->second.getGid());
      edgeMapItem.setLid(fLid);
      
      fLid++;
      
      elements1d.push_back(edgeGeoItem,edgeMapItem);
    }
  }
  
  //Communicate back to processors
  edgesCommunicator.vectorPid(elements1d);
  
  //Local indexing
  pGraphManip<GEOELEMENT1D,ELMAP,NODEMAP> local1dManip;
  local1dManip.setNormalIndexing(elements1d);
  
  //Push to local the elements2d
  elements1d.setColMap(grid2d->getNodes().getMapRef());
  elements1d.pushToLocal();
  
  //Create the grid2d
  grid1d = rcp(new mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>);
  grid1d->setNodes(grid2d->getNodes());
  grid1d->setElements(elements1d);
  grid1d->transferMap();
  grid1d->setMeshStandard(STDA);
  
  //Remove unused points
  meshDoctor1d<GEOSHAPE1D,ELMAP,NODEMAP> gridDoctor1d(commDev);
  gridDoctor1d.removeUnusedPoints(grid1d);
  
  //Grid2d checking
  mesh1dGlobalManip<GEOSHAPE1D,ELMAP,NODEMAP> gridManip1d(commDev);
  assert(gridManip1d.check(grid1d));
 
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Connecting1d-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting surface 1d" << " "; time(&start); cout << endl;}
  
  connectGrid1d = rcp(new connect1d<GEOSHAPE1D,ELMAP,NODEMAP>(commDev));
  connectGrid1d->setMesh1d(grid1d);
  connectGrid1d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}  
  
  
  //Boundary2d connecting--------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Boundary connecting 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d->setMesh1d(grid1d);
  connectGrid2d->buildBoundaryConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Flags------------------------------------------------------------------------------------------
  mesh2dCreated    = true;
  mesh1dCreated    = true;
  connect2dCreated = true;
  connect1dCreated = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
gmMesh_to_stdB(const string & meshfile, const bool & verbose)
{
  assert(staticAssert< GEOSHAPE::geoName == LINEARTRIANGLE >::returnValue);
  assert(commDevLoaded);
  
  //Allocations------------------------------------------------------------------------------------
  UInt pid      = commDev->rank();
  UInt numPids  = commDev->size();
  UInt printPid = 0;
  time_t start, end;
  
  //Allocations
  NODESVECT nodes;
  GRAPH2D   elements2d;
  GRAPH1D   elements1d;
  
  
  //Mesh loading-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Start mesh loading" << " "; time(&start); cout << endl;}
  loadMesh<GEOSHAPE,ELMAP,NODEMAP> meshLoader;
  meshLoader.setParam(pid,printPid,numPids);
  meshLoader.gmMesh(meshfile,nodes,elements2d,elements1d);
  
  elements1d.colIsLocal() = false;
  elements2d.colIsLocal() = false;
  
  commDev->barrier();
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //The grid3d-------------------------------------------------------------------------------------
  grid2d = Teuchos::rcp(new mesh2d<GEOSHAPE,ELMAP,NODEMAP>);
  grid2d->setNodes(nodes);
  grid2d->setElements(elements2d);
  grid2d->transferMap();
  grid2d->setMeshStandard(STDL);
  
  
  //Mesh partition---------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh partition" << " "; time(&start); cout << endl;}
  
  mesh2dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip2d(commDev);
  gridManip2d.meshPartition(grid2d);
  assert(gridManip2d.check(grid2d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh overlapping-------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh cheking" << " "; time(&start); cout << endl;}
  
  meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP> gridDoctor(commDev);
  gridDoctor.fixGeoIds(grid2d);
  assert(gridManip2d.check(grid2d));
  assert(gridDoctor.checkGeoIds(grid2d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Surface 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d = Teuchos::rcp(new connect2d<GEOSHAPE2D,ELMAP,NODEMAP>(commDev));
  connectGrid2d->setMesh2d(grid2d);
  connectGrid2d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Edges cheking---------------------------------------------------------------------------------- 
  if((pid == printPid) && verbose) {cout << "Check Edges" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(commDev);
  assert(checkEdges.check(grid2d->getEdges()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Distribute the elements2d----------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Distribute elements 1d" << " "; time(&start); cout << endl;}
    
  //Find local cadidates - the ones that has only one element connected
  GEOELEMENT1D  edgeGeoItem(true);
  ELMAP         edgeMapItem;
  
  GRAPH1D tempEdges2d = grid2d->getEdges();
  tempEdges2d.pushToGlobal();
  
  
  //Normal communication
  GEOELEMENT1D maxEdge(true), minEdge(true);
  GEOELEMENT1D maxEl1d(true), minEl1d(true);
  GEOELEMENT1D maxItem(true), minItem(true);
  
  pVectGlobalManip<GEOELEMENT1D,ELMAP> facesManipulator(commDev);
  facesManipulator.dataMinMax(tempEdges2d,minEdge,maxEdge);
  facesManipulator.dataMinMax(elements1d,minEl1d,maxEl1d);
  
  minItem = std::min(minEdge,minEl1d);
  maxItem = std::max(maxEdge,maxEl1d);  
  
  pVectComm<GEOELEMENT1D,ELMAP> edgesCommunicator(commDev);
  edgesCommunicator.vectorData(tempEdges2d,minItem,maxItem);
  edgesCommunicator.vectorData(elements1d,minItem,maxItem);
  
  //Ids matching -> elements2d creation
  typedef pair<GEOELEMENT1D,ELMAP>                    VALUETYPE;
  typedef typename map<GEOELEMENT1D,ELMAP>::iterator  ITERATOR1D;
  typedef pair<ITERATOR1D,bool>                       PAIR1D;
  
  map<GEOELEMENT1D,ELMAP>  elements1dSet;
  VALUETYPE                value1d;
  ITERATOR1D               iter;
  PAIR1D                   ret;
  
  for(UInt i=1; i <= elements1d.size(); ++i)
  {
    value1d.first  = elements1d.getItemL(i);
    value1d.second = elements1d.getRowMapL(i);
    
    ret = elements1dSet.insert(value1d);
    assert(ret.second == true); //The elements1d should be unique
  }
  
  elements1d.clear();
  UInt fLid = 1;
  
  for(UInt i=1; i <= tempEdges2d.size(); ++i)
  {
    edgeGeoItem = tempEdges2d.getItemL(i);
    edgeMapItem = tempEdges2d.getRowMapL(i);
    
    if(elements1dSet.count(edgeGeoItem) == 1)
    {     
      iter = elements1dSet.find(edgeGeoItem);
      
      edgeGeoItem.setGeoId(iter->first.getGeoId());
      edgeMapItem.setGid(iter->second.getGid());
      edgeMapItem.setLid(fLid);
      
      fLid++;
      
      elements1d.push_back(edgeGeoItem,edgeMapItem);
    }
  }
  
  //Communicate back to processors
  edgesCommunicator.vectorPid(elements1d);
  
  //Local indexing
  pGraphManip<GEOELEMENT1D,ELMAP,NODEMAP> local1dManip;
  local1dManip.setNormalIndexing(elements1d);
  
  //Push to local the elements2d
  elements1d.setColMap(grid2d->getNodes().getMapRef());
  elements1d.pushToLocal();
  
  //Create the grid2d
  grid1d = Teuchos::rcp(new mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>);
  grid1d->setNodes(grid2d->getNodes());
  grid1d->setElements(elements1d);
  grid1d->transferMap();
  grid1d->setMeshStandard(STDB);
  
  //Remove unused points
  meshDoctor1d<GEOSHAPE1D,ELMAP,NODEMAP> gridDoctor1d(commDev);
  gridDoctor1d.removeUnusedPoints(grid1d);
  
  //Grid2d checking
  mesh1dGlobalManip<GEOSHAPE1D,ELMAP,NODEMAP> gridManip1d(commDev);
  assert(gridManip1d.check(grid1d));
 
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Connecting1d-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting surface 1d" << " "; time(&start); cout << endl;}
  
  connectGrid1d = Teuchos::rcp(new connect1d<GEOSHAPE1D,ELMAP,NODEMAP>(commDev));
  connectGrid1d->setMesh1d(grid1d);
  connectGrid1d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}  
  
  
  //Boundary2d connecting--------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Boundary connecting 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d->setMesh1d(grid1d);
  connectGrid2d->buildBoundaryConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Flags------------------------------------------------------------------------------------------
  mesh2dCreated    = true;
  mesh1dCreated    = true;
  connect2dCreated = true;
  connect1dCreated = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
femap_to_stdA(const string & meshfile, const string & colorfile, const bool & verbose)
{
  assert(commDevLoaded);
  assert(staticAssert< GEOSHAPE::geoName == LINEARQUAD >::returnValue);
  
  //Allocations------------------------------------------------------------------------------------
  UInt pid      = commDev->rank();
  UInt numPids  = commDev->size();
  UInt printPid = 0;
  time_t start, end;
  
  //Load coloring----------------------------------------------------------------------------------
  load readColor;
  UInt numColors, color;
  sVect<UInt> powerColor;
  
  //File check
  ifstream file(colorfile.c_str());
  
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Reading file  
  readColor.leggiriga(file);
  readColor.assegnavalore(0,numColors);
  
  powerColor.reserve(numColors);
  
  for(UInt i=1; i <= numColors; ++i)
  {
    readColor.leggiriga(file);
    readColor.assegnavalore(0,color);
    
    powerColor.push_back(color);
  }
  
  //Closing file
  file.close();
  
  //Testing color file
  typedef set<UInt>::iterator ITERATOR;
 
  UInt maxColor=1, minColor=1;
  set<UInt> iterList;
  pair<ITERATOR,bool> pairColor;
  
  for(UInt i=1; i <= powerColor.size(); ++i)
  {
    pairColor = iterList.insert(powerColor(i));
    assert(pairColor.second);
    
    maxColor = max(maxColor, powerColor(i));
    minColor = min(minColor, powerColor(i));
  }
  
  assert(maxColor == numColors);
  assert(minColor == 1);
  
  
  //Mesh loading-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Start mesh loading" << " "; time(&start); cout << endl;}
  
  typedef pVect<UInt,NODEMAP> COLORVECT;
  
  NODESVECT nodes;
  GRAPH2D   elements2d;
  COLORVECT nodesColor;
  
  loadMesh<GEOSHAPE,ELMAP,NODEMAP> meshLoader;
  meshLoader.setParam(pid,printPid,numPids);
  meshLoader.nastran(meshfile,nodes,elements2d,nodesColor);
  
  elements2d.colIsLocal() = false;
  
  commDev->barrier();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //The grid2d-------------------------------------------------------------------------------------
  grid2d = rcp(new mesh2d<GEOSHAPE,ELMAP,NODEMAP>);
  grid2d->setNodes(nodes);
  grid2d->setElements(elements2d);
  grid2d->transferMap();
  grid2d->setMeshStandard(STDL);
  
  
  //Mesh partition---------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh partition" << " "; time(&start); cout << endl;}
  
  mesh2dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip2d(commDev);
  gridManip2d.meshPartition(grid2d);
  assert(gridManip2d.check(grid2d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh overlapping-------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh overlap" << " "; time(&start); cout << endl;}
  
  meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP> gridDoctor(commDev);
  gridManip2d.meshOverlapping(grid2d);
  assert(gridManip2d.check(grid2d));
  assert(gridDoctor.checkGeoIds(grid2d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Surface 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d = rcp(new connect2d<GEOSHAPE2D,ELMAP,NODEMAP>(commDev));
  connectGrid2d->setMesh2d(grid2d);
  connectGrid2d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Edges cheking----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Check Edges" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(commDev);
  assert(checkEdges.check(grid2d->getEdges()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Cheking colors---------------------------------------------------------------------------------
  UInt locMax = 0;
  
  for(UInt i=1; i <= nodesColor.size(); ++i)
  { locMax = max(locMax,nodesColor(i)); }
  
  sVect<UInt> vectMax(numPids);
  all_gather(*commDev,locMax,vectMax);
  
  for(UInt i=1; i <= numPids; ++i)
  { locMax = max(locMax,vectMax(i));}
  
  assert(locMax == numColors);
  
  
  //Generate 1d elements---------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Generate elements 1d" << " "; time(&start); cout << endl;}  
  
  //Remap of the nodesColor
  pVectGlobalManip<UInt,NODEMAP> nodesColorComm(commDev);
  
  pMap<NODEMAP> newColorMap = grid2d->getNodes().getMapRef();
  nodesColorComm.changeMap(nodesColor,newColorMap);
  
  //Global edges generation
  pGraph<pGraphItem,ELMAP,ELMAP> elementsToEdges = connectGrid2d->getElementToEdge();
  pGraph<pGraphItem,ELMAP,ELMAP> globalEdgesToElments;
  
  pGraphGlobalManip<pGraphItem,ELMAP,ELMAP>  tempGraphManip(commDev);
  tempGraphManip.inversion(elementsToEdges,globalEdgesToElments);
  
  //Get one-connected faces
  typedef typename MESH2D::GRAPH1D GRAPH1D;
  GRAPH1D      elements1d;
  GEOELEMENT1D el1d;
  ELMAP        el1dMap;
  UInt minPower, nodeLid;
  UInt tot = 1;
  
  assert(grid2d->getEdges().size() == globalEdgesToElments.size());
  
  for(UInt i=1; i <= grid2d->getEdges().size(); ++i)
  {
    assert((globalEdgesToElments(i).size() == 1) || (globalEdgesToElments(i).size() == 2));
    
    if(globalEdgesToElments(i).size() == 1)
    {
      el1d    = grid2d->getEdges().getItemL(i);
      el1dMap = grid2d->getEdges().getRowMapL(i);
      
      assert(el1d.size() == 2);
      
      //GeoId determination
      nodeLid  = el1d.getCid(1);      
      color    = nodesColor(nodeLid);
      minColor = color;
      minPower = powerColor(color);
      
      for(UInt j=1; j <= el1d.size(); ++j)
      {
	nodeLid  = el1d.getCid(j);
        color    = nodesColor(nodeLid);
	
	if(powerColor(color) < minPower)
	{
	  minColor = color;
          minPower = powerColor(color);
	}
      }
      
      el1d.setGeoId(minColor);
      
      //Set local Id
      el1dMap.setPid(pid);
      el1dMap.setLid(tot);
      tot++;
      
      //Add the 2d element
      elements1d.push_back(el1dMap,el1d);
    }
  }
  
  elements1d.setColMap(grid2d->getNodes().getMapRef());
  
  //Fixing the ownership and global nombering of the 2d elements
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> BelementsManip(commDev);
  BelementsManip.buildGlobalNumbering(elements1d);
  
  //Create the grid2d
  grid1d = rcp(new mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>);
  grid1d->setNodes(grid2d->getNodes());
  grid1d->setElements(elements1d);
  grid1d->transferMap();
  grid1d->setMeshStandard(STDA);
  
  //Remove unused points
  meshDoctor1d<GEOSHAPE1D,ELMAP,NODEMAP> gridDoctor1d(commDev);
  gridDoctor1d.removeUnusedPoints(grid1d);
 
  //Grid2d checking
  mesh1dGlobalManip<GEOSHAPE1D,ELMAP,NODEMAP> gridManip1d(commDev);
  assert(gridManip1d.check(grid1d));
 
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Connecting1d-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting line 1d" << " "; time(&start); cout << endl;}
  
  connectGrid1d = rcp(new connect1d<GEOSHAPE1D,ELMAP,NODEMAP>(commDev));
  connectGrid1d->setMesh1d(grid1d);
  connectGrid1d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}  
  
  
  //Boundary2d connecting--------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Boundary connecting 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d->setMesh1d(grid1d);
  connectGrid2d->buildBoundaryConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Flags------------------------------------------------------------------------------------------
  mesh2dCreated    = true;
  mesh1dCreated    = true;
  connect2dCreated = true;
  connect1dCreated = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
femap_to_stdB(const string & meshfile, const string & colorfile, const bool & verbose)
{
  assert(commDevLoaded);
  assert(staticAssert< GEOSHAPE::geoName == LINEARQUAD >::returnValue);
  
  //Allocations------------------------------------------------------------------------------------
  UInt pid      = commDev->rank();
  UInt numPids  = commDev->size();
  UInt printPid = 0;
  time_t start, end;
  
  //Load coloring----------------------------------------------------------------------------------
  load readColor;
  UInt numColors, color;
  sVect<UInt> powerColor;
  
  //File check
  ifstream file(colorfile.c_str());
  
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Reading file  
  readColor.leggiriga(file);
  readColor.assegnavalore(0,numColors);
  
  powerColor.reserve(numColors);
  
  for(UInt i=1; i <= numColors; ++i)
  {
    readColor.leggiriga(file);
    readColor.assegnavalore(0,color);
    
    powerColor.push_back(color);
  }
  
  //Closing file
  file.close();
  
  //Testing color file
  typedef set<UInt>::iterator ITERATOR;
 
  UInt maxColor=1, minColor=1;
  set<UInt> iterList;
  pair<ITERATOR,bool> pairColor;
  
  for(UInt i=1; i <= powerColor.size(); ++i)
  {
    pairColor = iterList.insert(powerColor(i));
    assert(pairColor.second);
    
    maxColor = max(maxColor, powerColor(i));
    minColor = min(minColor, powerColor(i));
  }
  
  assert(maxColor == numColors);
  assert(minColor == 1);
  
  
  //Mesh loading-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Start mesh loading" << " "; time(&start); cout << endl;}
  
  typedef pVect<UInt,NODEMAP> COLORVECT;
  
  NODESVECT nodes;
  GRAPH2D   elements2d;
  COLORVECT nodesColor;
  
  loadMesh<GEOSHAPE,ELMAP,NODEMAP> meshLoader;
  meshLoader.setParam(pid,printPid,numPids);
  meshLoader.nastran(meshfile,nodes,elements2d,nodesColor);
  
  elements2d.colIsLocal() = false;
  
  commDev->barrier();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //The grid2d-------------------------------------------------------------------------------------
  grid2d = Teuchos::rcp(new mesh2d<GEOSHAPE,ELMAP,NODEMAP>);
  grid2d->setNodes(nodes);
  grid2d->setElements(elements2d);
  grid2d->transferMap();
  grid2d->setMeshStandard(STDL);
  
  
  //Mesh partition---------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh partition" << " "; time(&start); cout << endl;}
  
  mesh2dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip2d(commDev);
  gridManip2d.meshPartition(grid2d);
  assert(gridManip2d.check(grid2d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Surface 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d = Teuchos::rcp(new connect2d<GEOSHAPE2D,ELMAP,NODEMAP>(commDev));
  connectGrid2d->setMesh2d(grid2d);
  connectGrid2d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Edges cheking----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Check Edges" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(commDev);
  assert(checkEdges.check(grid2d->getEdges()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Cheking colors---------------------------------------------------------------------------------
  UInt locMax = 0;
  
  for(UInt i=1; i <= nodesColor.size(); ++i)
  { locMax = max(locMax,nodesColor(i)); }
  
  sVect<UInt> vectMax(numPids);
  all_gather(*commDev,locMax,vectMax);
  
  for(UInt i=1; i <= numPids; ++i)
  { locMax = max(locMax,vectMax(i));}
  
  assert(locMax == numColors);
  
  
  //Generate 1d elements---------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Generate elements 1d" << " "; time(&start); cout << endl;}  
  
  //Remap of the nodesColor
  pVectGlobalManip<UInt,NODEMAP> nodesColorComm(commDev);
  
  pMap<NODEMAP> newColorMap = grid2d->getNodes().getMapRef();
  nodesColorComm.changeMap(nodesColor,newColorMap);
  
  //Global edges generation
  pGraph<pGraphItem,ELMAP,ELMAP> elementsToEdges = connectGrid2d->getElementToEdge();
  pGraph<pGraphItem,ELMAP,ELMAP> globalEdgesToElments;
  
  pGraphGlobalManip<pGraphItem,ELMAP,ELMAP>  tempGraphManip(commDev);
  tempGraphManip.inversion(elementsToEdges,globalEdgesToElments);
  
  //Get one-connected faces
  typedef typename MESH2D::GRAPH1D GRAPH1D;
  GRAPH1D      elements1d;
  GEOELEMENT1D el1d;
  ELMAP        el1dMap;
  UInt minPower, nodeLid;
  UInt tot = 1;
  
  assert(grid2d->getEdges().size() == globalEdgesToElments.size());
  
  for(UInt i=1; i <= grid2d->getEdges().size(); ++i)
  {
    assert((globalEdgesToElments(i).size() == 1) || (globalEdgesToElments(i).size() == 2));
    
    if(globalEdgesToElments(i).size() == 1)
    {
      el1d    = grid2d->getEdges().getItemL(i);
      el1dMap = grid2d->getEdges().getRowMapL(i);
      
      assert(el1d.size() == 2);
      
      //GeoId determination
      nodeLid  = el1d.getCid(1);      
      color    = nodesColor(nodeLid);
      minColor = color;
      minPower = powerColor(color);
      
      for(UInt j=1; j <= el1d.size(); ++j)
      {
	nodeLid  = el1d.getCid(j);
        color    = nodesColor(nodeLid);
	
	if(powerColor(color) < minPower)
	{
	  minColor = color;
          minPower = powerColor(color);
	}
      }
      
      el1d.setGeoId(minColor);
      
      //Set local Id
      el1dMap.setPid(pid);
      el1dMap.setLid(tot);
      tot++;
      
      //Add the 2d element
      elements1d.push_back(el1dMap,el1d);
    }
  }
  
  elements1d.setColMap(grid2d->getNodes().getMapRef());
  
  //Fixing the ownership and global nombering of the 2d elements
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> BelementsManip(commDev);
  BelementsManip.buildGlobalNumbering(elements1d);
  
  //Create the grid2d
  grid1d = Teuchos::rcp(new mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>);
  grid1d->setNodes(grid2d->getNodes());
  grid1d->setElements(elements1d);
  grid1d->transferMap();
  grid1d->setMeshStandard(STDB);
  
  //Remove unused points
  meshDoctor1d<GEOSHAPE1D,ELMAP,NODEMAP> gridDoctor1d(commDev);
  gridDoctor1d.removeUnusedPoints(grid1d);
 
  //Grid2d checking
  mesh1dGlobalManip<GEOSHAPE1D,ELMAP,NODEMAP> gridManip1d(commDev);
  assert(gridManip1d.check(grid1d));
 
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Connecting1d-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting line 1d" << " "; time(&start); cout << endl;}
  
  connectGrid1d = Teuchos::rcp(new connect1d<GEOSHAPE1D,ELMAP,NODEMAP>(commDev));
  connectGrid1d->setMesh1d(grid1d);
  connectGrid1d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}  
  
  
  //Boundary2d connecting--------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Boundary connecting 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d->setMesh1d(grid1d);
  connectGrid2d->buildBoundaryConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Flags------------------------------------------------------------------------------------------
  mesh2dCreated    = true;
  mesh1dCreated    = true;
  connect2dCreated = true;
  connect1dCreated = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
neutral_to_stdA(const string & meshfile, const string & colorfile, const bool & verbose)
{
  assert(commDevLoaded);
  assert(staticAssert< GEOSHAPE::geoName == LINEARQUAD >::returnValue);
  
  //Allocations------------------------------------------------------------------------------------
  UInt pid      = commDev->rank();
  UInt numPids  = commDev->size();
  UInt printPid = 0;
  time_t start, end;
  
  //Load coloring----------------------------------------------------------------------------------
  load readColor;
  UInt numColors, color;
  sVect<UInt> powerColor;
  
  //File check
  ifstream file(colorfile.c_str());
  
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Reading file  
  readColor.leggiriga(file);
  readColor.assegnavalore(0,numColors);
  
  powerColor.reserve(numColors);
  
  for(UInt i=1; i <= numColors; ++i)
  {
    readColor.leggiriga(file);
    readColor.assegnavalore(0,color);
    
    powerColor.push_back(color);
  }
  
  //Closing file
  file.close();
  
  //Testing color file
  typedef set<UInt>::iterator ITERATOR;
 
  UInt maxColor=1, minColor=1;
  set<UInt> iterList;
  pair<ITERATOR,bool> pairColor;
  
  for(UInt i=1; i <= powerColor.size(); ++i)
  {
    pairColor = iterList.insert(powerColor(i));
    assert(pairColor.second);
    
    maxColor = max(maxColor, powerColor(i));
    minColor = min(minColor, powerColor(i));
  }
  
  assert(maxColor == numColors);
  assert(minColor == 1);
  
  
  //Mesh loading-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Start mesh loading" << " "; time(&start); cout << endl;}
  
  typedef pVect<UInt,NODEMAP> COLORVECT;
  
  NODESVECT nodes;
  GRAPH2D   elements2d;
  COLORVECT nodesColor;
  
  loadMesh<GEOSHAPE,ELMAP,NODEMAP> meshLoader;
  meshLoader.setParam(pid,printPid,numPids);
  meshLoader.neutral(meshfile,nodes,elements2d,nodesColor);
  
  elements2d.colIsLocal() = false;
  
  commDev->barrier();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //The grid2d-------------------------------------------------------------------------------------
  grid2d = rcp(new mesh2d<GEOSHAPE,ELMAP,NODEMAP>);
  grid2d->setNodes(nodes);
  grid2d->setElements(elements2d);
  grid2d->transferMap();
  grid2d->setMeshStandard(STDL);
  
  
  //Mesh partition---------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh partition" << " "; time(&start); cout << endl;}
  
  mesh2dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip2d(commDev);
  gridManip2d.meshPartition(grid2d);
  assert(gridManip2d.check(grid2d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh overlapping-------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh overlap" << " "; time(&start); cout << endl;}
  
  meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP> gridDoctor(commDev);
  gridManip2d.meshOverlapping(grid2d);
  assert(gridManip2d.check(grid2d));
  assert(gridDoctor.checkGeoIds(grid2d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Surface 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d = rcp(new connect2d<GEOSHAPE2D,ELMAP,NODEMAP>(commDev));
  connectGrid2d->setMesh2d(grid2d);
  connectGrid2d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Edges cheking----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Check Edges" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(commDev);
  assert(checkEdges.check(grid2d->getEdges()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Cheking colors---------------------------------------------------------------------------------
  UInt locMax = 0;
  
  for(UInt i=1; i <= nodesColor.size(); ++i)
  { locMax = max(locMax,nodesColor(i)); }
  
  sVect<UInt> vectMax(numPids);
  all_gather(*commDev,locMax,vectMax);
  
  for(UInt i=1; i <= numPids; ++i)
  { locMax = max(locMax,vectMax(i));}
  
  assert(locMax == numColors);
  
  
  //Generate 1d elements---------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Generate elements 1d" << " "; time(&start); cout << endl;}  
  
  //Remap of the nodesColor
  pVectGlobalManip<UInt,NODEMAP> nodesColorComm(commDev);
  
  pMap<NODEMAP> newColorMap = grid2d->getNodes().getMapRef();
  nodesColorComm.changeMap(nodesColor,newColorMap);
  
  //Global edges generation
  pGraph<pGraphItem,ELMAP,ELMAP> elementsToEdges = connectGrid2d->getElementToEdge();
  pGraph<pGraphItem,ELMAP,ELMAP> globalEdgesToElments;
  
  pGraphGlobalManip<pGraphItem,ELMAP,ELMAP>  tempGraphManip(commDev);
  tempGraphManip.inversion(elementsToEdges,globalEdgesToElments);
  
  //Get one-connected faces
  typedef typename MESH2D::GRAPH1D GRAPH1D;
  GRAPH1D      elements1d;
  GEOELEMENT1D el1d;
  ELMAP        el1dMap;
  UInt minPower, nodeLid;
  UInt tot = 1;
  
  assert(grid2d->getEdges().size() == globalEdgesToElments.size());
  
  for(UInt i=1; i <= grid2d->getEdges().size(); ++i)
  {
    assert((globalEdgesToElments(i).size() == 1) || (globalEdgesToElments(i).size() == 2));
    
    if(globalEdgesToElments(i).size() == 1)
    {
      el1d    = grid2d->getEdges().getItemL(i);
      el1dMap = grid2d->getEdges().getRowMapL(i);
      
      assert(el1d.size() == 2);
      
      //GeoId determination
      nodeLid  = el1d.getCid(1);      
      color    = nodesColor(nodeLid);
      minColor = color;
      minPower = powerColor(color);
      
      for(UInt j=1; j <= el1d.size(); ++j)
      {
	nodeLid  = el1d.getCid(j);
        color    = nodesColor(nodeLid);
	
	if(powerColor(color) < minPower)
	{
	  minColor = color;
          minPower = powerColor(color);
	}
      }
      
      el1d.setGeoId(minColor);
      
      //Set local Id
      el1dMap.setPid(pid);
      el1dMap.setLid(tot);
      tot++;
      
      //Add the 2d element
      elements1d.push_back(el1dMap,el1d);
    }
  }
  
  elements1d.setColMap(grid2d->getNodes().getMapRef());
  
  //Fixing the ownership and global nombering of the 2d elements
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> BelementsManip(commDev);
  BelementsManip.buildGlobalNumbering(elements1d);
  
  //Create the grid2d
  grid1d = rcp(new mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>);
  grid1d->setNodes(grid2d->getNodes());
  grid1d->setElements(elements1d);
  grid1d->transferMap();
  grid1d->setMeshStandard(STDA);
  
  //Remove unused points
  meshDoctor1d<GEOSHAPE1D,ELMAP,NODEMAP> gridDoctor1d(commDev);
  gridDoctor1d.removeUnusedPoints(grid1d);
 
  //Grid2d checking
  mesh1dGlobalManip<GEOSHAPE1D,ELMAP,NODEMAP> gridManip1d(commDev);
  assert(gridManip1d.check(grid1d));
 
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Connecting1d-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting line 1d" << " "; time(&start); cout << endl;}
  
  connectGrid1d = rcp(new connect1d<GEOSHAPE1D,ELMAP,NODEMAP>(commDev));
  connectGrid1d->setMesh1d(grid1d);
  connectGrid1d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}  
  
  
  //Boundary2d connecting--------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Boundary connecting 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d->setMesh1d(grid1d);
  connectGrid2d->buildBoundaryConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Flags------------------------------------------------------------------------------------------
  mesh2dCreated    = true;
  mesh1dCreated    = true;
  connect2dCreated = true;
  connect1dCreated = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit2d<GEOSHAPE,ELMAP,NODEMAP>::
neutral_to_stdB(const string & meshfile, const string & colorfile, const bool & verbose)
{
  assert(commDevLoaded);
  assert(staticAssert< GEOSHAPE::geoName == LINEARQUAD >::returnValue);
  
  //Allocations------------------------------------------------------------------------------------
  UInt pid      = commDev->rank();
  UInt numPids  = commDev->size();
  UInt printPid = 0;
  time_t start, end;
  
  //Load coloring----------------------------------------------------------------------------------
  load readColor;
  UInt numColors, color;
  sVect<UInt> powerColor;
  
  //File check
  ifstream file(colorfile.c_str());
  
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Reading file  
  readColor.leggiriga(file);
  readColor.assegnavalore(0,numColors);
  
  powerColor.reserve(numColors);
  
  for(UInt i=1; i <= numColors; ++i)
  {
    readColor.leggiriga(file);
    readColor.assegnavalore(0,color);
    
    powerColor.push_back(color);
  }
  
  //Closing file
  file.close();
  
  //Testing color file
  typedef set<UInt>::iterator ITERATOR;
 
  UInt maxColor=1, minColor=1;
  set<UInt> iterList;
  pair<ITERATOR,bool> pairColor;
  
  for(UInt i=1; i <= powerColor.size(); ++i)
  {
    pairColor = iterList.insert(powerColor(i));
    assert(pairColor.second);
    
    maxColor = max(maxColor, powerColor(i));
    minColor = min(minColor, powerColor(i));
  }
  
  assert(maxColor == numColors);
  assert(minColor == 1);
  
  
  //Mesh loading-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Start mesh loading" << " "; time(&start); cout << endl;}
  
  typedef pVect<UInt,NODEMAP> COLORVECT;
  
  NODESVECT nodes;
  GRAPH2D   elements2d;
  COLORVECT nodesColor;
  
  loadMesh<GEOSHAPE,ELMAP,NODEMAP> meshLoader;
  meshLoader.setParam(pid,printPid,numPids);
  meshLoader.neutral(meshfile,nodes,elements2d,nodesColor);
  
  elements2d.colIsLocal() = false;
  
  commDev->barrier();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //The grid2d-------------------------------------------------------------------------------------
  grid2d = rcp(new mesh2d<GEOSHAPE,ELMAP,NODEMAP>);
  grid2d->setNodes(nodes);
  grid2d->setElements(elements2d);
  grid2d->transferMap();
  grid2d->setMeshStandard(STDL);
  
  
  //Mesh partition---------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh partition" << " "; time(&start); cout << endl;}
  
  mesh2dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip2d(commDev);
  gridManip2d.meshPartition(grid2d);
  assert(gridManip2d.check(grid2d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Surface 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d = rcp(new connect2d<GEOSHAPE2D,ELMAP,NODEMAP>(commDev));
  connectGrid2d->setMesh2d(grid2d);
  connectGrid2d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Edges cheking----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Check Edges" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(commDev);
  assert(checkEdges.check(grid2d->getEdges()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Cheking colors---------------------------------------------------------------------------------
  UInt locMax = 0;
  
  for(UInt i=1; i <= nodesColor.size(); ++i)
  { locMax = max(locMax,nodesColor(i)); }
  
  sVect<UInt> vectMax(numPids);
  all_gather(*commDev,locMax,vectMax);
  
  for(UInt i=1; i <= numPids; ++i)
  { locMax = max(locMax,vectMax(i));}
  
  assert(locMax == numColors);
  
  
  //Generate 1d elements---------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Generate elements 1d" << " "; time(&start); cout << endl;}  
  
  //Remap of the nodesColor
  pVectGlobalManip<UInt,NODEMAP> nodesColorComm(commDev);
  
  pMap<NODEMAP> newColorMap = grid2d->getNodes().getMapRef();
  nodesColorComm.changeMap(nodesColor,newColorMap);
  
  //Global edges generation
  pGraph<pGraphItem,ELMAP,ELMAP> elementsToEdges = connectGrid2d->getElementToEdge();
  pGraph<pGraphItem,ELMAP,ELMAP> globalEdgesToElments;
  
  pGraphGlobalManip<pGraphItem,ELMAP,ELMAP>  tempGraphManip(commDev);
  tempGraphManip.inversion(elementsToEdges,globalEdgesToElments);
  
  //Get one-connected faces
  typedef typename MESH2D::GRAPH1D GRAPH1D;
  GRAPH1D      elements1d;
  GEOELEMENT1D el1d;
  ELMAP        el1dMap;
  UInt minPower, nodeLid;
  UInt tot = 1;
  
  assert(grid2d->getEdges().size() == globalEdgesToElments.size());
  
  for(UInt i=1; i <= grid2d->getEdges().size(); ++i)
  {
    assert((globalEdgesToElments(i).size() == 1) || (globalEdgesToElments(i).size() == 2));
    
    if(globalEdgesToElments(i).size() == 1)
    {
      el1d    = grid2d->getEdges().getItemL(i);
      el1dMap = grid2d->getEdges().getRowMapL(i);
      
      assert(el1d.size() == 2);
      
      //GeoId determination
      nodeLid  = el1d.getCid(1);      
      color    = nodesColor(nodeLid);
      minColor = color;
      minPower = powerColor(color);
      
      for(UInt j=1; j <= el1d.size(); ++j)
      {
	nodeLid  = el1d.getCid(j);
        color    = nodesColor(nodeLid);
	
	if(powerColor(color) < minPower)
	{
	  minColor = color;
          minPower = powerColor(color);
	}
      }
      
      el1d.setGeoId(minColor);
      
      //Set local Id
      el1dMap.setPid(pid);
      el1dMap.setLid(tot);
      tot++;
      
      //Add the 2d element
      elements1d.push_back(el1dMap,el1d);
    }
  }
  
  elements1d.setColMap(grid2d->getNodes().getMapRef());
  
  //Fixing the ownership and global nombering of the 2d elements
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> BelementsManip(commDev);
  BelementsManip.buildGlobalNumbering(elements1d);
  
  //Create the grid2d
  grid1d = rcp(new mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>);
  grid1d->setNodes(grid2d->getNodes());
  grid1d->setElements(elements1d);
  grid1d->transferMap();
  grid1d->setMeshStandard(STDB);
  
  //Remove unused points
  meshDoctor1d<GEOSHAPE1D,ELMAP,NODEMAP> gridDoctor1d(commDev);
  gridDoctor1d.removeUnusedPoints(grid1d);
 
  //Grid2d checking
  mesh1dGlobalManip<GEOSHAPE1D,ELMAP,NODEMAP> gridManip1d(commDev);
  assert(gridManip1d.check(grid1d));
 
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Connecting1d-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting line 1d" << " "; time(&start); cout << endl;}
  
  connectGrid1d = rcp(new connect1d<GEOSHAPE1D,ELMAP,NODEMAP>(commDev));
  connectGrid1d->setMesh1d(grid1d);
  connectGrid1d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}  
  
  
  //Boundary2d connecting--------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Boundary connecting 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d->setMesh1d(grid1d);
  connectGrid2d->buildBoundaryConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Flags------------------------------------------------------------------------------------------
  mesh2dCreated    = true;
  mesh1dCreated    = true;
  connect2dCreated = true;
  connect1dCreated = true;
}


#endif
