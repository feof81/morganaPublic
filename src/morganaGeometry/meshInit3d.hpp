/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MESHINIT3D_HPP
#define MESHINIT3D_HPP

#include <map>
#include <set>

#include "time.h"
#include <boost/iterator/iterator_concepts.hpp>

#include "pMapItem.h"
#include "pMapGlobalManip.h"

#include "pGraph.hpp"
#include "pGraphItem.h"
#include "pGraphItemOriented.h"
#include "pGraphItemSubLoc.h"
#include "pGraphGlobalManip.hpp"

#include "loadMesh.hpp"

#include "connect2d.hpp"
#include "connect3d.hpp"

#include "pGraphManip.hpp"

#include "mesh2dGlobalManip.hpp"
#include "mesh3dGlobalManip.hpp"
#include "meshDoctor2d.hpp"
#include "meshDoctor3d.hpp"


/*! Mesh initialization 3d */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class meshInit3d
{
    /*! @name Typedefs */ //@{
  public:
    typedef GEOSHAPE                             GEOSHAPE3D;
    typedef typename   GEOSHAPE::GEOBSHAPE       GEOSHAPE2D;
    typedef typename GEOSHAPE2D::GEOBSHAPE       GEOSHAPE1D;
    
    typedef geoElement<GEOSHAPE3D>               GEOELEMENT3D;
    typedef geoElement<GEOSHAPE2D>               GEOELEMENT2D;
    typedef geoElement<GEOSHAPE1D>               GEOELEMENT1D;
    
    typedef mesh3d<GEOSHAPE3D,ELMAP,NODEMAP>     MESH3D;
    typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>     MESH2D;
    
    typedef connect3d<GEOSHAPE3D,ELMAP,NODEMAP>  CONNECT3D;
    typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>  CONNECT2D;
    
    typedef typename MESH3D::NODESVECT           NODESVECT;
    typedef typename MESH3D::GRAPH3D             GRAPH3D;
    typedef typename MESH3D::GRAPH2D             GRAPH2D;
    //@}
    
    /*! @name Control flags */ //@{
  public:
    bool commDevLoaded;
    bool mesh3dCreated, mesh2dCreated;
    bool connect3dCreated, connect2dCreated;
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<communicator>  commDev;
    Teuchos::RCP<MESH3D>        grid3d;
    Teuchos::RCP<MESH2D>        grid2d;
    Teuchos::RCP<CONNECT3D>     connectGrid3d;
    Teuchos::RCP<CONNECT2D>     connectGrid2d;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    meshInit3d();
    meshInit3d(const Teuchos::RCP<communicator> & CommDev);
    meshInit3d(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Get functions */ //@{
  public:
    Teuchos::RCP<MESH3D>    getGrid3d();
    Teuchos::RCP<MESH2D>    getGrid2d();
    Teuchos::RCP<CONNECT3D> getConnectGrid3d();
    Teuchos::RCP<CONNECT2D> getConnectGrid2d();
    
    Teuchos::RCP<MESH3D>    getCopyGrid3d();
    Teuchos::RCP<MESH2D>    getCopyGrid2d();
    Teuchos::RCP<CONNECT3D> getCopyConnectGrid3d();
    Teuchos::RCP<CONNECT2D> getCopyConnectGrid2d();
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
    
    /*! @name Mesh reInit */ //@{
  public:
    void reinit_stdA(const Teuchos::RCP<MESH2D> & inGrid2d,
                     const Teuchos::RCP<MESH3D> & inGrid3d,
                     const set<UInt>            & activeGeoIds,
                     const bool                 & verbose = true);
    
    void reinit_stdB(const Teuchos::RCP<MESH2D> & inGrid2d,
                     const Teuchos::RCP<MESH3D> & inGrid3d,
                     const set<UInt>            & activeGeoIds,
                     const bool                 & verbose = true);
    //@}
};


//_______________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
meshInit3d()
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded     = false;
  mesh3dCreated     = false;
  mesh2dCreated     = false;
  connect3dCreated  = false;
  connect2dCreated  = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
meshInit3d(const Teuchos::RCP<communicator> & CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDev = CommDev;
  
  commDevLoaded     = true;
  mesh3dCreated     = false;
  mesh2dCreated     = false;
  connect3dCreated  = false;
  connect2dCreated  = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
meshInit3d(communicator & CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDev = Teuchos::rcpFromRef(CommDev);
  
  commDevLoaded     = true;
  mesh3dCreated     = false;
  mesh2dCreated     = false;
  connect3dCreated  = false;
  connect2dCreated  = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}



//_______________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::MESH3D>
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
getGrid3d()
{
  assert(mesh3dCreated);
  return(grid3d);
}
 
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::MESH2D>
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
getGrid2d()
{
  assert(mesh2dCreated);
  return(grid2d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::CONNECT3D>
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
getConnectGrid3d()
{
  assert(connect3dCreated);
  return(connectGrid3d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::CONNECT2D>
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
getConnectGrid2d()
{
  assert(connect2dCreated);
  return(connectGrid2d);
}



template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::MESH3D>
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
getCopyGrid3d()
{
  assert(mesh3dCreated);
  Teuchos::RCP<MESH3D>  newGrid3d(new MESH3D(*grid3d));
  return(newGrid3d);
}
 
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::MESH2D>
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
getCopyGrid2d()
{
  assert(mesh2dCreated);
  Teuchos::RCP<MESH2D>  newGrid2d(new MESH2D(*grid2d));
  return(newGrid2d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::CONNECT3D>
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
getCopyConnectGrid3d()
{
  assert(connect3dCreated);
  Teuchos::RCP<CONNECT3D>  newConnectGrid3d(new CONNECT3D(*connectGrid3d));
  return(newConnectGrid3d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Teuchos::RCP<typename meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::CONNECT2D>
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
getCopyConnectGrid2d()
{
  assert(connect2dCreated);
  Teuchos::RCP<CONNECT2D>  newConnectGrid2d(new CONNECT2D(*connectGrid2d));
  return(newConnectGrid2d);
}


//_______________________________________________________________________________________________________
// LOADING FUNCTIONS
//-------------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
gmMesh_to_stdA(const string & meshfile, const bool & verbose)
{
  assert(staticAssert< GEOSHAPE::geoName == LINEARTETRA >::returnValue);
  assert(commDevLoaded);
  
  //Allocations------------------------------------------------------------------------------------
  UInt pid      = commDev->rank();
  UInt numPids  = commDev->size();
  UInt printPid = 0;
  time_t start, end;
  
  //Allocations
  NODESVECT nodes;
  GRAPH3D   elements3d;
  GRAPH2D   elements2d;
  
  //Mesh loading-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Start mesh loading" << " "; time(&start); cout << endl;}
  loadMesh<GEOSHAPE,ELMAP,NODEMAP> meshLoader;
  meshLoader.setParam(pid,printPid,numPids);
  meshLoader.gmMeshParallel(meshfile,commDev,nodes,elements3d,elements2d);
  
  elements2d.colIsLocal() = false;
  elements3d.colIsLocal() = false;
  
  commDev->barrier();
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //The grid3d-------------------------------------------------------------------------------------
  grid3d = Teuchos::rcp(new mesh3d<GEOSHAPE,ELMAP,NODEMAP>);
  grid3d->setNodes(nodes);
  grid3d->setElements(elements3d);
  grid3d->transferMap();
  grid3d->setMeshStandard(STDL);
  
  //Mesh partition---------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh partition" << " "; time(&start); cout << endl;}
  
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip3d(commDev);
  gridManip3d.meshPartition(grid3d);
  assert(gridManip3d.check(grid3d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Mesh overlapping-------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh overlap" << " "; time(&start); cout << endl;}
  
  gridManip3d.meshOverlapping(grid3d);
  
  meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP> gridDoctor(commDev);
  gridDoctor.fixGeoIds(grid3d);
  assert(gridManip3d.check(grid3d));
  assert(gridDoctor.checkGeoIds(grid3d));
  assert(gridDoctor.checkJacobian(grid3d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Volume 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d = Teuchos::rcp(new connect3d<GEOSHAPE3D,ELMAP,NODEMAP>(commDev));
  connectGrid3d->setMesh3d(grid3d);
  connectGrid3d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Faces and edges cheking------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Check Faces" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> checkFaces(commDev);
  assert(checkFaces.check(grid3d->getFaces()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  if((pid == printPid) && verbose) {cout << "Check Edges" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(commDev);
  assert(checkEdges.check(grid3d->getEdges()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Distribute the elements2d----------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Distribute elements 2d" << " "; time(&start); cout << endl;}
    
  //Find local cadidates - the ones that has only one element connected
  GEOELEMENT2D  faceGeoItem(true);
  ELMAP         faceMapItem;
  
  GRAPH2D tempFaces3d = grid3d->getFaces();
  tempFaces3d.pushToGlobal();
  
  //Normal communication
  GEOELEMENT2D maxFace(true), minFace(true);
  GEOELEMENT2D maxEl2d(true), minEl2d(true);
  GEOELEMENT2D maxItem(true), minItem(true);
  
  pVectGlobalManip<GEOELEMENT2D,ELMAP> facesManipulator(commDev);
  facesManipulator.dataMinMax(tempFaces3d,minFace,maxFace);
  facesManipulator.dataMinMax(elements2d,minEl2d,maxEl2d);
  
  minItem = std::min(minFace,minEl2d);
  maxItem = std::max(maxFace,maxEl2d);  
  
  pVectComm<GEOELEMENT2D,ELMAP> facesCommunicator(commDev);
  facesCommunicator.vectorData(tempFaces3d,minItem,maxItem);
  facesCommunicator.vectorData(elements2d,minItem,maxItem);  
  
  //Ids matching -> elements2d creation
  typedef pair<GEOELEMENT2D,ELMAP>                    VALUETYPE;
  typedef typename map<GEOELEMENT2D,ELMAP>::iterator  ITERATOR2D;
  typedef pair<ITERATOR2D,bool>                       PAIR2D;
  
  map<GEOELEMENT2D,ELMAP>  elements2dSet;
  VALUETYPE                value2d;
  ITERATOR2D               iter;
  PAIR2D                   ret;
  
  for(UInt i=1; i <= elements2d.size(); ++i)
  {
    value2d.first  = elements2d.getItemL(i);
    value2d.second = elements2d.getRowMapL(i);
    
    ret = elements2dSet.insert(value2d);
    assert(ret.second == true); //The elements2d should be unique
  }
  
  elements2d.clear();
  UInt fLid = 1;
  
  for(UInt i=1; i <= tempFaces3d.size(); ++i)
  {
    faceGeoItem = tempFaces3d.getItemL(i);
    faceMapItem = tempFaces3d.getRowMapL(i);
    
    if(elements2dSet.count(faceGeoItem) == 1) //If not means that the face is not a boundary face
    {     
      iter = elements2dSet.find(faceGeoItem);
      
      faceGeoItem.setGeoId(iter->first.getGeoId());
      faceMapItem.setGid(iter->second.getGid());
      faceMapItem.setLid(fLid);
      
      fLid++;
      
      elements2d.push_back(faceGeoItem,faceMapItem);
    }
  }
  
  //Communicate back to processors
  facesCommunicator.vectorPid(elements2d);
  
  //Local indexing
  pGraphManip<GEOELEMENT2D,ELMAP,NODEMAP> local2dManip;
  local2dManip.setNormalIndexing(elements2d);
  
  //Push to local the elements2d
  elements2d.setColMap(grid3d->getNodes().getMapRef());
  elements2d.pushToLocal();
  
  
  //Create the grid2d
  grid2d = Teuchos::rcp(new mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>);
  grid2d->setNodes(grid3d->getNodes());
  grid2d->setElements(elements2d);
  grid2d->transferMap();
  grid2d->setMeshStandard(STDA);
  
  //Remove unused points
  meshDoctor2d<GEOSHAPE2D,ELMAP,NODEMAP> gridDoctor2d(commDev);
  gridDoctor2d.removeUnusedPoints(grid2d);
  
  //Grid2d checking
  mesh2dGlobalManip<GEOSHAPE2D,ELMAP,NODEMAP> gridManip2d(commDev);
  assert(gridManip2d.check(grid2d));
 
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Connecting2d-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting surface 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d = Teuchos::rcp(new connect2d<GEOSHAPE2D,ELMAP,NODEMAP>(commDev));
  connectGrid2d->setMesh2d(grid2d);
  connectGrid2d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}  
  
  //Boundary3d connecting--------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Boundary connecting 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d->setMesh2d(grid2d);
  connectGrid3d->buildBoundaryConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Flags------------------------------------------------------------------------------------------
  mesh3dCreated    = true;
  mesh2dCreated    = true;
  connect3dCreated = true;
  connect2dCreated = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
gmMesh_to_stdB(const string & meshfile, const bool & verbose)
{
  assert(staticAssert< GEOSHAPE::geoName == LINEARTETRA >::returnValue);
  assert(commDevLoaded);
  
  //Allocations------------------------------------------------------------------------------------
  UInt pid      = commDev->rank();
  UInt numPids  = commDev->size();
  UInt printPid = 0;
  time_t start, end;
  
  //Allocations
  NODESVECT nodes;
  GRAPH3D   elements3d;
  GRAPH2D   elements2d;
  
  
  //Mesh loading-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Start mesh loading" << " "; time(&start); cout << endl;}
  loadMesh<GEOSHAPE,ELMAP,NODEMAP> meshLoader;
  meshLoader.setParam(pid,printPid,numPids);
  meshLoader.gmMeshParallel(meshfile,commDev,nodes,elements3d,elements2d);
  
  elements2d.colIsLocal() = false;
  elements3d.colIsLocal() = false;
  
  commDev->barrier();
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //The grid3d-------------------------------------------------------------------------------------
  grid3d = Teuchos::rcp(new mesh3d<GEOSHAPE,ELMAP,NODEMAP>);
  grid3d->setNodes(nodes);
  grid3d->setElements(elements3d);
  grid3d->transferMap();
  grid3d->setMeshStandard(STDL);
  
  
  //Mesh partition---------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh partition" << " "; time(&start); cout << endl;}
  
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip3d(commDev);
  gridManip3d.meshPartition(grid3d);
  assert(gridManip3d.check(grid3d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh overlapping-------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh cheking" << " "; time(&start); cout << endl;}
  
  meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP> gridDoctor(commDev);
  gridDoctor.fixGeoIds(grid3d);
  assert(gridManip3d.check(grid3d));
  assert(gridDoctor.checkGeoIds(grid3d));
  assert(gridDoctor.checkJacobian(grid3d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Volume 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d = Teuchos::rcp(new connect3d<GEOSHAPE3D,ELMAP,NODEMAP>(commDev));
  connectGrid3d->setMesh3d(grid3d);
  connectGrid3d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Faces and edges cheking------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Check Faces" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> checkFaces(commDev);
  assert(checkFaces.check(grid3d->getFaces()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  if((pid == printPid) && verbose) {cout << "Check Edges" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(commDev);
  assert(checkEdges.check(grid3d->getEdges()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
    
  
  //Distribute the elements2d----------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Distribute elements 2d" << " "; time(&start); cout << endl;}
    
  //Find local cadidates - the ones that has only one element connected
  GEOELEMENT2D  faceGeoItem(true);
  ELMAP         faceMapItem;
  
  GRAPH2D tempFaces3d = grid3d->getFaces();
  tempFaces3d.pushToGlobal();
  
  //Normal communication
  GEOELEMENT2D maxFace(true), minFace(true);
  GEOELEMENT2D maxEl2d(true), minEl2d(true);
  GEOELEMENT2D maxItem(true), minItem(true);
  
  pVectGlobalManip<GEOELEMENT2D,ELMAP> facesManipulator(commDev);
  facesManipulator.dataMinMax(tempFaces3d,minFace,maxFace);
  facesManipulator.dataMinMax(elements2d,minEl2d,maxEl2d);
  
  minItem = std::min(minFace,minEl2d);
  maxItem = std::max(maxFace,maxEl2d);  
  
  pVectComm<GEOELEMENT2D,ELMAP> facesCommunicator(commDev);
  facesCommunicator.vectorData(tempFaces3d,minItem,maxItem);
  facesCommunicator.vectorData(elements2d,minItem,maxItem);  
  
  //Ids matching -> elements2d creation
  typedef pair<GEOELEMENT2D,ELMAP>                    VALUETYPE;
  typedef typename map<GEOELEMENT2D,ELMAP>::iterator  ITERATOR2D;
  typedef pair<ITERATOR2D,bool>                       PAIR2D;
  
  map<GEOELEMENT2D,ELMAP>  elements2dSet;
  VALUETYPE                value2d;
  ITERATOR2D               iter;
  PAIR2D                   ret;
  
  for(UInt i=1; i <= elements2d.size(); ++i)
  {
    value2d.first  = elements2d.getItemL(i);
    value2d.second = elements2d.getRowMapL(i);
    
    ret = elements2dSet.insert(value2d);
    assert(ret.second == true); //The elements2d should be unique
  }
  
  elements2d.clear();
  UInt fLid = 1;
  
  for(UInt i=1; i <= tempFaces3d.size(); ++i)
  {
    faceGeoItem = tempFaces3d.getItemL(i);
    faceMapItem = tempFaces3d.getRowMapL(i);
    
    if(elements2dSet.count(faceGeoItem) == 1)
    {     
      iter = elements2dSet.find(faceGeoItem);
      
      faceGeoItem.setGeoId(iter->first.getGeoId());
      faceMapItem.setGid(iter->second.getGid());
      faceMapItem.setLid(fLid);
      
      fLid++;
      
      elements2d.push_back(faceGeoItem,faceMapItem);
    }
  }
  
  //Communicate back to processors
  facesCommunicator.vectorPid(elements2d);
  
  //Local indexing
  pGraphManip<GEOELEMENT2D,ELMAP,NODEMAP> local2dManip;
  local2dManip.setNormalIndexing(elements2d);
  
  //Push to local the elements2d
  elements2d.setColMap(grid3d->getNodes().getMapRef());
  elements2d.pushToLocal();
  
  
  //Create the grid2d
  grid2d = Teuchos::rcp(new mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>);
  grid2d->setNodes(grid3d->getNodes());
  grid2d->setElements(elements2d);
  grid2d->transferMap();
  grid2d->setMeshStandard(STDB);
  
  //Remove unused points
  meshDoctor2d<GEOSHAPE2D,ELMAP,NODEMAP> gridDoctor2d(commDev);
  gridDoctor2d.removeUnusedPoints(grid2d);
  
  //Grid2d checking
  mesh2dGlobalManip<GEOSHAPE2D,ELMAP,NODEMAP> gridManip2d(commDev);
  assert(gridManip2d.check(grid2d));
 
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Connecting2d-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting surface 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d = Teuchos::rcp(new connect2d<GEOSHAPE2D,ELMAP,NODEMAP>(commDev));
  connectGrid2d->setMesh2d(grid2d);
  connectGrid2d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}  
  
  
  //Boundary3d connecting--------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Boundary connecting 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d->setMesh2d(grid2d);
  connectGrid3d->buildBoundaryConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Flags------------------------------------------------------------------------------------------
  mesh3dCreated    = true;
  mesh2dCreated    = true;
  connect3dCreated = true;
  connect2dCreated = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
femap_to_stdA(const string & meshfile, const string & colorfile, const bool & verbose)
{
  assert(staticAssert< GEOSHAPE::geoName == LINEARHEXA >::returnValue);
  assert(commDevLoaded);
  
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
  GRAPH3D   elements3d;
  COLORVECT nodesColor;
  
  loadMesh<GEOSHAPE,ELMAP,NODEMAP> meshLoader;
  meshLoader.setParam(pid,printPid,numPids);
  meshLoader.nastran(meshfile,nodes,elements3d,nodesColor);
  
  elements3d.colIsLocal() = false;
  
  commDev->barrier();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //The grid3d-------------------------------------------------------------------------------------
  grid3d = rcp(new mesh3d<GEOSHAPE,ELMAP,NODEMAP>);
  grid3d->setNodes(nodes);
  grid3d->setElements(elements3d);
  grid3d->transferMap();
  grid3d->setMeshStandard(STDL);
  
  
  //Mesh partition---------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh partition" << " "; time(&start); cout << endl;}
  
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip3d(commDev);
  gridManip3d.meshPartition(grid3d);
  assert(gridManip3d.check(grid3d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh overlapping-------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh overlap" << " "; time(&start); cout << endl;}
  
  meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP> gridDoctor(commDev);
  gridManip3d.meshOverlapping(grid3d);
  assert(gridManip3d.check(grid3d));
  assert(gridDoctor.checkGeoIds(grid3d));
  assert(gridDoctor.checkJacobian(grid3d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Volume 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d = rcp(new connect3d<GEOSHAPE3D,ELMAP,NODEMAP>(commDev));
  connectGrid3d->setMesh3d(grid3d);
  connectGrid3d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Faces and edges cheking------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Check Faces" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> checkFaces(commDev);
  assert(checkFaces.check(grid3d->getFaces()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  if((pid == printPid) && verbose) {cout << "Check Edges" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(commDev);
  assert(checkEdges.check(grid3d->getEdges()));
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
  
  
  //Generate 2d elements---------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Generate elements 2d" << " "; time(&start); cout << endl;}  
  
  //Remap of the nodesColor
  pVectGlobalManip<UInt,NODEMAP> nodesColorComm(commDev);
  
  pMap<NODEMAP> newColorMap = grid3d->getNodes().getMapRef();
  nodesColorComm.changeMap(nodesColor,newColorMap);
  
  //Global faces generation
  pGraph<pGraphItem,ELMAP,ELMAP> elementsToFaces = connectGrid3d->getElementToFace();
  pGraph<pGraphItem,ELMAP,ELMAP> globalFacesToElments;
  
  pGraphGlobalManip<pGraphItem,ELMAP,ELMAP>  tempGraphManip(commDev);
  tempGraphManip.inversion(elementsToFaces,globalFacesToElments);
  
  //Get one-connected faces
  typedef typename MESH3D::GRAPH2D GRAPH2D;
  GRAPH2D      elements2d;
  GEOELEMENT2D el2d;
  ELMAP        el2dMap;
  UInt minPower, nodeLid;
  UInt tot = 1;
  
  assert(grid3d->getFaces().size() == globalFacesToElments.size());
  
  for(UInt i=1; i <= grid3d->getFaces().size(); ++i)
  {
    assert((globalFacesToElments(i).size() == 1) || (globalFacesToElments(i).size() == 2));
    
    if(globalFacesToElments(i).size() == 1)
    {
      el2d    = grid3d->getFaces().getItemL(i);
      el2dMap = grid3d->getFaces().getRowMapL(i);
      
      assert(el2d.size() == 4);
      
      //GeoId determination
      nodeLid  = el2d.getCid(1);      
      color    = nodesColor(nodeLid);
      minColor = color;
      minPower = powerColor(color);
      
      for(UInt j=1; j <= el2d.size(); ++j)
      {
        nodeLid  = el2d.getCid(j);
        color    = nodesColor(nodeLid);
  
        if(powerColor(color) < minPower)
        {
          minColor = color;
          minPower = powerColor(color);
        }
      }
      
      el2d.setGeoId(minColor);
      
      //Set local Id
      el2dMap.setPid(pid);
      el2dMap.setLid(tot);
      tot++;
      
      //Add the 2d element
      elements2d.push_back(el2dMap,el2d);
    }
  }
  
  elements2d.setColMap(grid3d->getNodes().getMapRef());
  
  //Fixing the ownership and global nombering of the 2d elements
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> BelementsManip(commDev);
  BelementsManip.buildGlobalNumbering(elements2d);
  
  //Create the grid2d
  grid2d = rcp(new mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>);
  grid2d->setNodes(grid3d->getNodes());
  grid2d->setElements(elements2d);
  grid2d->transferMap();
  grid2d->setMeshStandard(STDA);
  
  //Remove unused points
  meshDoctor2d<GEOSHAPE2D,ELMAP,NODEMAP> gridDoctor2d(commDev);
  gridDoctor2d.removeUnusedPoints(grid2d);
 
  //Grid2d checking
  mesh2dGlobalManip<GEOSHAPE2D,ELMAP,NODEMAP> gridManip2d(commDev);
  assert(gridManip2d.check(grid2d));
 
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Connecting2d-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting surface 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d = rcp(new connect2d<GEOSHAPE2D,ELMAP,NODEMAP>(commDev));
  connectGrid2d->setMesh2d(grid2d);
  connectGrid2d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}  
  
  
  //Boundary3d connecting--------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Boundary connecting 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d->setMesh2d(grid2d);
  connectGrid3d->buildBoundaryConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Flags------------------------------------------------------------------------------------------
  mesh3dCreated    = true;
  mesh2dCreated    = true;
  connect3dCreated = true;
  connect2dCreated = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
femap_to_stdB(const string & meshfile, const string & colorfile, const bool & verbose)
{
  assert(staticAssert< GEOSHAPE::geoName == LINEARHEXA >::returnValue);
  assert(commDevLoaded);
  
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
  GRAPH3D   elements3d;
  COLORVECT nodesColor;
  
  loadMesh<GEOSHAPE,ELMAP,NODEMAP> meshLoader;
  meshLoader.setParam(pid,printPid,numPids);
  meshLoader.nastran(meshfile,nodes,elements3d,nodesColor);
  
  elements3d.colIsLocal() = false;
  
  commDev->barrier();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //The grid3d-------------------------------------------------------------------------------------
  grid3d = rcp(new mesh3d<GEOSHAPE,ELMAP,NODEMAP>);
  grid3d->setNodes(nodes);
  grid3d->setElements(elements3d);
  grid3d->transferMap();
  grid3d->setMeshStandard(STDL);
  
  
  //Mesh partition---------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh partition" << " "; time(&start); cout << endl;}
  
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip3d(commDev);
  gridManip3d.meshPartition(grid3d);
  assert(gridManip3d.check(grid3d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Volume 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d = rcp(new connect3d<GEOSHAPE3D,ELMAP,NODEMAP>(commDev));
  connectGrid3d->setMesh3d(grid3d);
  connectGrid3d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Faces and edges cheking------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Check Faces" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> checkFaces(commDev);
  assert(checkFaces.check(grid3d->getFaces()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  if((pid == printPid) && verbose) {cout << "Check Edges" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(commDev);
  assert(checkEdges.check(grid3d->getEdges()));
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
  
  
  //Generate 2d elements---------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Generate elements 2d" << " "; time(&start); cout << endl;}
  
  //Remap of the nodesColor
  pVectGlobalManip<UInt,NODEMAP> nodesColorComm(commDev);
  
  pMap<NODEMAP> newColorMap = grid3d->getNodes().getMapRef();
  nodesColorComm.changeMap(nodesColor,newColorMap);
  
  //Global faces generation
  pGraph<pGraphItem,ELMAP,ELMAP> elementsToFaces = connectGrid3d->getElementToFace();
  pGraph<pGraphItem,ELMAP,ELMAP> globalFacesToElments;
  
  pGraphGlobalManip<pGraphItem,ELMAP,ELMAP>  tempGraphManip(commDev);
  tempGraphManip.inversion(elementsToFaces,globalFacesToElments);
  
  //Get one-connected faces
  typedef typename MESH3D::GRAPH2D GRAPH2D;
  GRAPH2D      elements2d;
  GEOELEMENT2D el2d;
  ELMAP        el2dMap;
  UInt minPower, nodeLid;
  UInt tot = 1;
  
  assert(grid3d->getFaces().size() == globalFacesToElments.size());
  
  for(UInt i=1; i <= grid3d->getFaces().size(); ++i)
  {
    assert((globalFacesToElments(i).size() == 1) || (globalFacesToElments(i).size() == 2));
    
    if(globalFacesToElments(i).size() == 1)
    {
      el2d    = grid3d->getFaces().getItemL(i);
      el2dMap = grid3d->getFaces().getRowMapL(i);
      
      assert(el2d.size() == 4);
      
      //GeoId determination
      nodeLid  = el2d.getCid(1);
      color    = nodesColor(nodeLid);
      minColor = color;
      minPower = powerColor(color);
      
      for(UInt j=1; j <= el2d.size(); ++j)
      {
	nodeLid  = el2d.getCid(j);
        color    = nodesColor(nodeLid);
	
	if(powerColor(color) < minPower)
	{
	  minColor = color;
          minPower = powerColor(color);
	}
      }
      
      el2d.setGeoId(minColor);
      
      //Set local Id
      el2dMap.setPid(pid);
      el2dMap.setLid(tot);
      tot++;
      
      //Add the 2d element
      elements2d.push_back(el2dMap,el2d);
    }
  }
  
  elements2d.setColMap(grid3d->getNodes().getMapRef());
  
  //Fixing the ownership and global nombering of the 2d elements
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> BelementsManip(commDev);
  BelementsManip.buildGlobalNumbering(elements2d);
  
  //Create the grid2d
  grid2d = rcp(new mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>);
  grid2d->setNodes(grid3d->getNodes());
  grid2d->setElements(elements2d);
  grid2d->transferMap();
  grid2d->setMeshStandard(STDB);
  
  //Remove unused points
  meshDoctor2d<GEOSHAPE2D,ELMAP,NODEMAP> gridDoctor2d(commDev);
  gridDoctor2d.removeUnusedPoints(grid2d);
 
  //Grid2d checking
  mesh2dGlobalManip<GEOSHAPE2D,ELMAP,NODEMAP> gridManip2d(commDev);
  assert(gridManip2d.check(grid2d));
 
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Connecting2d-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting surface 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d = rcp(new connect2d<GEOSHAPE2D,ELMAP,NODEMAP>(commDev));
  connectGrid2d->setMesh2d(grid2d);
  connectGrid2d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}  
  
  
  //Boundary3d connecting--------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Boundary connecting 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d->setMesh2d(grid2d);
  connectGrid3d->buildBoundaryConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Flags------------------------------------------------------------------------------------------
  mesh3dCreated    = true;
  mesh2dCreated    = true;
  connect3dCreated = true;
  connect2dCreated = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
neutral_to_stdA(const string & meshfile, const string & colorfile, const bool & verbose)
{
  assert(staticAssert< GEOSHAPE::geoName == LINEARHEXA >::returnValue);
  assert(commDevLoaded);
  
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
  GRAPH3D   elements3d;
  COLORVECT nodesColor;
  
  loadMesh<GEOSHAPE,ELMAP,NODEMAP> meshLoader;
  meshLoader.setParam(pid,printPid,numPids);
  meshLoader.neutral(meshfile,nodes,elements3d,nodesColor);
  
  elements3d.colIsLocal() = false;
  
  commDev->barrier();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //The grid3d-------------------------------------------------------------------------------------
  grid3d = rcp(new mesh3d<GEOSHAPE,ELMAP,NODEMAP>);
  grid3d->setNodes(nodes);
  grid3d->setElements(elements3d);
  grid3d->transferMap();
  grid3d->setMeshStandard(STDL);
  
  
  //Mesh partition---------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh partition" << " "; time(&start); cout << endl;}
  
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip3d(commDev);
  gridManip3d.meshPartition(grid3d);
  assert(gridManip3d.check(grid3d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh overlapping-------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh overlap" << " "; time(&start); cout << endl;}
  
  meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP> gridDoctor(commDev);
  gridManip3d.meshOverlapping(grid3d);
  assert(gridManip3d.check(grid3d));
  assert(gridDoctor.checkGeoIds(grid3d));
  assert(gridDoctor.checkJacobian(grid3d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Volume 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d = rcp(new connect3d<GEOSHAPE3D,ELMAP,NODEMAP>(commDev));
  connectGrid3d->setMesh3d(grid3d);
  connectGrid3d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Faces and edges cheking------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Check Faces" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> checkFaces(commDev);
  assert(checkFaces.check(grid3d->getFaces()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  if((pid == printPid) && verbose) {cout << "Check Edges" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(commDev);
  assert(checkEdges.check(grid3d->getEdges()));
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
  
  
  //Generate 2d elements---------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Generate elements 2d" << " "; time(&start); cout << endl;}  
  
  //Remap of the nodesColor
  pVectGlobalManip<UInt,NODEMAP> nodesColorComm(commDev);
  
  pMap<NODEMAP> newColorMap = grid3d->getNodes().getMapRef();
  nodesColorComm.changeMap(nodesColor,newColorMap);
  
  //Global faces generation
  pGraph<pGraphItem,ELMAP,ELMAP> elementsToFaces = connectGrid3d->getElementToFace();
  pGraph<pGraphItem,ELMAP,ELMAP> globalFacesToElments;
  
  pGraphGlobalManip<pGraphItem,ELMAP,ELMAP>  tempGraphManip(commDev);
  tempGraphManip.inversion(elementsToFaces,globalFacesToElments);
  
  //Get one-connected faces
  typedef typename MESH3D::GRAPH2D GRAPH2D;
  GRAPH2D      elements2d;
  GEOELEMENT2D el2d;
  ELMAP        el2dMap;
  UInt minPower, nodeLid;
  UInt tot = 1;
  
  assert(grid3d->getFaces().size() == globalFacesToElments.size());
  
  for(UInt i=1; i <= grid3d->getFaces().size(); ++i)
  {
    assert((globalFacesToElments(i).size() == 1) || (globalFacesToElments(i).size() == 2));
    
    if(globalFacesToElments(i).size() == 1)
    {
      el2d    = grid3d->getFaces().getItemL(i);
      el2dMap = grid3d->getFaces().getRowMapL(i);
      
      assert(el2d.size() == 4);
      
      //GeoId determination
      nodeLid  = el2d.getCid(1);      
      color    = nodesColor(nodeLid);
      minColor = color;
      minPower = powerColor(color);
      
      for(UInt j=1; j <= el2d.size(); ++j)
      {
	nodeLid  = el2d.getCid(j);
        color    = nodesColor(nodeLid);
	
	if(powerColor(color) < minPower)
	{
	  minColor = color;
          minPower = powerColor(color);
	}
      }
      
      el2d.setGeoId(minColor);
      
      //Set local Id
      el2dMap.setPid(pid);
      el2dMap.setLid(tot);
      tot++;
      
      //Add the 2d element
      elements2d.push_back(el2dMap,el2d);
    }
  }
  
  elements2d.setColMap(grid3d->getNodes().getMapRef());
  
  //Fixing the ownership and global nombering of the 2d elements
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> BelementsManip(commDev);
  BelementsManip.buildGlobalNumbering(elements2d);
  
  //Create the grid2d
  grid2d = rcp(new mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>);
  grid2d->setNodes(grid3d->getNodes());
  grid2d->setElements(elements2d);
  grid2d->transferMap();
  grid2d->setMeshStandard(STDA);
  
  //Remove unused points
  meshDoctor2d<GEOSHAPE2D,ELMAP,NODEMAP> gridDoctor2d(commDev);
  gridDoctor2d.removeUnusedPoints(grid2d);
 
  //Grid2d checking
  mesh2dGlobalManip<GEOSHAPE2D,ELMAP,NODEMAP> gridManip2d(commDev);
  assert(gridManip2d.check(grid2d));
 
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Connecting2d-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting surface 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d = rcp(new connect2d<GEOSHAPE2D,ELMAP,NODEMAP>(commDev));
  connectGrid2d->setMesh2d(grid2d);
  connectGrid2d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}  
  
  
  //Boundary3d connecting--------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Boundary connecting 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d->setMesh2d(grid2d);
  connectGrid3d->buildBoundaryConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Flags------------------------------------------------------------------------------------------
  mesh3dCreated    = true;
  mesh2dCreated    = true;
  connect3dCreated = true;
  connect2dCreated = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
neutral_to_stdB(const string & meshfile, const string & colorfile, const bool & verbose)
{
  assert(staticAssert< GEOSHAPE::geoName == LINEARHEXA >::returnValue);
  assert(commDevLoaded);
  
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
  GRAPH3D   elements3d;
  COLORVECT nodesColor;
  
  loadMesh<GEOSHAPE,ELMAP,NODEMAP> meshLoader;
  meshLoader.setParam(pid,printPid,numPids);
  meshLoader.neutral(meshfile,nodes,elements3d,nodesColor);
  
  elements3d.colIsLocal() = false;
  
  commDev->barrier();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //The grid3d-------------------------------------------------------------------------------------
  grid3d = rcp(new mesh3d<GEOSHAPE,ELMAP,NODEMAP>);
  grid3d->setNodes(nodes);
  grid3d->setElements(elements3d);
  grid3d->transferMap();
  grid3d->setMeshStandard(STDL);
  
  
  //Mesh partition---------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh partition" << " "; time(&start); cout << endl;}
  
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip3d(commDev);
  gridManip3d.meshPartition(grid3d);
  assert(gridManip3d.check(grid3d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Volume 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d = rcp(new connect3d<GEOSHAPE3D,ELMAP,NODEMAP>(commDev));
  connectGrid3d->setMesh3d(grid3d);
  connectGrid3d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Faces and edges cheking------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Check Faces" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> checkFaces(commDev);
  assert(checkFaces.check(grid3d->getFaces()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  if((pid == printPid) && verbose) {cout << "Check Edges" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(commDev);
  assert(checkEdges.check(grid3d->getEdges()));
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
  
  
  //Generate 2d elements---------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Generate elements 2d" << " "; time(&start); cout << endl;}
  
  //Remap of the nodesColor
  pVectGlobalManip<UInt,NODEMAP> nodesColorComm(commDev);
  
  pMap<NODEMAP> newColorMap = grid3d->getNodes().getMapRef();
  nodesColorComm.changeMap(nodesColor,newColorMap);
  
  //Global faces generation
  pGraph<pGraphItem,ELMAP,ELMAP> elementsToFaces = connectGrid3d->getElementToFace();
  pGraph<pGraphItem,ELMAP,ELMAP> globalFacesToElments;
  
  pGraphGlobalManip<pGraphItem,ELMAP,ELMAP>  tempGraphManip(commDev);
  tempGraphManip.inversion(elementsToFaces,globalFacesToElments);
  
  //Get one-connected faces
  typedef typename MESH3D::GRAPH2D GRAPH2D;
  GRAPH2D      elements2d;
  GEOELEMENT2D el2d;
  ELMAP        el2dMap;
  UInt minPower, nodeLid;
  UInt tot = 1;
  
  assert(grid3d->getFaces().size() == globalFacesToElments.size());
  
  for(UInt i=1; i <= grid3d->getFaces().size(); ++i)
  {
    assert((globalFacesToElments(i).size() == 1) || (globalFacesToElments(i).size() == 2));
    
    if(globalFacesToElments(i).size() == 1)
    {
      el2d    = grid3d->getFaces().getItemL(i);
      el2dMap = grid3d->getFaces().getRowMapL(i);
      
      assert(el2d.size() == 4);
      
      //GeoId determination
      nodeLid  = el2d.getCid(1);
      color    = nodesColor(nodeLid);
      minColor = color;
      minPower = powerColor(color);
      
      for(UInt j=1; j <= el2d.size(); ++j)
      {
	nodeLid  = el2d.getCid(j);
        color    = nodesColor(nodeLid);
	
	if(powerColor(color) < minPower)
	{
	  minColor = color;
          minPower = powerColor(color);
	}
      }
      
      el2d.setGeoId(minColor);
      
      //Set local Id
      el2dMap.setPid(pid);
      el2dMap.setLid(tot);
      tot++;
      
      //Add the 2d element
      elements2d.push_back(el2dMap,el2d);
    }
  }
  
  elements2d.setColMap(grid3d->getNodes().getMapRef());
  
  //Fixing the ownership and global nombering of the 2d elements
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> BelementsManip(commDev);
  BelementsManip.buildGlobalNumbering(elements2d);
  
  //Create the grid2d
  grid2d = rcp(new mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>);
  grid2d->setNodes(grid3d->getNodes());
  grid2d->setElements(elements2d);
  grid2d->transferMap();
  grid2d->setMeshStandard(STDB);
  
  //Remove unused points
  meshDoctor2d<GEOSHAPE2D,ELMAP,NODEMAP> gridDoctor2d(commDev);
  gridDoctor2d.removeUnusedPoints(grid2d);
 
  //Grid2d checking
  mesh2dGlobalManip<GEOSHAPE2D,ELMAP,NODEMAP> gridManip2d(commDev);
  assert(gridManip2d.check(grid2d));
 
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Connecting2d-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting surface 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d = rcp(new connect2d<GEOSHAPE2D,ELMAP,NODEMAP>(commDev));
  connectGrid2d->setMesh2d(grid2d);
  connectGrid2d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}  
  
  
  //Boundary3d connecting--------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Boundary connecting 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d->setMesh2d(grid2d);
  connectGrid3d->buildBoundaryConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Flags------------------------------------------------------------------------------------------
  mesh3dCreated    = true;
  mesh2dCreated    = true;
  connect3dCreated = true;
  connect2dCreated = true;
}


//_______________________________________________________________________________________________________
// MESH RE-INIT
//-------------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
reinit_stdA(const Teuchos::RCP<MESH2D> & inGrid2d,
            const Teuchos::RCP<MESH3D> & inGrid3d,
            const set<UInt>            & activeGeoIds,
            const bool                 & verbose)
{
  //Assert-----------------------------------------------------------------------------------------
  assert(inGrid3d->getMeshStandard() == STDA);
  
  //Allocations------------------------------------------------------------------------------------
  UInt pid      = commDev->rank();
  UInt numPids  = commDev->size();
  UInt printPid = 0;
  time_t start, end;
  
  //Create elements2d list-------------------------------------------------------------------------
  typedef std::pair<point3d,UInt>          PAIR;
  typedef std::map<point3d,UInt>::iterator ITER;
  
  GRAPH2D elements2d = inGrid2d->getElements();
  std::map<point3d,UInt> nodes3dList;
  
  PAIR pair;
  ITER iter2d;
  point3d P;
  UInt cid;
  
  for(UInt i=1; i <= inGrid3d->getNumNodes(); ++i)
  { 
    pair.first  = inGrid3d->getNodeL(i);
    pair.second = inGrid3d->getNodes().getMapL(i).getGid();    
    nodes3dList.insert(pair);
  }
  
  for(UInt i=1; i <= elements2d.size(); ++i)
  {
    for(UInt j=1; j <= elements2d.getItemL(i).size(); ++j)
    {
      cid    = elements2d.getItemL(i).getCid(j);
      P      = inGrid2d->getNodeL(cid);
      iter2d = nodes3dList.find(P);
      
      elements2d.getItemL(i).setCid(j, iter2d->second);
    }
  }
  
  elements2d.colIsLocal() = false;
  
  
  //Destroy overlap--------------------------------------------------------------------------------
  *grid3d = *inGrid3d;
  
  if((pid == printPid) && verbose) {cout << "Destroy overlap" << " "; time(&start); cout << endl;}
  
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip3d(commDev);
  gridManip3d.destroyOverlap(grid3d);
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Compute weights--------------------------------------------------------------------------------
  sVect<UInt> elWeights(grid3d->getNumElements());
  UInt geoId, count;
  
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    geoId = grid3d->getElementL(i).getGeoId();
    count = activeGeoIds.count(geoId);
    
    elWeights(i) = count * 999 + 1;
  }
  
  //Mesh re-partition------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh balancing" << " "; time(&start); cout << endl;}
  
  gridManip3d.meshBalancing(grid3d,elWeights);
  assert(gridManip3d.check(grid3d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh overlapping-------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh overlap" << " "; time(&start); cout << endl;}
  
  gridManip3d.meshOverlapping(grid3d);
  
  meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP> gridDoctor(commDev);
  gridDoctor.fixGeoIds(grid3d);
  assert(gridManip3d.check(grid3d));
  assert(gridDoctor.checkGeoIds(grid3d));
  assert(gridDoctor.checkJacobian(grid3d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Volume 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d = Teuchos::rcp(new connect3d<GEOSHAPE3D,ELMAP,NODEMAP>(commDev));
  connectGrid3d->setMesh3d(grid3d);
  connectGrid3d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Faces and edges cheking------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Check Faces" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> checkFaces(commDev);
  assert(checkFaces.check(grid3d->getFaces()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  if((pid == printPid) && verbose) {cout << "Check Edges" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(commDev);
  assert(checkEdges.check(grid3d->getEdges()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Distribute the elements2d----------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Distribute elements 2d" << " "; time(&start); cout << endl;}
    
  //Find local cadidates - the ones that has only one element connected
  GEOELEMENT2D  faceGeoItem(true);
  ELMAP         faceMapItem;
  
  GRAPH2D tempFaces3d = grid3d->getFaces();
  tempFaces3d.pushToGlobal();
  
  //Normal communication
  GEOELEMENT2D maxFace(true), minFace(true);
  
  pVectGlobalManip<GEOELEMENT2D,ELMAP> facesManipulator(commDev);
  facesManipulator.dataMinMax(tempFaces3d,minFace,maxFace); 
  
  pVectComm<GEOELEMENT2D,ELMAP> facesCommunicator(commDev);
  facesCommunicator.vectorData(tempFaces3d,minFace,maxFace);
  facesCommunicator.vectorData(elements2d,minFace,maxFace);  
  
  //Ids matching -> elements2d creation
  typedef std::pair<GEOELEMENT2D,ELMAP>               VALUETYPE;
  typedef typename map<GEOELEMENT2D,ELMAP>::iterator ITERATOR2D;
  typedef std::pair<ITERATOR2D,bool>                     PAIR2D;
  
  map<GEOELEMENT2D,ELMAP>  elements2dSet;
  VALUETYPE                value2d;
  ITERATOR2D               iter;
  PAIR2D                   ret;
  
  for(UInt i=1; i <= elements2d.size(); ++i)
  {
    value2d.first  = elements2d.getItemL(i);
    value2d.second = elements2d.getRowMapL(i);
    
    ret = elements2dSet.insert(value2d);
  }
  
  elements2d.clear();
  UInt fLid = 1;
  
  for(UInt i=1; i <= tempFaces3d.size(); ++i)
  {
    faceGeoItem = tempFaces3d.getItemL(i);
    faceMapItem = tempFaces3d.getRowMapL(i);
    
    if(elements2dSet.count(faceGeoItem) == 1) //If not means that the face is not a boundary face
    {     
      iter = elements2dSet.find(faceGeoItem);
      
      faceGeoItem.setGeoId(iter->first.getGeoId());
      faceMapItem.setGid(iter->second.getGid());
      faceMapItem.setLid(fLid);
      
      fLid++;
      
      elements2d.push_back(faceGeoItem,faceMapItem);
    }
  }
  
  //Communicate back to processors
  facesCommunicator.vectorPid(elements2d);
  
  //Local indexing
  pGraphManip<GEOELEMENT2D,ELMAP,NODEMAP> local2dManip;
  local2dManip.setNormalIndexing(elements2d);
  
  //Push to local the elements2d
  elements2d.setColMap(grid3d->getNodes().getMapRef());
  elements2d.pushToLocal();
  
  
  //Create the grid2d
  grid2d = Teuchos::rcp(new mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>);
  grid2d->setNodes(grid3d->getNodes());
  grid2d->setElements(elements2d);
  grid2d->transferMap();
  grid2d->setMeshStandard(STDA);
  
  //Remove unused points
  meshDoctor2d<GEOSHAPE2D,ELMAP,NODEMAP> gridDoctor2d(commDev);
  gridDoctor2d.removeUnusedPoints(grid2d);
  
  //Grid2d checking
  mesh2dGlobalManip<GEOSHAPE2D,ELMAP,NODEMAP> gridManip2d(commDev);
  assert(gridManip2d.check(grid2d));
 
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Connecting2d-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting surface 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d = Teuchos::rcp(new connect2d<GEOSHAPE2D,ELMAP,NODEMAP>(commDev));
  connectGrid2d->setMesh2d(grid2d);
  connectGrid2d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}  
  
  
  //Boundary3d connecting--------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Boundary connecting 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d->setMesh2d(grid2d);
  connectGrid3d->buildBoundaryConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Flags------------------------------------------------------------------------------------------
  mesh3dCreated    = true;
  mesh2dCreated    = true;
  connect3dCreated = true;
  connect2dCreated = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshInit3d<GEOSHAPE,ELMAP,NODEMAP>::
reinit_stdB(const Teuchos::RCP<MESH2D> & inGrid2d,
            const Teuchos::RCP<MESH3D> & inGrid3d,
            const set<UInt>            & activeGeoIds,
            const bool                 & verbose)
{
  //Assert-----------------------------------------------------------------------------------------
  assert(inGrid3d->getMeshStandard() == STDB);
  
  //Allocations------------------------------------------------------------------------------------
  UInt pid      = commDev->rank();
  UInt numPids  = commDev->size();
  UInt printPid = 0;
  time_t start, end;
  
  //Create elements2d list-------------------------------------------------------------------------
  typedef std::pair<point3d,UInt>          PAIR;
  typedef std::map<point3d,UInt>::iterator ITER;
  
  GRAPH2D elements2d = inGrid2d->getElements();
  std::map<point3d,UInt> nodes3dList;
  
  PAIR pair;
  ITER iter2d;
  point3d P;
  UInt cid;
  
  for(UInt i=1; i <= inGrid3d->getNumNodes(); ++i)
  { 
    pair.first  = inGrid3d->getNodeL(i);
    pair.second = inGrid3d->getNodes().getMapL(i).getGid();    
    nodes3dList.insert(pair);
  }
  
  for(UInt i=1; i <= elements2d.size(); ++i)
  {
    for(UInt j=1; j <= elements2d.getItemL(i).size(); ++j)
    {
      cid    = elements2d.getItemL(i).getCid(j);
      P      = inGrid2d->getNodeL(cid);
      iter2d = nodes3dList.find(P);
      
      elements2d.getItemL(i).setCid(j, iter2d->second);
    }
  }
  
  elements2d.colIsLocal() = false;
  
  //Compute weights--------------------------------------------------------------------------------
  sVect<UInt> elWeights(grid3d->getNumElements());
  UInt geoId, count;
  
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    geoId = grid3d->getElementL(i).getGeoId();
    count = activeGeoIds.count(geoId);
    
    elWeights(i) = count * 999 + 1;
  }  
  
  //Mesh re-partition------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Mesh balancing" << " "; time(&start); cout << endl;}
  
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> gridManip3d(commDev);
  gridManip3d.meshBalancing(grid3d,elWeights);
  assert(gridManip3d.check(grid3d));
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Mesh connecting--------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting Volume 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d = Teuchos::rcp(new connect3d<GEOSHAPE3D,ELMAP,NODEMAP>(commDev));
  connectGrid3d->setMesh3d(grid3d);
  connectGrid3d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Faces and edges cheking------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Check Faces" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> checkFaces(commDev);
  assert(checkFaces.check(grid3d->getFaces()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  if((pid == printPid) && verbose) {cout << "Check Edges" << " "; time(&start); cout << endl;}
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(commDev);
  assert(checkEdges.check(grid3d->getEdges()));
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Distribute the elements2d----------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Distribute elements 2d" << " "; time(&start); cout << endl;}
    
  //Find local cadidates - the ones that has only one element connected
  GEOELEMENT2D  faceGeoItem(true);
  ELMAP         faceMapItem;
  
  GRAPH2D tempFaces3d = grid3d->getFaces();
  tempFaces3d.pushToGlobal();
  
  //Normal communication
  GEOELEMENT2D maxFace(true), minFace(true);
  
  pVectGlobalManip<GEOELEMENT2D,ELMAP> facesManipulator(commDev);
  facesManipulator.dataMinMax(tempFaces3d,minFace,maxFace); 
  
  pVectComm<GEOELEMENT2D,ELMAP> facesCommunicator(commDev);
  facesCommunicator.vectorData(tempFaces3d,minFace,maxFace);
  facesCommunicator.vectorData(elements2d,minFace,maxFace);  
  
  //Ids matching -> elements2d creation
  typedef std::pair<GEOELEMENT2D,ELMAP>               VALUETYPE;
  typedef typename map<GEOELEMENT2D,ELMAP>::iterator ITERATOR2D;
  typedef std::pair<ITERATOR2D,bool>                     PAIR2D;
  
  map<GEOELEMENT2D,ELMAP>  elements2dSet;
  VALUETYPE                value2d;
  ITERATOR2D               iter;
  PAIR2D                   ret;
  
  for(UInt i=1; i <= elements2d.size(); ++i)
  {
    value2d.first  = elements2d.getItemL(i);
    value2d.second = elements2d.getRowMapL(i);
    
    ret = elements2dSet.insert(value2d);
  }
  
  elements2d.clear();
  UInt fLid = 1;
  
  for(UInt i=1; i <= tempFaces3d.size(); ++i)
  {
    faceGeoItem = tempFaces3d.getItemL(i);
    faceMapItem = tempFaces3d.getRowMapL(i);
    
    if(elements2dSet.count(faceGeoItem) == 1) //If not means that the face is not a boundary face
    {     
      iter = elements2dSet.find(faceGeoItem);
      
      faceGeoItem.setGeoId(iter->first.getGeoId());
      faceMapItem.setGid(iter->second.getGid());
      faceMapItem.setLid(fLid);
      
      fLid++;
      
      elements2d.push_back(faceGeoItem,faceMapItem);
    }
  }
  
  //Communicate back to processors
  facesCommunicator.vectorPid(elements2d);
  
  //Local indexing
  pGraphManip<GEOELEMENT2D,ELMAP,NODEMAP> local2dManip;
  local2dManip.setNormalIndexing(elements2d);
  
  //Push to local the elements2d
  elements2d.setColMap(grid3d->getNodes().getMapRef());
  elements2d.pushToLocal();
  
  
  //Create the grid2d
  grid2d = Teuchos::rcp(new mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>);
  grid2d->setNodes(grid3d->getNodes());
  grid2d->setElements(elements2d);
  grid2d->transferMap();
  grid2d->setMeshStandard(STDB);
  
  //Remove unused points
  meshDoctor2d<GEOSHAPE2D,ELMAP,NODEMAP> gridDoctor2d(commDev);
  gridDoctor2d.removeUnusedPoints(grid2d);
  
  //Grid2d checking
  mesh2dGlobalManip<GEOSHAPE2D,ELMAP,NODEMAP> gridManip2d(commDev);
  assert(gridManip2d.check(grid2d));
 
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  //Connecting2d-----------------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Connecting surface 2d" << " "; time(&start); cout << endl;}
  
  connectGrid2d = Teuchos::rcp(new connect2d<GEOSHAPE2D,ELMAP,NODEMAP>(commDev));
  connectGrid2d->setMesh2d(grid2d);
  connectGrid2d->buildConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}  
  
  
  //Boundary3d connecting--------------------------------------------------------------------------
  if((pid == printPid) && verbose) {cout << "Boundary connecting 3d" << " "; time(&start); cout << endl;}
  
  connectGrid3d->setMesh2d(grid2d);
  connectGrid3d->buildBoundaryConnectivity();
  
  if((pid == printPid) && verbose) {time(&end); cout << "done (" << difftime(end, start) << " s)" << endl << endl;}
  
  
  //Flags------------------------------------------------------------------------------------------
  mesh3dCreated    = true;
  mesh2dCreated    = true;
  connect3dCreated = true;
  connect2dCreated = true;
}

#endif
