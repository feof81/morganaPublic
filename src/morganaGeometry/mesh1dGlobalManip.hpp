/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef MESH1DGLOBALMANIP_HPP
#define MESH1DGLOBALMANIP_HPP

#include <map>
#include "parmetis.h"

#include "pVectGlobalManip.hpp"
#include "pGraphGlobalManip.hpp"
#include "mesh1d.hpp"
#include "traitsGeometry.hpp"


//_________________________________________________________________________________________________
// UNSPECIALIZED
//-------------------------------------------------------------------------------------------------

/*! Global mesh manipulator 1d - unspecialized, empty */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP> class mesh1dGlobalManip
{ };


//_________________________________________________________________________________________________
// PMAPITEM
//-------------------------------------------------------------------------------------------------

/*! Global mesh manipulator 1d - \c pMapItem */
template<typename GEOSHAPE, typename NODEMAP>
class mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItem ELMAP;
    typedef mesh1d<GEOSHAPE,ELMAP,NODEMAP> MESH1D;
    
    typedef typename GEOSHAPE::GEOBSHAPE   GEOSHAPE0D;
    
    typedef geoElement<GEOSHAPE>   GEOELEMENT;
    typedef geoElement<GEOSHAPE0D> GEOELEMENT0D;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
  
    /*! @name Constructors */ //@{
  public:
    mesh1dGlobalManip();
    mesh1dGlobalManip(const Teuchos::RCP<const communicator> & CommDev);
    mesh1dGlobalManip(const communicator & CommDev);
    //@}
    
    /*! @name Information functions */ //@{
  public:
    /*! Check the graphs, the equality of the col maps of the elements, faces and edges with the map of the points */
    bool check(const Teuchos::RCP<const MESH1D> & Grid1d) const;
    bool check(const MESH1D & Grid1d) const;
    
    /*! Return the total number of vertices on the processors */
    UInt getNumGlobalVertices(const Teuchos::RCP<const MESH1D> & Grid1d) const;
    UInt getNumGlobalVertices(const MESH1D & Grid1d) const;
    
    /*! Return the total numer of elements on the processors */
    UInt getNumGlobalElements(const Teuchos::RCP<const MESH1D> & Grid1d) const;
    UInt getNumGlobalElements(const MESH1D & Grid1d) const;
    //@}
    
    /*! @name Operation functions */ //@{
  public:
    /*! The mesh should be in global numbering and a new mesh in local numbering is built,
    no nodes or elements should be repeated across the processors */
    void meshPartition(const Teuchos::RCP<MESH1D> & Grid1d) const;
    void meshPartition(MESH1D & Grid1d) const;
    
    /*! The mesh should be in local numbering, no element repetition is allowed */
    void meshBalancing(const Teuchos::RCP<MESH1D> & Grid1d) const;
    void meshBalancing(MESH1D & grid1d) const;
    
    /*! The mesh should be in local numbering. The appendended elements and nodes are listed below the original ones */
    void meshOverlapping(const Teuchos::RCP<MESH1D> & Grid1d) const;
    void meshOverlapping(MESH1D & grid1d) const;
    
    /*! Reduces the mesh to a single process and pushes the \c cids to the global numbering.
    Duplicated elements or points are eliminated. */
    void gather(const UInt & gatherGid, const Teuchos::RCP<MESH1D> & targetGrid1d, const Teuchos::RCP<const MESH1D> & sourceGrid1d);
    void gather(const UInt & gatherGid, MESH1D & targetGrid1d, const MESH1D & sourceGrid1d);
    
    /*! Eliminates the overlapped regions. ATTENTION: since there is no way history of the way
    the mesh was overlapped the shared elements are put in a particular domain in an
    un-predictable manner */
    void destroyOverlap(const Teuchos::RCP<MESH1D> & grid1d) const;
    void destroyOverlap(MESH1D & grid1d) const;
    //@}
    
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MESH1D>       & OldMesh,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MESH1D>             & NewMesh);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MESH1D       & OldMesh,
                            const communicator & NewCommDev,
                                  MESH1D       & NewMesh);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MESH1D>       & OldMesh,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MESH1D>             & NewMesh);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MESH1D       & OldMesh,
                            const communicator & NewCommDev,
                                  MESH1D       & NewMesh);
    //@}
};


template<typename GEOSHAPE, typename NODEMAP>
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
mesh1dGlobalManip()
{
  assert(GEOSHAPE::nDim == 1);
  commDevLoaded = false;
}

template<typename GEOSHAPE, typename NODEMAP>
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
mesh1dGlobalManip(const Teuchos::RCP<const communicator> & CommDev)
{
  assert(GEOSHAPE::nDim == 1);
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename NODEMAP>
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
mesh1dGlobalManip(const communicator & CommDev)
{
  assert(GEOSHAPE::nDim == 1);
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename NODEMAP>
bool
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
check(const Teuchos::RCP<const MESH1D> & Grid1d) const
{
  return(check(*Grid1d));
}

template<typename GEOSHAPE, typename NODEMAP>
bool
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
check(const MESH1D & Grid1d) const
{
  assert(Grid1d.getMapTransferred());
  
  //Allocations
  bool flag;
  UInt pid  = commDev->rank();
  UInt size = commDev->size();
  
  //Nodes vect check
  pVectGlobalManip<point3d,NODEMAP> nodesChecking(commDev);
  bool flagN = nodesChecking.check(Grid1d.getNodes());
  
  if(!flagN)
  { cout << "ERROR: Nodes check failed. Pid: " << pid << endl; }
  
  //Elements check
  pGraphGlobalManip<GEOELEMENT,pMapItem,NODEMAP> elementsChecking(commDev);
  bool flagE = elementsChecking.check(Grid1d.getElements());
  
  if(!flagE)
  { cout << "ERROR: Elements check failed. Pid: " << pid << endl; }
  
  flag = flagN & flagE;
  
  //Cross map checking
  pMapManip<NODEMAP> crossMapChecker;
  
  //Checking against elements col map
  sVect<UInt> allFlags(size);
  UInt myFlag = crossMapChecker.isEqual(Grid1d.getNodes().getMapRef(),Grid1d.getElements().getColMap());
  
  if(myFlag != 1)
  { cout << "ERROR: Elements col map differs from Nodes col map. Pid: " << pid << endl; }
  
  all_gather(*commDev, myFlag, allFlags);
  
  for(UInt i=1; i <= size; ++i)
  { myFlag = myFlag & allFlags(i); }
  
  flag = flag & myFlag;
  
  //Cross check among processors
  sVect<UInt> flags;
  all_gather(*commDev,UInt(flag),flags);
  
  for(UInt i=1; i <= size; ++i)
  { flag = flag & bool(flags(i));}
  
  return(flag);
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalVertices(const Teuchos::RCP<const MESH1D> & Grid1d) const
{
  pVectGlobalManip<point3d,NODEMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid1d->getNodes()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalVertices(const MESH1D & Grid1d) const
{
  pVectGlobalManip<point3d,NODEMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid1d.getNodes()));
} 

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalElements(const Teuchos::RCP<const MESH1D> & Grid1d) const
{
  pVectGlobalManip<GEOELEMENT,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid1d->getElements()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalElements(const MESH1D & Grid1d) const
{
  pVectGlobalManip<GEOELEMENT,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid1d.getElements()));
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshPartition(const Teuchos::RCP<MESH1D> & Grid1d) const
{
  assert(commDevLoaded);
  meshPartition(*Grid1d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshPartition(mesh1d<GEOSHAPE,ELMAP,NODEMAP> & Grid1d) const
{
  //Asserts
  assert(commDevLoaded);
  assert(Grid1d.getElements().colIsLocal() == false);
  assert(Grid1d.getMeshStandard() == STDL);
  
  //Clearing
  Grid1d.transferMap();
  
  //PARMETIS_______________________________________________________________________________________
  UInt commSize = commDev->size();
  
  //Sizes communication
  UInt mySize = Grid1d.getNumElements();
  sVect<UInt> allSizes(commSize);
  
  all_gather(*commDev, mySize, allSizes);
  
  //Elements distribution
  idx_t * elmdist = new idx_t[commSize + 1];
  
  elmdist[0] = 0;
  
  for(UInt i=1 ; i <= commSize; ++i)
  {
    elmdist[i] = elmdist[i-1] + allSizes(i);
  }  
  
  //Element pointer
  idx_t * eptr = new idx_t[mySize + 1];
  
  eptr[0] = 0;
  
  for(UInt i=1 ; i <= mySize; ++i)
  {
    eptr[i] = eptr[i-1] + GEOSHAPE::numPoints;
  }
   
  //Element indices
  UInt  node, k = 0;
  idx_t * eind = new idx_t[mySize * GEOSHAPE::numPoints];
  
  for(UInt i=1; i <= mySize; ++i)
  {
    for(UInt j=1; j <= GEOSHAPE::numPoints; ++j)
    {
      node    = Grid1d.getElementL(i).getCid(j) - 1; //Change numbering      
      eind[k] = node;
      k++;
    }
    
    assert(k == UInt(eptr[i]));
  }
  
  //Some tandard flags
  idx_t * elmwgt = NULL;
  idx_t wgtflag = 0;
  idx_t numflag = 0;
  idx_t ncon    = 1;
  
  //Other flags
  idx_t ncommonnodes = 1;
  idx_t nparts       = commSize;
  
  //Weights
  real_t * tpwgts = new real_t[ncon * nparts];
  
  for(int i=0; i < nparts; ++i)
  {
    tpwgts[i] = Real(1.0) / Real(nparts);
  }
  
  //Balancing index 
  real_t ubvec = 1.05;
  
  //Options Vector
  idx_t * options = new idx_t[3];
  options[0] = 0;
  options[1] = 0;
  options[2] = 0;
  
  //Mpi comm 
  MPI_Comm mpiComm(*commDev);
  
  //Output
  idx_t edgecut;
  idx_t * part = new idx_t[mySize];
  
  //Parmetis call  
  int flag = ParMETIS_V3_PartMeshKway( elmdist, eptr, eind, elmwgt, 
    &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts, 
    tpwgts, &ubvec, options,
    &edgecut, part, 
    &mpiComm );
    
  assert(flag == METIS_OK);
  
  //Download the elements graph
  typedef geoElement<GEOSHAPE>              GEOELEMENT;
  typedef pGraph<GEOELEMENT,ELMAP,NODEMAP>  GRAPH;
  
  GRAPH elGraph = Grid1d.getElements();
    
  //Downloading the data on the graph map
  for(UInt i=1; i <= mySize; ++i)
  {
    elGraph.getRowMapL(i).setPid(part[i-1]);
  }
  
  //Clearing
  delete elmdist;
  delete eptr;
  delete eind;
  delete tpwgts;
  delete options;
  delete part;
  
  //COMMUNICATION ACROSS THE COMUNICATOR___________________________________________________________
  pVectComm<GEOELEMENT,ELMAP> graphTransfer(commDev);
  graphTransfer.vectorPid(elGraph);
  
  //Elements - Row to local
  pVectManip<GEOELEMENT,ELMAP> graphManipulator;
  graphManipulator.setNormalIndexing(elGraph);
  
  //Elements - Col to local mapping
  typedef typename map<UInt, UInt>::iterator ITER;
  
  pair<UInt, UInt> glPair;
  map<UInt, UInt>  glMapping;
  
  glPair.second = 0;  
  
  for(UInt i=1; i <= elGraph.rowSize(); ++i)
  {
    for(UInt j=1; j <= elGraph.rowSizeL(i); ++j)
    {
      glPair.first  = elGraph.getItemL(i).getCid(j);
      glMapping.insert(glPair);
    }
  }
  
  k = 1;
  for(ITER iter = glMapping.begin(); iter != glMapping.end(); ++iter)
  {
    iter->second = k;
    ++k;
  }
  
  //Insert the new mapping
  UInt gid;
  ITER iterator;
  
  elGraph.colIsLocal() = true;
  
  for(UInt i=1; i <= elGraph.rowSize(); ++i)
  {
    for(UInt j=1; j <= elGraph.rowSizeL(i); ++j)
    {
      gid                    = elGraph.getCid_LL(i,j);
      iterator               = glMapping.find(gid);
      elGraph.getCid_LL(i,j) = iterator->second;
    }
  }
  
  //Insert the new colMap
  NODEMAP mapItem;
  pMap<NODEMAP> colMap;
  
  for(ITER iter = glMapping.begin(); iter != glMapping.end(); ++iter)
  {
    mapItem.setLid(iter->second);
    mapItem.setGid(iter->first);
    mapItem.setPid(commDev->rank());
    
    colMap.push_back(mapItem);
  }
  
  //ColMapFixing
  colMapFixer_changeMap(colMap,*commDev);
  
  //Col map-pre cheking
  pMapGlobalManip<NODEMAP> hangingNodesCheck(commDev);
  assert(hangingNodesCheck.check(colMap));
  
  //Inserting in the graph
  elGraph.setColMap(colMap);
  elGraph.updateRowFinder();
  elGraph.updateColFinder();
  
  //Insert the new element-graph in the mesh
  Grid1d.setElements(elGraph);

  //Change the map to the nodes
  pVect<point3d,NODEMAP> nodesVect = Grid1d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Insert the nodes
  Grid1d.setNodes(nodesVect);
  
  //Internal methods
  Grid1d.computeNumVertices();
  Grid1d.transferMap();
  Grid1d.setMeshStandard(STDB);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshBalancing(const Teuchos::RCP<MESH1D> & Grid1d)  const
{
  assert(commDevLoaded);
  meshBalancing(*Grid1d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshBalancing(MESH1D & Grid1d) const
{
  assert(commDevLoaded);
  assert(Grid1d.getElements().colIsLocal());
  assert(Grid1d.getMeshStandard() == STDB);
  
  //Clearing
  Grid1d.transferMap();
  
  //PARMETIS_______________________________________________________________________________________
  UInt commSize = commDev->size();
  
  //Sizes communication
  UInt mySize = Grid1d.getNumElements();
  sVect<UInt> allSizes(commSize);
  
  all_gather(*commDev, mySize, allSizes);
  
  //Elements distribution
  idx_t * elmdist = new idx_t[commSize + 1];
  
  elmdist[0] = 0;
  
  for(UInt i=1 ; i <= commSize; ++i)
  {
    elmdist[i] = elmdist[i-1] + allSizes(i);
  }  
  
  //Element pointer
  idx_t * eptr = new idx_t[mySize + 1];
  
  eptr[0] = 0;
  
  for(UInt i=1 ; i <= mySize; ++i)
  {
    eptr[i] = eptr[i-1] + GEOSHAPE::numPoints;
  }
  
  //Element indices
  UInt  node, k = 0;
  idx_t * eind = new idx_t[mySize * GEOSHAPE::numPoints];
  
  for(UInt i=1; i <= mySize; ++i)
  {
    for(UInt j=1; j <= GEOSHAPE::numPoints; ++j)
    {
      node    = Grid1d.getElements().getCid_LG(i,j) - 1; //Change numbering      
      eind[k] = node;
      k++;
    }
    
    assert(k == UInt(eptr[i]));
  }
  
  //Some standard flags
  idx_t * elmwgt = NULL;
  idx_t wgtflag = 0;
  idx_t numflag = 0;
  idx_t ncon    = 1;
  
  //Other flags
  idx_t ncommonnodes = 1;
  idx_t nparts       = commSize;
  
  //Weights
  real_t * tpwgts = new real_t[ncon * nparts];
  
  for(int i=0; i < nparts; ++i)
  {
    tpwgts[i] = Real(1.0) / Real(nparts);
  }
  
  //Balancing index 
  real_t ubvec = 1.05;
  
  //Options Vector
  idx_t * options = new idx_t[3];
  options[0] = 0;
  options[1] = 0;
  options[2] = 0;
  
  //Mpi comm 
  MPI_Comm mpiComm(*commDev);
  
  //Output
  idx_t edgecut;
  idx_t * part = new idx_t[mySize];
  
  //Parmetis call  
  int flag = ParMETIS_V3_PartMeshKway( elmdist, eptr, eind, elmwgt, 
    &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts, 
    tpwgts, &ubvec, options,
    &edgecut, part, 
    &mpiComm );
    
  assert(flag == METIS_OK);
  
  //Download the elements graph
  typedef geoElement<GEOSHAPE>              GEOELEMENT;
  typedef pGraph<GEOELEMENT,ELMAP,NODEMAP>  GRAPH;
  
  GRAPH elGraph = Grid1d.getElements();
    
  //Downloading the data on the graph map
  for(UInt i=1; i <= mySize; ++i)
  {
    elGraph.getRowMapL(i).setPid(part[i-1]);
  }
  
  //Clearing
  delete elmdist;
  delete eptr;
  delete eind;
  delete tpwgts;
  delete options;
  delete part;
  
  
  //COMMUNICATION ACROSS THE COMUNICATOR___________________________________________________________
  elGraph.pushToGlobal();
  
  //Graph transfer
  pVectComm<GEOELEMENT,ELMAP> graphTransfer(commDev);
  graphTransfer.vectorPid(elGraph);
  
  //Elements - Row to local
  pVectManip<GEOELEMENT,ELMAP> graphManipulator;
  graphManipulator.setNormalIndexing(elGraph);
  
  //Elements - Col to local mapping
  typedef typename map<UInt, UInt>::iterator ITER;
  
  pair<UInt, UInt> glPair;
  map<UInt, UInt>  glMapping;
  
  glPair.second = 0;
  
  for(UInt i=1; i <= elGraph.rowSize(); ++i)
  {
    for(UInt j=1; j <= elGraph.rowSizeL(i); ++j)
    {
      glPair.first  = elGraph.getItemL(i).getCid(j);
      glMapping.insert(glPair);
    }
  }
  
  k = 1;
  for(ITER iter = glMapping.begin(); iter != glMapping.end(); ++iter)
  {
    iter->second = k;
    ++k;
  }
  
  //Insert the new mapping
  UInt gid;
  ITER iterator;
  
  elGraph.colIsLocal() = true;
  
  for(UInt i=1; i <= elGraph.rowSize(); ++i)
  {
    for(UInt j=1; j <= elGraph.rowSizeL(i); ++j)
    {
      gid                    = elGraph.getCid_LL(i,j);
      iterator               = glMapping.find(gid);
      elGraph.getCid_LL(i,j) = iterator->second;
    }
  }
  
  //Insert the new colMap
  NODEMAP mapItem;
  pMap<NODEMAP> colMap;
  
  for(ITER iter = glMapping.begin(); iter != glMapping.end(); ++iter)
  {
    mapItem.setLid(iter->second);
    mapItem.setGid(iter->first);
    mapItem.setPid(commDev->rank());
    
    colMap.push_back(mapItem);
  }
  
  //ColMapFixing
  colMapFixer_changeMap(colMap,*commDev);
  
  //Inserting in the graph
  elGraph.setColMap(colMap);
  elGraph.updateRowFinder();
  elGraph.updateColFinder();
  
  //Insert the new element-graph in the mesh
  Grid1d.setElements(elGraph);

  //Change the map to the nodes
  pVect<point3d,NODEMAP> nodesVect = Grid1d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Insert the nodes
  Grid1d.setNodes(nodesVect);
  
  //Internal methods
  Grid1d.computeNumVertices();
  Grid1d.transferMap();
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshOverlapping(const Teuchos::RCP<MESH1D> & Grid1d) const
{
  assert(commDevLoaded);
  meshOverlapping(*Grid1d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshOverlapping(MESH1D & Grid1d) const
{
  assert(commDevLoaded);
  assert(Grid1d.getElements().colIsLocal());
  assert(Grid1d.getMeshStandard() == STDB);
  
  //Clearing
  Grid1d.transferMap();
  
  //Elements overlapping
  typedef pGraph<GEOELEMENT,ELMAP,NODEMAP>  GRAPH;
  
  GRAPH elGraph = Grid1d.getElements();
  
  pGraphGlobalManip<GEOELEMENT,pMapItem,NODEMAP> elGraphManip(commDev);
  elGraphManip.overlap(elGraph);
  
  //Nodes overlapping
  pMap<NODEMAP>          colMap    = elGraph.getColMap();
  pVect<point3d,NODEMAP> nodesVect = Grid1d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Grid fill
  Grid1d.setElements(elGraph);
  Grid1d.setNodes(nodesVect);
  Grid1d.transferMap();
  Grid1d.setMeshStandard(STDA);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
gather(const UInt & gatherGid, const Teuchos::RCP<MESH1D> & targetGrid1d, const Teuchos::RCP<const MESH1D> & sourceGrid1d)
{
  assert(commDevLoaded);
  gather(gatherGid,*targetGrid1d,*sourceGrid1d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
gather(const UInt & gatherGid, MESH1D & targetGrid1d, const MESH1D & sourceGrid1d)
{
  typedef typename MESH1D::NODESVECT    NODESVECT;
  typedef typename MESH1D::GRAPH1D      GRAPH1D;
  typedef typename MESH1D::GEOELEMENT1D GEOELEMENT1D;
  
  //Mesh cheking
  assert(commDevLoaded);
  assert(sourceGrid1d.getElements().colIsLocal());
  
  //Nodes and graph extraction
  NODESVECT nodes    = sourceGrid1d.getNodes();
  GRAPH1D   elements = sourceGrid1d.getElements();
  
  //Pids renumbering
  for(UInt i=1; i <= nodes.size(); ++i)
  {
    nodes.getMapL(i).setPid(gatherGid);
  }
  
  for(UInt i=1; i <= elements.size(); ++i)
  {
    elements.getRowMapL(i).setPid(gatherGid);
  }
  
  
  //Pushing to global
  elements.pushToGlobal();
  
  
  //Nodes - Elements communication
  pVectComm<point2d,NODEMAP> nodesComm(commDev);
  nodesComm.vectorPid(nodes);

  pVectComm<GEOELEMENT1D,ELMAP> elementsComm(commDev);
  elementsComm.vectorPid(elements);
  
  
  //Nodes - Elements exclusive merging
  pVectManip<point2d,NODEMAP> nodesManip;
  nodesManip.unionExclusive(nodes);  
  
  pVectManip<GEOELEMENT1D,ELMAP> elementsManip;
  elementsManip.unionExclusive(elements);
  
  
  //Set notmal indexing
  pMapManip<NODEMAP> nodesMapManip;
  nodesMapManip.setNormalIndexing(nodes.getMapRef());
  
  pMapManip<ELMAP> elementsMapManip;
  elementsMapManip.setNormalIndexing(elements.getRowMap());
  
  //Insert
  elements.setColMap(nodes.getMapRef());
  elements.updateRowFinder();
  elements.updateColFinder();
  elements.colIsLocal() = true;
  
  targetGrid1d.clear();
  targetGrid1d.setNodes(nodes);
  targetGrid1d.setElements(elements);
  targetGrid1d.transferMap();
  targetGrid1d.setMeshStandard(STDU);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
destroyOverlap(const Teuchos::RCP<MESH1D> & grid1d) const
{
  assert(commDevLoaded);
  destroyOverlap(*grid1d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
destroyOverlap(MESH1D & grid1d) const
{
  //Mesh cheking-----------------------------------------------------------------------------------
  assert(commDevLoaded);
  assert(grid1d.getElements().colIsLocal());
  
  //Typedefs
  typedef typename MESH1D::GRAPH1D    GRAPH1D;
  typedef typename MESH1D::NODESVECT  NODESVECT;
  
  //Extract data
  GRAPH1D elements = grid1d.getElements();
  NODESVECT  nodes = grid1d.getNodes();
  
  
  //Find shared elements---------------------------------------------------------------------------
  pMapItemShare sharedItem;
  pMap<pMapItemShare> tempMap(elements.size());
  
  for(UInt i=1; i <= elements.size(); ++i)
  {
    sharedItem.setPid(elements.getMapL(i).getPid());
    sharedItem.setLid(elements.getMapL(i).getLid());
    sharedItem.setGid(elements.getMapL(i).getGid());
    
    tempMap(i) = sharedItem;
  }
  
  pMapGlobalManip<pMapItemShare> sharedMapManip(commDev);
  sharedMapManip.updateOwningSharing(tempMap);
  
  //Eliminate the unused elements
  UInt tot = 1;
  GRAPH1D  newElements;
  pMapItem elMapItem;
  
  for(UInt i=1; i <= elements.size(); ++i)
  {
    if(tempMap(i).getOwned())
    {
      elMapItem = elements.getRowMapL(i);
      elMapItem.setLid(tot);
      newElements.push_back(elMapItem,elements.getItemL(i));
      
      ++tot;
    }
  }

  //Eliminate the unused points--------------------------------------------------------------------
  UInt id;
  
  sVect<bool> active(nodes.size());
  sVect<UInt> newOrder(nodes.size());
  
  for(UInt i=1; i <= active.size(); ++i)
  { active(i) = false; }
  
  //Find active nodes
  for(UInt i=1; i <= newElements.size(); ++i)
  {
    for(UInt j=1; j <= grid1d.getNumPoints(); ++j)  //I points sono locali
    {
      id         = newElements.getCid_LL(i,j);
      active(id) = true;
    }
  }
  
  //Nodes renumbering and reduction
  tot = 1;
  NODESVECT newNodes;
  NODEMAP   nodeMap;
  
  for(UInt i=1; i <= active.size(); ++i)
  {
    if(active(i))
    {
      newOrder(i) = tot;
      
      nodeMap = nodes.getMapL(i);
      nodeMap.setLid(tot);
      
      newNodes.push_back(nodes.getDataL(i), nodeMap);
      tot++;
    }
  }
  
  newNodes.updateFinder();
  
  //Elements correction
  geoElement<GEOSHAPE> geoElement1d(true);
  
  for(UInt i=1; i <= newElements.size(); ++i)
  {
    geoElement1d = newElements.getItemL(i);
    
    for(UInt j=1; j <= grid1d.getNumPoints(); ++j)
    {
      id = newElements.getCid_LL(i,j);
      geoElement1d.setCid(j,newOrder(id));
    }
    
    newElements(i) = geoElement1d;
  } 
  
  
  //Rebuild the owning-sharing for the nodes-------------------------------------------------------
  nodesMapFixer(newNodes.getMapRef(),*commDev);
  
  //Build new grid---------------------------------------------------------------------------------
  newNodes.updateFinder();
  
  newElements.setColMap(newNodes.getMapRef());
  newElements.updateRowFinder();
  newElements.updateColFinder();
  
  grid1d.clearEdges();
  
  grid1d.setNodes(newNodes);
  grid1d.setElements(newElements);
  grid1d.computeNumVertices();
  grid1d.transferMap();
  grid1d.setMeshStandard(STDB);
  
  //Final check------------------------------------------------------------------------------------
  assert(check(grid1d));
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MESH1D>       & OldMesh,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MESH1D>             & NewMesh)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldMesh,
                    *NewCommDev,
                    *NewMesh);
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
reduceCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MESH1D       & OldMesh,
                   const communicator & NewCommDev,
                         MESH1D       & NewMesh)
{
  //Typedefs---------------------------------------------------------
  typedef typename MESH1D::NODESVECT NODESVECT;
  typedef typename MESH1D::BOOLVECT  BOOLVECT;
  typedef typename MESH1D::GRAPH1D   GRAPH1D;
  
  typedef pVectGlobalManip<point3d,NODEMAP>                MANIP_POINT;
  typedef pVectGlobalManip<UInt,NODEMAP>                   MANIP_ISVERTEX;
  typedef pGraphGlobalManip<GEOELEMENT,pMapItem,NODEMAP>   MANIP_ELEMENT;
  
  //Manipulators-----------------------------------------------------
  NODESVECT newNodes;
  BOOLVECT  newIsVertex;
  GRAPH1D   newElements;
  
  MANIP_POINT    manipPoint;
  MANIP_ISVERTEX manipIsVertex;
  MANIP_ELEMENT  manipElement;
  
  //Copy to new communicator-----------------------------------------
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() < OldCommDev.size());
    
    //Alloc
    NewMesh.clear();
    NewMesh = OldMesh;
    
    //Reduce comm
    manipPoint.reduceCommunicator(isActive,
                                  OldCommDev,
                                  OldMesh.getNodes(),
                                  NewCommDev,
                                  newNodes);
    
    manipIsVertex.reduceCommunicator(isActive,
                                     OldCommDev,
                                     OldMesh.getIsVertex(),
                                     NewCommDev,
                                     newIsVertex);
    
    manipElement.reduceCommunicator(isActive,
                                    OldCommDev,
                                    OldMesh.getElements(),
                                    NewCommDev,
                                    newElements);
    
    //Fix the mesh
    NewMesh.setNodes(newNodes,newIsVertex);
    NewMesh.setElements(newElements);
    
    NewMesh.computeNumVertices();
    NewMesh.transferMap();
  }
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MESH1D>       & OldMesh,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MESH1D>             & NewMesh)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldMesh,
                    *NewCommDev,
                    *NewMesh);
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
expandCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MESH1D       & OldMesh,
                   const communicator & NewCommDev,
                         MESH1D       & NewMesh)
{
  //Typedefs---------------------------------------------------------
  typedef typename MESH1D::NODESVECT NODESVECT;
  typedef typename MESH1D::BOOLVECT  BOOLVECT;
  typedef typename MESH1D::GRAPH1D   GRAPH1D;
  
  typedef pVectGlobalManip<point3d,NODEMAP>                MANIP_POINT;
  typedef pVectGlobalManip<UInt,NODEMAP>                   MANIP_ISVERTEX;
  typedef pGraphGlobalManip<GEOELEMENT,pMapItem,NODEMAP>   MANIP_ELEMENT;
  
  //Manipulators-----------------------------------------------------
  NODESVECT newNodes;
  BOOLVECT  newIsVertex;
  GRAPH1D   newElements;
  
  MANIP_POINT    manipPoint;
  MANIP_ISVERTEX manipIsVertex;
  MANIP_ELEMENT  manipElement;
  
  //Copy to new communicator-----------------------------------------
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() > OldCommDev.size());
    
    //Alloc
    NewMesh.clear();
    NewMesh = OldMesh;
    
    //Reduce comm
    manipPoint.expandCommunicator(isActive,
                                  OldCommDev,
                                  OldMesh.getNodes(),
                                  NewCommDev,
                                  newNodes);
    
    manipIsVertex.expandCommunicator(isActive,
                                     OldCommDev,
                                     OldMesh.getIsVertex(),
                                     NewCommDev,
                                     newIsVertex);
    
    manipElement.expandCommunicator(isActive,
                                    OldCommDev,
                                    OldMesh.getElements(),
                                    NewCommDev,
                                    newElements);
    
    //Fix the mesh
    NewMesh.setNodes(newNodes,newIsVertex);
    NewMesh.setElements(newElements);
  }
  
  NewMesh.computeNumVertices();
  NewMesh.transferMap();
}



//_________________________________________________________________________________________________
// PMAPITEMSHARE
//-------------------------------------------------------------------------------------------------

/*! Global mesh manipulator 1d - \c pMapItemShare */
template<typename GEOSHAPE, typename NODEMAP>
class mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItemShare ELMAP;
    typedef mesh1d<GEOSHAPE,ELMAP,NODEMAP> MESH1D;
    
    typedef typename GEOSHAPE::GEOBSHAPE   GEOSHAPE0D;
    
    typedef geoElement<GEOSHAPE>   GEOELEMENT;
    typedef geoElement<GEOSHAPE0D> GEOELEMENT0D;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
  
    /*! @name Constructors */ //@{
  public:
    mesh1dGlobalManip();
    mesh1dGlobalManip(const Teuchos::RCP<const communicator> & CommDev);
    mesh1dGlobalManip(const communicator & CommDev);
    //@}
    
    /*! @name Information functions */ //@{
  public:
    /*! Check the graphs, the equality of the col maps of the elements, faces and edges with the map of the points */
    bool check(const Teuchos::RCP<const MESH1D> & Grid1d) const;
    bool check(const MESH1D & Grid1d) const;
    
    /*! Return the total number of vertices on the processors */
    UInt getNumGlobalVertices(const Teuchos::RCP<const MESH1D> & Grid1d) const;
    UInt getNumGlobalVertices(const MESH1D & Grid1d) const;
    
    /*! Return the total numer of elements on the processors */
    UInt getNumGlobalElements(const Teuchos::RCP<const MESH1D> & Grid1d) const;
    UInt getNumGlobalElements(const MESH1D & Grid1d) const;
    //@}
    
    /*! @name Operation functions */ //@{
  public:
    /*! The mesh should be in global numbering and a new mesh in local numbering is built */
    void meshPartition(const Teuchos::RCP<MESH1D> & Grid1d) const;
    void meshPartition(MESH1D & Grid1d) const;
    
    /*! The mesh should be in local numbering, the not-owned elements are discarded */
    void meshBalancing(const Teuchos::RCP<MESH1D> & Grid1d) const;
    void meshBalancing(MESH1D & grid1d) const;
    
    /*! The mesh should be in local numbering. The appendended elements and nodes are listed below the original ones */
    void meshOverlapping(const Teuchos::RCP<MESH1D> & Grid1d) const;
    void meshOverlapping(MESH1D & Grid1d) const;
    
    /*! Reduces the mesh to a single process and pushes the \c cids to the global numbering.
    Duplicated elements or points are eliminated. */
    void gather(const UInt & gatherGid, const Teuchos::RCP<MESH1D> & targetGrid1d, const Teuchos::RCP<const MESH1D> & sourceGrid1d);
    void gather(const UInt & gatherGid, MESH1D & targetGrid1d, const MESH1D & sourceGrid1d);
    
    /*! Eliminates the overlapped regions. ATTENTION: since there is no way history of the way
    the mesh was overlapped the shared elements are put in a particular domain in an
    un-predictable manner */
    void destroyOverlap(const Teuchos::RCP<MESH1D> & grid1d) const;
    void destroyOverlap(MESH1D & grid1d) const;
    //@}
    
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MESH1D>       & OldMesh,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MESH1D>             & NewMesh);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MESH1D       & OldMesh,
                            const communicator & NewCommDev,
                                  MESH1D       & NewMesh);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MESH1D>       & OldMesh,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MESH1D>             & NewMesh);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MESH1D       & OldMesh,
                            const communicator & NewCommDev,
                                  MESH1D       & NewMesh);
    //@}
};

template<typename GEOSHAPE, typename NODEMAP>
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
mesh1dGlobalManip()
{
  assert(GEOSHAPE::nDim == 1);
  commDevLoaded = false;
}

template<typename GEOSHAPE, typename NODEMAP>
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
mesh1dGlobalManip(const Teuchos::RCP<const communicator> & CommDev)
{
  assert(GEOSHAPE::nDim == 1);
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename NODEMAP>
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
mesh1dGlobalManip(const communicator & CommDev)
{
  assert(GEOSHAPE::nDim == 1);
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename NODEMAP>
bool
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
check(const Teuchos::RCP<const MESH1D> & Grid1d) const
{
  return(check(*Grid1d));
}

template<typename GEOSHAPE, typename NODEMAP>
bool
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
check(const MESH1D & Grid1d) const
{
  assert(Grid1d.getMapTransferred());
  
  //Allocations
  bool flag;
  UInt pid  = commDev->rank();
  UInt size = commDev->size();
  
  //Nodes vect check
  pVectGlobalManip<point3d,NODEMAP> nodesChecking(commDev);
  bool flagN = nodesChecking.check(Grid1d.getNodes());
  
  if(!flagN)
  { cout << "ERROR: Nodes check failed. Pid: " << pid << endl; }
  
  //Elements check
  pGraphGlobalManip<GEOELEMENT,pMapItemShare,NODEMAP> elementsChecking(commDev);
  bool flagE = elementsChecking.check(Grid1d.getElements());
  
  if(!flagE)
  { cout << "ERROR: Elements check failed. Pid: " << pid << endl; }
  
  flag = flagN & flagE;
  
  //Cross map checking
  pMapManip<NODEMAP> crossMapChecker;
  
  //Checking against elements col map
  sVect<UInt> allFlags(size);
  UInt myFlag = crossMapChecker.isEqual(Grid1d.getNodes().getMapRef(),Grid1d.getElements().getColMap());
  
  if(myFlag != 1)
  { cout << "ERROR: Elements col map differs from Nodes col map. Pid: " << pid << endl; }
  
  all_gather(*commDev, myFlag, allFlags);
  
  for(UInt i=1; i <= size; ++i)
  { myFlag = myFlag & allFlags(i); }
  
  flag = flag & myFlag;
  
  
  //Cross checking for ownership - elements
  auxiliaryGraphCheck1d<GEOSHAPE,ELMAP,NODEMAP> elOwCheker;
  bool flagWE = elOwCheker.check(Grid1d.getElements());
  
  flag = flag && flagWE;
  
  if(!flagWE)
  { cout << "ERROR: Elements cross ownership failed. Pid: " << pid << endl; }  
  
  
  //Cross check among processors
  sVect<UInt> flags;
  all_gather(*commDev,UInt(flag),flags);
  
  for(UInt i=1; i <= size; ++i)
  { flag = flag & bool(flags(i));}
  
  return(flag);
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalVertices(const Teuchos::RCP<const MESH1D> & Grid1d) const
{
  pVectGlobalManip<point3d,NODEMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid1d->getNodes()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalVertices(const MESH1D & Grid1d) const
{
  pVectGlobalManip<point3d,NODEMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid1d.getNodes()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalElements(const Teuchos::RCP<const MESH1D> & Grid1d) const
{
  pVectGlobalManip<GEOELEMENT,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid1d->getElements()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalElements(const MESH1D & Grid1d) const
{
  pVectGlobalManip<GEOELEMENT,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid1d.getElements()));
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshPartition(const Teuchos::RCP<MESH1D> & Grid1d) const
{
  assert(commDevLoaded);
  meshPartition(*Grid1d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshPartition(MESH1D & Grid1d) const
{
  //Asserts
  assert(commDevLoaded);
  assert(Grid1d.getElements().colIsLocal() == false);
  assert(Grid1d.getMeshStandard() == STDL);
  
  //Clearing
  Grid1d.transferMap();
  
  //Sharing-owningMap______________________________________________________________________________
  typedef typename MESH1D::GRAPH1D GRAPH1D;
  
  GRAPH1D tempElements = Grid1d.getElements();
  
  for(UInt i=1; i <= tempElements.getRowMap().size(); ++i)
  {
    tempElements.getRowMap().get(i).setShared(false);
    tempElements.getRowMap().get(i).setOwned(true);
  }
  
  Grid1d.setElements(tempElements);
  Grid1d.transferMap();
  
  //PARMETIS_______________________________________________________________________________________
  UInt commSize = commDev->size();
  
  //Sizes communication
  UInt mySize = Grid1d.getNumElements();
  sVect<UInt> allSizes(commSize);
  
  all_gather(*commDev, mySize, allSizes);
  
  //Elements distribution
  idx_t * elmdist = new idx_t[commSize + 1];
  
  elmdist[0] = 0;
  
  for(UInt i=1 ; i <= commSize; ++i)
  {
    elmdist[i] = elmdist[i-1] + allSizes(i);
  }  
  
  //Element pointer
  idx_t * eptr = new idx_t[mySize + 1];
  
  eptr[0] = 0;
  
  for(UInt i=1 ; i <= mySize; ++i)
  {
    eptr[i] = eptr[i-1] + GEOSHAPE::numPoints;
  }
   
  //Element indices
  UInt  node, k = 0;
  idx_t * eind = new idx_t[mySize * GEOSHAPE::numPoints];
  
  for(UInt i=1; i <= mySize; ++i)
  {
    for(UInt j=1; j <= GEOSHAPE::numPoints; ++j)
    {
      node    = Grid1d.getElementL(i).getCid(j) - 1; //Change numbering      
      eind[k] = node;
      k++;
    }
    
    assert(k == UInt(eptr[i]));
  }
  
  //Some tandard flags
  idx_t * elmwgt = NULL;
  idx_t wgtflag = 0;
  idx_t numflag = 0;
  idx_t ncon    = 1;
  
  //Other flags
  idx_t ncommonnodes = 1;
  idx_t nparts       = commSize;
  
  //Weights
  real_t * tpwgts = new real_t[ncon * nparts];
  
  for(int i=0; i < nparts; ++i)
  {
    tpwgts[i] = Real(1.0) / Real(nparts);
  }
  
  //Balancing index 
  real_t ubvec = 1.05;
  
  //Options Vector
  idx_t * options = new idx_t[3];
  options[0] = 0;
  options[1] = 0;
  options[2] = 0;
  
  //Mpi comm 
  MPI_Comm mpiComm(*commDev);
  
  //Output
  idx_t edgecut;
  idx_t * part = new idx_t[mySize];
  
  //Parmetis call  
  ParMETIS_V3_PartMeshKway( elmdist, eptr, eind, elmwgt, 
    &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts, 
    tpwgts, &ubvec, options,
    &edgecut, part, 
    &mpiComm );
  
  //Download the elements graph
  typedef geoElement<GEOSHAPE>              GEOELEMENT;
  typedef pGraph<GEOELEMENT,ELMAP,NODEMAP>  GRAPH;
  
  GRAPH elGraph = Grid1d.getElements();
    
  //Downloading the data on the graph map
  for(UInt i=1; i <= mySize; ++i)
  {
    elGraph.getRowMapL(i).setPid(part[i-1]);
  }
  
  //Clearing
  delete elmdist;
  delete eptr;
  delete eind;
  delete tpwgts;
  delete options;
  delete part;
  
  //COMMUNICATION ACROSS THE COMUNICATOR___________________________________________________________
  pVectComm<GEOELEMENT,ELMAP> graphTransfer(commDev);
  graphTransfer.vectorPid(elGraph);
  
  //Elements - Row to local
  pVectManip<GEOELEMENT,ELMAP> graphManipulator;
  graphManipulator.setNormalIndexing(elGraph);
  
  //Elements - Col to local mapping
  typedef typename map<UInt, UInt>::iterator ITER;
  
  pair<UInt, UInt> glPair;
  map<UInt, UInt>  glMapping;
  
  glPair.second = 0;
  
  for(UInt i=1; i <= elGraph.rowSize(); ++i)
  {
    for(UInt j=1; j <= elGraph.rowSizeL(i); ++j)
    {
      glPair.first  = elGraph.getItemL(i).getCid(j);
      glMapping.insert(glPair);
    }
  }
  
  k = 1;
  for(ITER iter = glMapping.begin(); iter != glMapping.end(); ++iter)
  {
    iter->second = k;
    ++k;
  }
  
  //Insert the new mapping
  UInt gid;
  ITER iterator;
  
  elGraph.colIsLocal() = true;
  
  for(UInt i=1; i <= elGraph.rowSize(); ++i)
  {
    for(UInt j=1; j <= elGraph.rowSizeL(i); ++j)
    {
      gid                    = elGraph.getCid_LL(i,j);
      iterator               = glMapping.find(gid);
      elGraph.getCid_LL(i,j) = iterator->second;
    }
  }
  
  //Insert the new colMap
  NODEMAP mapItem;
  pMap<NODEMAP> colMap;
  
  for(ITER iter = glMapping.begin(); iter != glMapping.end(); ++iter)
  {
    mapItem.setLid(iter->second);
    mapItem.setGid(iter->first);
    mapItem.setPid(commDev->rank());
    
    colMap.push_back(mapItem);
  }
  
  //ColMapFixing
  colMapFixer_changeMap(colMap,*commDev);
  
  //Col map-pre cheking
  pMapGlobalManip<NODEMAP> hangingNodesCheck(commDev);
  assert(hangingNodesCheck.check(colMap));
  
  //Inserting in the graph
  elGraph.setColMap(colMap);
  elGraph.updateRowFinder();
  elGraph.updateColFinder();
  
  //Insert the new element-graph in the mesh
  Grid1d.setElements(elGraph);

  //Change the map to the nodes
  pVect<point3d,NODEMAP> nodesVect = Grid1d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Insert the nodes
  Grid1d.setNodes(nodesVect);
  
  //Internal methods
  Grid1d.computeNumVertices();
  Grid1d.transferMap();
  Grid1d.setMeshStandard(STDB);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshBalancing(const Teuchos::RCP<MESH1D> & Grid1d) const
{
  assert(commDevLoaded);
  meshBalancing(*Grid1d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshBalancing(MESH1D & Grid1d) const
{
  assert(commDevLoaded);
  assert(Grid1d.getMeshStandard() == STDB);
  
  typedef typename MESH1D::GRAPH1D      ELGRAPH;
  typedef typename MESH1D::GEOELEMENT1D ELEMENT;
  typedef typename MESH1D::NODESVECT    NODESVECT;
  
  //Clearing
  Grid1d.transferMap();
  
  //Grid copy
  ELGRAPH   elementsCopy = Grid1d.getElements();
  ELGRAPH   elList;
  NODESVECT nodesCopy = Grid1d.getNodes();;
  NODESVECT nodesList;
  
  ELEMENT       element;
  point3d       node;
  pMapItemShare elMapItem;
  NODEMAP       nodeMapItem;
  
  elementsCopy.pushToGlobal();
  
  //Shared elements cut-out
  UInt k=1;
  for(UInt i=1; i <= elementsCopy.getRowMap().size(); ++i)
  {
    elMapItem = elementsCopy.getRowMapL(i);
    
    if(elMapItem.getOwned())
    {
      element = elementsCopy.getItemL(i);
      elMapItem.setLid(k);
      elMapItem.setShared(false);
      elMapItem.setOwned(true);
      k++;
      
      elList.push_back(element,elMapItem);
    }
  }
  
  //Shared nodes cut-out
  k=1;
  for(UInt i=1; i <= nodesCopy.getMapRef().size(); ++i)
  {
    nodeMapItem = nodesCopy.getMapL(i);
    
    if(nodeMapItem.getOwned())
    {
      node = nodesCopy.getDataL(i);
      nodeMapItem.setLid(k);
      
      nodesList.push_back(node,nodeMapItem);
    }
  }
  
  //Spoiled mesh finalization
  MESH1D tempGrid;
  
  elList.colIsLocal() = false;
  nodesList.updateFinder();
  
  tempGrid.setElements(elList);
  tempGrid.setNodes(nodesList);
  tempGrid.transferMap();
  
  //Mesh repartition
  tempGrid.setMeshStandard(STDL);
  meshPartition(tempGrid);
  tempGrid.setMeshStandard(STDB);
  
  //Mesh loading
  Grid1d = tempGrid;
  Grid1d.transferMap();
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshOverlapping(const Teuchos::RCP<MESH1D> & Grid1d) const
{
  assert(commDevLoaded);
  meshOverlapping(*Grid1d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshOverlapping(MESH1D & Grid1d) const
{
  assert(commDevLoaded);
  assert(Grid1d.getElements().colIsLocal());
  assert(Grid1d.getMeshStandard() == STDB);
  
  //Clearing
  Grid1d.transferMap();
  
  //Elements overlapping
  typedef pGraph<GEOELEMENT,ELMAP,NODEMAP>  GRAPH;
  
  GRAPH elGraph = Grid1d.getElements();
  
  pGraphGlobalManip<GEOELEMENT,pMapItemShare,NODEMAP> elGraphManip(commDev);
  elGraphManip.overlap(elGraph);
  
  //Nodes overlapping
  pMap<NODEMAP>          colMap    = elGraph.getColMap();
  pVect<point3d,NODEMAP> nodesVect = Grid1d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Grid fill
  Grid1d.setElements(elGraph);
  Grid1d.setNodes(nodesVect);
  Grid1d.transferMap();
  Grid1d.setMeshStandard(STDA);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
gather(const UInt & gatherGid, const Teuchos::RCP<MESH1D> & targetGrid1d, const Teuchos::RCP<const MESH1D> & sourceGrid1d)
{
  assert(commDevLoaded);
  gather(gatherGid,*targetGrid1d,*sourceGrid1d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
gather(const UInt & gatherGid, MESH1D & targetGrid1d, const MESH1D & sourceGrid1d)
{
  typedef typename MESH1D::NODESVECT    NODESVECT;
  typedef typename MESH1D::GRAPH1D      GRAPH1D;
  typedef typename MESH1D::GEOELEMENT1D GEOELEMENT1D;
  
  //Mesh cheking
  assert(commDevLoaded);
  assert(sourceGrid1d.getElements().colIsLocal());
  
  //Nodes and graph extraction
  NODESVECT nodes    = sourceGrid1d.getNodes();
  GRAPH1D   elements = sourceGrid1d.getElements();
  
  //Pids renumbering
  for(UInt i=1; i <= nodes.size(); ++i)
  {
    nodes.getMapL(i).setPid(gatherGid);
    nodes.getMapL(i).setOwned(true);
    nodes.getMapL(i).setShared(false);
  }
  
  for(UInt i=1; i <= elements.size(); ++i)
  {
    elements.getRowMapL(i).setPid(gatherGid);
    elements.getRowMapL(i).setOwned(true);
    elements.getRowMapL(i).setShared(false);
  }
  
  
  //Pushing to global
  elements.pushToGlobal();
  
  
  //Nodes - Elements communication
  pVectComm<point3d,NODEMAP> nodesComm(commDev);
  nodesComm.vectorPid(nodes);

  pVectComm<GEOELEMENT1D,ELMAP> elementsComm(commDev);
  elementsComm.vectorPid(elements);
  
  
  //Nodes - Elements exclusive merging
  pVectManip<point3d,NODEMAP> nodesManip;
  nodesManip.unionExclusive(nodes);  
  
  pVectManip<GEOELEMENT1D,ELMAP> elementsManip;
  elementsManip.unionExclusive(elements);
  
  
  //Set notmal indexing
  pMapManip<NODEMAP> nodesMapManip;
  nodesMapManip.setNormalIndexing(nodes.getMapRef());
  
  pMapManip<ELMAP> elementsMapManip;
  elementsMapManip.setNormalIndexing(elements.getRowMap());
  
  //Insert
  elements.setColMap(nodes.getMapRef());
  elements.updateRowFinder();
  elements.updateColFinder();
  elements.colIsLocal() = true;
  
  targetGrid1d.clear();
  targetGrid1d.setNodes(nodes);
  targetGrid1d.setElements(elements);
  targetGrid1d.transferMap();
  targetGrid1d.setMeshStandard(STDU);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
destroyOverlap(const Teuchos::RCP<MESH1D> & grid1d) const
{
  assert(commDevLoaded);
  destroyOverlap(*grid1d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
destroyOverlap(MESH1D & grid1d) const
{
  assert(grid1d->getMeshStandard() == STDA);
  
  //Mesh cheking-----------------------------------------------------------------------------------
  assert(commDevLoaded);
  assert(grid1d.getElements().colIsLocal());
  
  //Typedefs
  typedef typename MESH1D::GRAPH1D    GRAPH1D;
  typedef typename MESH1D::NODESVECT  NODESVECT;
  
  //Extract data
  GRAPH1D elements = grid1d.getElements();
  NODESVECT  nodes = grid1d.getNodes();
  
  
  //Eliminate the unused elements------------------------------------------------------------------
  UInt tot = 1;
  GRAPH1D       newElements;
  pMapItemShare elMapItem;
  
  for(UInt i=1; i <= elements.size(); ++i)
  {
    elMapItem = elements.getRowMapL(i);
    
    if(elMapItem.getOwned())
    {
      elMapItem.setLid(tot);
      elMapItem.setShared(false);
      elMapItem.setOwned(true);
      newElements.push_back(elMapItem,elements.getItemL(i));
      
      ++tot;
    }
  }

  //Eliminate the unused points--------------------------------------------------------------------
  UInt id;
  
  sVect<bool> active(nodes.size());
  sVect<UInt> newOrder(nodes.size());
  
  for(UInt i=1; i <= active.size(); ++i)
  { active(i) = false; }
  
  //Find active nodes
  for(UInt i=1; i <= newElements.size(); ++i)
  {
    for(UInt j=1; j <= grid1d.getNumPoints(); ++j)  //I points sono locali
    {
      id         = newElements.getCid_LL(i,j);
      active(id) = true;
    }
  }
  
  //Nodes renumbering and reduction
  tot = 1;
  NODESVECT newNodes;
  NODEMAP   nodeMap;
  
  for(UInt i=1; i <= active.size(); ++i)
  {
    if(active(i))
    {
      newOrder(i) = tot;
      
      nodeMap = nodes.getMapL(i);
      nodeMap.setLid(tot);
      
      newNodes.push_back(nodes.getDataL(i), nodeMap);
      tot++;
    }
  }
  
  newNodes.updateFinder();
  
  //Elements correction
  geoElement<GEOSHAPE> geoElement1d(true);
  
  for(UInt i=1; i <= newElements.size(); ++i)
  {
    geoElement1d = newElements.getItemL(i);
    
    for(UInt j=1; j <= grid1d.getNumPoints(); ++j)
    {
      id = newElements.getCid_LL(i,j);
      geoElement1d.setCid(j,newOrder(id));
    }
    
    newElements(i) = geoElement1d;
  } 
  
  
  //Rebuild the owning-sharing for the nodes-------------------------------------------------------
  nodesMapFixer(newNodes.getMapRef(),*commDev);
  
  //Build new grid---------------------------------------------------------------------------------
  newNodes.updateFinder();
  
  newElements.setColMap(newNodes.getMapRef());
  newElements.updateRowFinder();
  newElements.updateColFinder();
  
  grid1d.clearFaces();
  grid1d.clearEdges();
  
  grid1d.setNodes(newNodes);
  grid1d.setElements(newElements);
  grid1d.computeNumVertices();
  grid1d.transferMap();
  grid1d.setMeshStandard(STDB);
  
  //Final check------------------------------------------------------------------------------------
  assert(check(grid1d));
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MESH1D>       & OldMesh,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MESH1D>             & NewMesh)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldMesh,
                    *NewCommDev,
                    *NewMesh);
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
reduceCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MESH1D       & OldMesh,
                   const communicator & NewCommDev,
                         MESH1D       & NewMesh)
{
  //Typedefs---------------------------------------------------------
  typedef typename MESH1D::NODESVECT NODESVECT;
  typedef typename MESH1D::BOOLVECT  BOOLVECT;
  typedef typename MESH1D::GRAPH1D   GRAPH1D;
  
  typedef pVectGlobalManip<point3d,NODEMAP>                     MANIP_POINT;
  typedef pVectGlobalManip<UInt,NODEMAP>                        MANIP_ISVERTEX;
  typedef pGraphGlobalManip<GEOELEMENT,pMapItemShare,NODEMAP>   MANIP_ELEMENT;
  
  //Manipulators-----------------------------------------------------
  NODESVECT newNodes;
  BOOLVECT  newIsVertex;
  GRAPH1D   newElements;
  
  MANIP_POINT    manipPoint;
  MANIP_ISVERTEX manipIsVertex;
  MANIP_ELEMENT  manipElement;
  
  //Copy to new communicator-----------------------------------------
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() < OldCommDev.size());
    
    //Alloc
    NewMesh.clear();
    NewMesh = OldMesh;
    
    //Reduce comm
    manipPoint.reduceCommunicator(isActive,
                                  OldCommDev,
                                  OldMesh.getNodes(),
                                  NewCommDev,
                                  newNodes);
    
    manipIsVertex.reduceCommunicator(isActive,
                                     OldCommDev,
                                     OldMesh.getIsVertex(),
                                     NewCommDev,
                                     newIsVertex);
    
    manipElement.reduceCommunicator(isActive,
                                    OldCommDev,
                                    OldMesh.getElements(),
                                    NewCommDev,
                                    newElements);
    
    //Fix the mesh
    NewMesh.setNodes(newNodes,newIsVertex);
    NewMesh.setElements(newElements);
    
    NewMesh.computeNumVertices();
    NewMesh.transferMap();
  }
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MESH1D>       & OldMesh,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MESH1D>             & NewMesh)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldMesh,
                    *NewCommDev,
                    *NewMesh);
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh1dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
expandCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MESH1D       & OldMesh,
                   const communicator & NewCommDev,
                         MESH1D       & NewMesh)
{
  //Typedefs---------------------------------------------------------
  typedef typename MESH1D::NODESVECT NODESVECT;
  typedef typename MESH1D::BOOLVECT  BOOLVECT;
  typedef typename MESH1D::GRAPH1D   GRAPH1D;
  
  typedef pVectGlobalManip<point3d,NODEMAP>                     MANIP_POINT;
  typedef pVectGlobalManip<UInt,NODEMAP>                        MANIP_ISVERTEX;
  typedef pGraphGlobalManip<GEOELEMENT,pMapItemShare,NODEMAP>   MANIP_ELEMENT;
  
  //Manipulators-----------------------------------------------------
  NODESVECT newNodes;
  BOOLVECT  newIsVertex;
  GRAPH1D   newElements;
  
  MANIP_POINT    manipPoint;
  MANIP_ISVERTEX manipIsVertex;
  MANIP_ELEMENT  manipElement;
  
  //Copy to new communicator-----------------------------------------
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() > OldCommDev.size());
    
    //Alloc
    NewMesh.clear();
    NewMesh = OldMesh;
    
    //Reduce comm
    manipPoint.expandCommunicator(isActive,
                                  OldCommDev,
                                  OldMesh.getNodes(),
                                  NewCommDev,
                                  newNodes);
    
    manipIsVertex.expandCommunicator(isActive,
                                     OldCommDev,
                                     OldMesh.getIsVertex(),
                                     NewCommDev,
                                     newIsVertex);
    
    manipElement.expandCommunicator(isActive,
                                    OldCommDev,
                                    OldMesh.getElements(),
                                    NewCommDev,
                                    newElements);
    
    //Fix the mesh
    NewMesh.setNodes(newNodes,newIsVertex);
    NewMesh.setElements(newElements);
  }
  
  NewMesh.computeNumVertices();
  NewMesh.transferMap();
}

#endif
