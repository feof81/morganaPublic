/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MESH2DGLOBALMANIP_HPP
#define MESH2DGLOBALMANIP_HPP

#include <map>
#include "parmetis.h"

#include "pVectGlobalManip.hpp"
#include "pGraphGlobalManip.hpp"
#include "mesh2d.hpp"
#include "traitsGeometry.hpp"


//_________________________________________________________________________________________________
// UNSPECIALIZED
//-------------------------------------------------------------------------------------------------

/*! Global mesh manipulator 2d - unspecialized, empty */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP> class mesh2dGlobalManip
{ };


//_________________________________________________________________________________________________
// PMAPITEM
//-------------------------------------------------------------------------------------------------

/*! Global mesh manipulator 2d - \c pMapItem */
template<typename GEOSHAPE, typename NODEMAP>
class mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItem ELMAP;
    typedef mesh2d<GEOSHAPE,ELMAP,NODEMAP> MESH2D;
    
    typedef typename GEOSHAPE::GEOBSHAPE   GEOSHAPE1D;
    
    typedef geoElement<GEOSHAPE>   GEOELEMENT;
    typedef geoElement<GEOSHAPE1D> GEOELEMENT1D;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
  
    /*! @name Constructors */ //@{
  public:
    mesh2dGlobalManip();
    mesh2dGlobalManip(const Teuchos::RCP<const communicator> & CommDev);
    mesh2dGlobalManip(const communicator & CommDev);
    //@}
    
    /*! @name Information functions */ //@{
  public:
    /*! Check the graphs, the equality of the col maps of the elements, faces and edges with the map of the points */
    bool check(const Teuchos::RCP<const MESH2D> & Grid2d) const;
    bool check(const MESH2D & Grid2d) const;
    
    /*! Return the total number of vertices on the processors */
    UInt getNumGlobalVertices(const Teuchos::RCP<const MESH2D> & Grid2d) const;
    UInt getNumGlobalVertices(const MESH2D & Grid2d) const;
    
    /*! Return the total numer of edges on the processors */
    UInt getNumGlobalEdges(const Teuchos::RCP<const MESH2D> & Grid2d) const;
    UInt getNumGlobalEdges(const MESH2D & Grid2d) const;
    
    /*! Return the total numer of elements on the processors */
    UInt getNumGlobalElements(const Teuchos::RCP<const MESH2D> & Grid2d) const;
    UInt getNumGlobalElements(const MESH2D & Grid2d) const;
    //@}
    
    /*! @name Operation functions */ //@{
  public:
    /*! The mesh should be in global numbering and a new mesh in local numbering is built,
    no nodes or elements should be repeated across the processors */
    void meshPartition(const Teuchos::RCP<MESH2D> & Grid2d) const;
    void meshPartition(MESH2D & Grid2d) const;
    
    /*! The mesh should be in local numbering, no element repetition is allowed */
    void meshBalancing(const Teuchos::RCP<MESH2D> & Grid2d) const;
    void meshBalancing(MESH2D & grid2d) const;
    
    /*! The mesh should be in local numbering, no element repetition is allowed */
    void meshBalancing(const Teuchos::RCP<MESH2D> & Grid2d, const sVect<UInt> & elWeights) const;
    void meshBalancing(MESH2D & grid2d, const sVect<UInt> & elWeights) const;
    
    /*! The mesh should be in local numbering. The appendended elements and nodes are listed below the original ones */
    void meshOverlapping(const Teuchos::RCP<MESH2D> & Grid2d) const;
    void meshOverlapping(MESH2D & Grid2d) const;
    
    /*! Reduces the mesh to a single process and pushes the \c cids to the global numbering.
    Duplicated elements or points are eliminated. */
    void gather(const UInt & gatherGid, const Teuchos::RCP<MESH2D> & targetGrid2d, const Teuchos::RCP<const MESH2D> & sourceGrid2d);
    void gather(const UInt & gatherGid, MESH2D & targetGrid2d, const MESH2D & sourceGrid2d);
    
    /*! Eliminates the overlapped regions. ATTENTION: since there is no history of the way
    the mesh was overlapped. The shared elements are put in a particular domain in an
    un-predictable manner */
    void destroyOverlap(const Teuchos::RCP<MESH2D> & grid2d) const;
    void destroyOverlap(MESH2D & grid2d) const;
    //@}
    
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MESH2D>       & OldMesh,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MESH2D>             & NewMesh);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MESH2D       & OldMesh,
                            const communicator & NewCommDev,
                                  MESH2D       & NewMesh);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MESH2D>       & OldMesh,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MESH2D>             & NewMesh);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MESH2D       & OldMesh,
                            const communicator & NewCommDev,
                                  MESH2D       & NewMesh);
    //@}
};

template<typename GEOSHAPE, typename NODEMAP>
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
mesh2dGlobalManip()
{
  assert(GEOSHAPE::nDim == 2);
  commDevLoaded = false;
}

template<typename GEOSHAPE, typename NODEMAP>
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
mesh2dGlobalManip(const Teuchos::RCP<const communicator> & CommDev)
{
  assert(GEOSHAPE::nDim == 2);
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename NODEMAP>
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
mesh2dGlobalManip(const communicator & CommDev)
{
  assert(GEOSHAPE::nDim == 2);
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename NODEMAP>
bool
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
check(const Teuchos::RCP<const MESH2D> & Grid2d) const
{
  return(check(*Grid2d));
}

template<typename GEOSHAPE, typename NODEMAP>
bool
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
check(const MESH2D & Grid2d) const
{
  assert(Grid2d.getMapTransferred());
  
  //Allocations
  bool flag;
  UInt pid  = commDev->rank();
  UInt size = commDev->size();
  
  //Nodes vect check
  pVectGlobalManip<point3d,NODEMAP> nodesChecking(commDev);
  bool flagN = nodesChecking.check(Grid2d.getNodes());
  
  if(!flagN)
  { cout << "ERROR: Nodes check failed. Pid: " << pid << endl; }
  
  //Elements check
  pGraphGlobalManip<GEOELEMENT,pMapItem,NODEMAP> elementsChecking(commDev);
  bool flagE = elementsChecking.check(Grid2d.getElements());
  
  if(!flagE)
  { cout << "ERROR: Elements check failed. Pid: " << pid << endl; }
  
  //Edges check
  pGraphGlobalManip<GEOELEMENT1D,pMapItem,NODEMAP> edgesChecking(commDev);
  bool flagD = edgesChecking.check(Grid2d.getEdges());
  
  if(!flagD)
  { cout << "ERROR: Edges check failed. Pid: " << pid << endl; }
  
  flag = flagN & flagE & flagD;
  
  //Cross map checking
  pMapManip<NODEMAP> crossMapChecker;
  
  //Checking against elements col map
  sVect<UInt> allFlags(size);
  UInt myFlag = crossMapChecker.isEqual(Grid2d.getNodes().getMapRef(),Grid2d.getElements().getColMap());
  
  if(myFlag != 1)
  { cout << "ERROR: Elements col map differs from Nodes col map. Pid: " << pid << endl; }
  
  all_gather(*commDev, myFlag, allFlags);
  
  for(UInt i=1; i <= size; ++i)
  { myFlag = myFlag & allFlags(i); }
  
  flag = flag & myFlag;
  
  //Checking against edges col map
  if(getNumGlobalEdges(Grid2d) != 0)
  { myFlag = crossMapChecker.isEqual(Grid2d.getNodes().getMapRef(),Grid2d.getEdges().getColMap()); }
  else
  { myFlag = true; }
    
  all_gather(*commDev, myFlag, allFlags);
  
  for(UInt i=1; i <= size; ++i)
  { myFlag = myFlag & allFlags(i); }
  
  flag = flag & myFlag;
    
  if(myFlag != 1)
  { cout << "ERROR: Edges col map differs from Nodes col map. Pid: " << pid << endl; }
  
  //Cross check among processors
  sVect<UInt> flags;
  all_gather(*commDev,UInt(flag),flags);
  
  for(UInt i=1; i <= size; ++i)
  { flag = flag & bool(flags(i));}
  
  return(flag);
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalVertices(const Teuchos::RCP<const MESH2D> & Grid2d) const
{
  pVectGlobalManip<point3d,NODEMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid2d->getNodes()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalVertices(const MESH2D & Grid2d) const
{
  pVectGlobalManip<point3d,NODEMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid2d.getNodes()));
}    

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalEdges(const Teuchos::RCP<const MESH2D> & Grid2d) const
{
  pVectGlobalManip<GEOELEMENT1D,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid2d->getEdges()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalEdges(const MESH2D & Grid2d) const
{
  pVectGlobalManip<GEOELEMENT1D,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid2d.getEdges()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalElements(const Teuchos::RCP<const MESH2D> & Grid2d) const
{
  pVectGlobalManip<GEOELEMENT,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid2d->getElements()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalElements(const MESH2D & Grid2d) const
{
  pVectGlobalManip<GEOELEMENT,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid2d.getElements()));
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshPartition(const Teuchos::RCP<MESH2D> & Grid2d) const
{
  assert(commDevLoaded);
  meshPartition(*Grid2d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshPartition(mesh2d<GEOSHAPE,ELMAP,NODEMAP> & Grid2d) const
{
  //Asserts
  assert(commDevLoaded);
  assert(Grid2d.getElements().colIsLocal() == false);
  assert(Grid2d.getMeshStandard() == STDL);
  
  //Clearing
  Grid2d.clearEdges();
  Grid2d.transferMap();
  
  //PARMETIS_______________________________________________________________________________________
  UInt commSize = commDev->size();
  
  //Sizes communication
  UInt mySize = Grid2d.getNumElements();
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
      node    = Grid2d.getElementL(i).getCid(j) - 1; //Change numbering      
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
  idx_t ncommonnodes = GEOSHAPE1D::numPoints;
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
  
  GRAPH elGraph = Grid2d.getElements();
    
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
  Grid2d.setElements(elGraph);

  //Change the map to the nodes
  pVect<point3d,NODEMAP> nodesVect = Grid2d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Insert the nodes
  Grid2d.setNodes(nodesVect);
  
  //Internal methods
  Grid2d.computeNumVertices();
  Grid2d.transferMap();
  Grid2d.setMeshStandard(STDB);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshBalancing(const Teuchos::RCP<MESH2D> & Grid2d)  const
{
  assert(commDevLoaded);
  meshBalancing(*Grid2d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshBalancing(mesh2d<GEOSHAPE,ELMAP,NODEMAP> & Grid2d) const
{
  assert(commDevLoaded);
  assert(Grid2d.getElements().colIsLocal());
  assert(Grid2d.getMeshStandard() == STDB);
  
  //Clearing
  Grid2d.clearEdges();
  Grid2d.transferMap();
  
  //PARMETIS_______________________________________________________________________________________
  UInt commSize = commDev->size();
  
  //Sizes communication
  UInt mySize = Grid2d.getNumElements();
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
      node    = Grid2d.getElements().getCid_LG(i,j) - 1; //Change numbering      
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
  idx_t ncommonnodes = GEOSHAPE1D::numPoints;
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
  
  GRAPH elGraph = Grid2d.getElements();
    
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
  Grid2d.clearEdges();  
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
  Grid2d.setElements(elGraph);

  //Change the map to the nodes
  pVect<point3d,NODEMAP> nodesVect = Grid2d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Insert the nodes
  Grid2d.setNodes(nodesVect);
  
  //Internal methods
  Grid2d.computeNumVertices();
  Grid2d.transferMap();
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshBalancing(const Teuchos::RCP<MESH2D> & Grid2d, const sVect<UInt> & elWeights) const
{
  assert(commDevLoaded);
  meshBalancing(*Grid2d,elWeights);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshBalancing(MESH2D & Grid2d, const sVect<UInt> & elWeights) const
{
  assert(commDevLoaded);
  assert(Grid2d.getElements().colIsLocal());
  assert(Grid2d.getMeshStandard() == STDB);
  
  //Clearing
  Grid2d.clearEdges();
  Grid2d.transferMap();
  
  //PARMETIS_______________________________________________________________________________________
  UInt commSize = commDev->size();
  
  //Sizes communication
  UInt mySize = Grid2d.getNumElements();
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
      node    = Grid2d.getElements().getCid_LG(i,j) - 1; //Change numbering      
      eind[k] = node;
      k++;
    }
    
    assert(k == UInt(eptr[i]));
  }
  
  //Weights flags
  idx_t * elmwgt = new idx_t[mySize];
  
  for(UInt i=0; i < mySize; ++i)
  { elmwgt[i] = elWeights(i+1); }
  
  //Some standard flags
  idx_t wgtflag = 2;
  idx_t numflag = 0;
  idx_t ncon    = 1;
  
  //Other flags
  idx_t ncommonnodes = GEOSHAPE1D::numPoints;
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
  
  GRAPH elGraph = Grid2d.getElements();
    
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
  Grid2d.clearEdges();  
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
  Grid2d.setElements(elGraph);

  //Change the map to the nodes
  pVect<point3d,NODEMAP> nodesVect = Grid2d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Insert the nodes
  Grid2d.setNodes(nodesVect);
  
  //Internal methods
  Grid2d.computeNumVertices();
  Grid2d.transferMap();
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshOverlapping(const Teuchos::RCP<MESH2D> & Grid2d) const
{
  assert(commDevLoaded);
  meshOverlapping(*Grid2d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshOverlapping(mesh2d<GEOSHAPE,ELMAP,NODEMAP> & Grid2d) const
{
  assert(commDevLoaded);
  assert(Grid2d.getElements().colIsLocal());
  assert(Grid2d.getMeshStandard() == STDB);
  
  //Clearing
  Grid2d.clearEdges();
  Grid2d.transferMap();
  
  //Elements overlapping
  typedef pGraph<GEOELEMENT,ELMAP,NODEMAP>  GRAPH;
  
  GRAPH elGraph = Grid2d.getElements();
  
  pGraphGlobalManip<GEOELEMENT,pMapItem,NODEMAP> elGraphManip(commDev);
  elGraphManip.overlap(elGraph);
  
  //Nodes overlapping
  pMap<NODEMAP>          colMap    = elGraph.getColMap();
  pVect<point3d,NODEMAP> nodesVect = Grid2d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Grid fill
  Grid2d.setElements(elGraph);
  Grid2d.setNodes(nodesVect);
  Grid2d.transferMap();
  Grid2d.setMeshStandard(STDA);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
gather(const UInt & gatherGid, const Teuchos::RCP<MESH2D> & targetGrid2d, const Teuchos::RCP<const MESH2D> & sourceGrid2d)
{
  assert(commDevLoaded);
  gather(gatherGid,*targetGrid2d,*sourceGrid2d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
gather(const UInt & gatherGid, MESH2D & targetGrid2d, const MESH2D & sourceGrid2d)
{
  typedef typename MESH2D::NODESVECT    NODESVECT;
  typedef typename MESH2D::GRAPH2D      GRAPH2D;
  typedef typename MESH2D::GEOELEMENT2D GEOELEMENT2D;
  
  //Mesh cheking
  assert(commDevLoaded);
  assert(sourceGrid2d.getElements().colIsLocal());
  
  //Nodes and graph extraction
  NODESVECT nodes    = sourceGrid2d.getNodes();
  GRAPH2D   elements = sourceGrid2d.getElements();
  
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
  pVectComm<point3d,NODEMAP> nodesComm(commDev);
  nodesComm.vectorPid(nodes);

  pVectComm<GEOELEMENT2D,ELMAP> elementsComm(commDev);
  elementsComm.vectorPid(elements);
  
  
  //Nodes - Elements exclusive merging
  pVectManip<point3d,NODEMAP> nodesManip;
  nodesManip.unionExclusive(nodes);  
  
  pVectManip<GEOELEMENT2D,ELMAP> elementsManip;
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
  
  targetGrid2d.clear();
  targetGrid2d.setNodes(nodes);
  targetGrid2d.setElements(elements);
  targetGrid2d.transferMap();
  targetGrid2d.setMeshStandard(STDU);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
destroyOverlap(const Teuchos::RCP<MESH2D> & grid2d) const
{
  assert(commDevLoaded);
  destroyOverlap(*grid2d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
destroyOverlap(MESH2D & grid2d) const
{
  assert(commDevLoaded);
  assert(grid2d.getMeshStandard() == STDA);
  
  //Mesh cheking-----------------------------------------------------------------------------------
  assert(commDevLoaded);
  assert(grid2d.getElements().colIsLocal());
  
  //Typedefs
  typedef typename MESH2D::GRAPH2D    GRAPH2D;
  typedef typename MESH2D::NODESVECT  NODESVECT;
  
  //Extract data
  GRAPH2D elements = grid2d.getElements();
  NODESVECT  nodes = grid2d.getNodes();
  
  
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
  GRAPH2D  newElements;
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
    for(UInt j=1; j <= grid2d.getNumPoints(); ++j)  //I points sono locali
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
  geoElement<GEOSHAPE> geoElement3d(true);
  
  for(UInt i=1; i <= newElements.size(); ++i)
  {
    geoElement3d = newElements.getItemL(i);
    
    for(UInt j=1; j <= grid2d.getNumPoints(); ++j)
    {
      id = newElements.getCid_LL(i,j);
      geoElement3d.setCid(j,newOrder(id));
    }
    
    newElements(i) = geoElement3d;
  } 
  
  
  //Rebuild the owning-sharing for the nodes-------------------------------------------------------
  nodesMapFixer(newNodes.getMapRef(),*commDev);
  
  //Build new grid---------------------------------------------------------------------------------
  newNodes.updateFinder();
  
  newElements.setColMap(newNodes.getMapRef());
  newElements.updateRowFinder();
  newElements.updateColFinder();
  
  grid2d.clearEdges();
  
  grid2d.setNodes(newNodes);
  grid2d.setElements(newElements);
  grid2d.computeNumVertices();
  grid2d.transferMap();
  grid2d.setMeshStandard(STDB);
  
  //Final check------------------------------------------------------------------------------------
  assert(check(grid2d));
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MESH2D>       & OldMesh,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MESH2D>             & NewMesh)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldMesh,
                    *NewCommDev,
                    *NewMesh);
}
  
template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
reduceCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MESH2D       & OldMesh,
                   const communicator & NewCommDev,
                         MESH2D       & NewMesh)
{
  //Typedefs---------------------------------------------------------
  typedef typename MESH2D::NODESVECT NODESVECT;
  typedef typename MESH2D::BOOLVECT  BOOLVECT;
  typedef typename MESH2D::GRAPH2D   GRAPH2D;
  typedef typename MESH2D::GRAPH1D   GRAPH1D;
  
  typedef pVectGlobalManip<point3d,NODEMAP>                MANIP_POINT;
  typedef pVectGlobalManip<UInt,NODEMAP>                   MANIP_ISVERTEX;
  typedef pGraphGlobalManip<GEOELEMENT,pMapItem,NODEMAP>   MANIP_ELEMENT;
  typedef pGraphGlobalManip<GEOELEMENT1D,pMapItem,NODEMAP> MANIP_EDGE;
  
  //Manipulators-----------------------------------------------------
  NODESVECT newNodes;
  BOOLVECT  newIsVertex;
  GRAPH2D   newElements;
  GRAPH1D   newEdges;
  
  MANIP_POINT    manipPoint;
  MANIP_ISVERTEX manipIsVertex;
  MANIP_ELEMENT  manipElement;
  MANIP_EDGE     manipEdge;
  
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
    
    manipEdge.reduceCommunicator(isActive,
                                 OldCommDev,
                                 OldMesh.getEdges(),
                                 NewCommDev,
                                 newEdges);
    
    //Fix the mesh
    NewMesh.setNodes(newNodes,newIsVertex);
    NewMesh.setElements(newElements);
    NewMesh.setEdges(newEdges);
    
    NewMesh.computeNumVertices();
    NewMesh.transferMap();
  }
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MESH2D>       & OldMesh,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MESH2D>             & NewMesh)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldMesh,
                    *NewCommDev,
                    *NewMesh);
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
expandCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MESH2D       & OldMesh,
                   const communicator & NewCommDev,
                         MESH2D       & NewMesh)
{
  //Typedefs---------------------------------------------------------
  typedef typename MESH2D::NODESVECT NODESVECT;
  typedef typename MESH2D::BOOLVECT  BOOLVECT;
  typedef typename MESH2D::GRAPH2D   GRAPH2D;
  typedef typename MESH2D::GRAPH1D   GRAPH1D;
  
  typedef pVectGlobalManip<point3d,NODEMAP>                MANIP_POINT;
  typedef pVectGlobalManip<UInt,NODEMAP>                   MANIP_ISVERTEX;
  typedef pGraphGlobalManip<GEOELEMENT,pMapItem,NODEMAP>   MANIP_ELEMENT;
  typedef pGraphGlobalManip<GEOELEMENT1D,pMapItem,NODEMAP> MANIP_EDGE;
  
  //Manipulators-----------------------------------------------------
  NODESVECT newNodes;
  BOOLVECT  newIsVertex;
  GRAPH2D   newElements;
  GRAPH1D   newEdges;
  
  MANIP_POINT    manipPoint;
  MANIP_ISVERTEX manipIsVertex;
  MANIP_ELEMENT  manipElement;
  MANIP_EDGE     manipEdge;
  
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
    
    manipEdge.expandCommunicator(isActive,
                                 OldCommDev,
                                 OldMesh.getEdges(),
                                 NewCommDev,
                                 newEdges);
    
    //Fix the mesh
    NewMesh.setNodes(newNodes,newIsVertex);
    NewMesh.setElements(newElements);
    NewMesh.setEdges(newEdges);
  }
  
  NewMesh.computeNumVertices();
  NewMesh.transferMap();
}



//_________________________________________________________________________________________________
// PMAPITEMSHARE
//-------------------------------------------------------------------------------------------------

/*! Global mesh manipulator 2d - \c pMapItemShare */
template<typename GEOSHAPE, typename NODEMAP>
class mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItemShare ELMAP;
    typedef mesh2d<GEOSHAPE,ELMAP,NODEMAP> MESH2D;
    
    typedef typename GEOSHAPE::GEOBSHAPE   GEOSHAPE1D;
    
    typedef geoElement<GEOSHAPE>   GEOELEMENT;
    typedef geoElement<GEOSHAPE1D> GEOELEMENT1D;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
  
    /*! @name Constructors */ //@{
  public:
    mesh2dGlobalManip();
    mesh2dGlobalManip(const Teuchos::RCP<const communicator> & CommDev);
    mesh2dGlobalManip(const communicator & CommDev);
    //@}
    
    /*! @name Information functions */ //@{
  public:
    /*! Check the graphs, the equality of the col maps of the elements, faces and edges with the map of the points */
    bool check(const Teuchos::RCP<const MESH2D> & Grid2d) const;
    bool check(const MESH2D & Grid2d) const;
    
    /*! Return the total number of vertices on the processors */
    UInt getNumGlobalVertices(const Teuchos::RCP<const MESH2D> & Grid2d) const;
    UInt getNumGlobalVertices(const MESH2D & Grid2d) const;
    
    /*! Return the total numer of edges on the processors */
    UInt getNumGlobalEdges(const Teuchos::RCP<const MESH2D> & Grid2d) const;
    UInt getNumGlobalEdges(const MESH2D & Grid2d) const;
    
    /*! Return the total numer of elements on the processors */
    UInt getNumGlobalElements(const Teuchos::RCP<const MESH2D> & Grid2d) const;
    UInt getNumGlobalElements(const MESH2D & Grid2d) const;
    //@}
    
    /*! @name Operation functions */ //@{
  public:
    /*! The mesh should be in global numbering and a new mesh in local numbering is built */
    void meshPartition(const Teuchos::RCP<MESH2D> & Grid2d) const;
    void meshPartition(MESH2D & Grid2d) const;
    
    /*! The mesh should be in local numbering, the not-owned elements are discarded */
    void meshBalancing(const Teuchos::RCP<MESH2D> & Grid2d) const;
    void meshBalancing(MESH2D & Grid2d) const;
    
    /*! The mesh should be in local numbering, no element repetition is allowed */
    void meshBalancing(const Teuchos::RCP<MESH2D> & Grid2d, const sVect<UInt> & elWeights) const;
    void meshBalancing(MESH2D & Grid2d, const sVect<UInt> & elWeights) const;
    
    /*! The mesh should be in local numbering. The appendended elements and nodes are listed below the original ones */
    void meshOverlapping(const Teuchos::RCP<MESH2D> & Grid2d) const;
    void meshOverlapping(MESH2D & Grid2d) const;
    
    /*! Reduces the mesh to a single process and pushes the \c cids to the global numbering.
    Duplicated elements or points are eliminated. */
    void gather(const UInt & gatherGid, const Teuchos::RCP<MESH2D> & targetGrid2d, const Teuchos::RCP<const MESH2D> & sourceGrid2d);
    void gather(const UInt & gatherGid, MESH2D & targetGrid2d, const MESH2D & sourceGrid2d);
    
    /*! Eliminates the overlapped regions. The shared elements are eliminated */
    void destroyOverlap(const Teuchos::RCP<MESH2D> & grid2d) const;
    void destroyOverlap(MESH2D & grid2d) const;
    //@}
    
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MESH2D>       & OldMesh,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MESH2D>             & NewMesh);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MESH2D       & OldMesh,
                            const communicator & NewCommDev,
                                  MESH2D       & NewMesh);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MESH2D>       & OldMesh,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MESH2D>             & NewMesh);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MESH2D       & OldMesh,
                            const communicator & NewCommDev,
                                  MESH2D       & NewMesh);
    //@}
};


template<typename GEOSHAPE, typename NODEMAP>
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
mesh2dGlobalManip()
{
  assert(GEOSHAPE::nDim == 2);
  commDevLoaded = false;
}

template<typename GEOSHAPE, typename NODEMAP>
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
mesh2dGlobalManip(const Teuchos::RCP<const communicator> & CommDev)
{
  assert(GEOSHAPE::nDim == 2);
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename NODEMAP>
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
mesh2dGlobalManip(const communicator & CommDev)
{
  assert(GEOSHAPE::nDim == 2);
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename NODEMAP>
bool
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
check(const Teuchos::RCP<const MESH2D> & Grid2d) const
{
  return(check(*Grid2d));
}

template<typename GEOSHAPE, typename NODEMAP>
bool
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
check(const MESH2D & Grid2d) const
{
  assert(Grid2d.getMapTransferred());
  
  //Allocations
  bool flag;
  UInt pid  = commDev->rank();
  UInt size = commDev->size();
  
  //Nodes vect check
  pVectGlobalManip<point3d,NODEMAP> nodesChecking(commDev);
  bool flagN = nodesChecking.check(Grid2d.getNodes());
  
  if(!flagN)
  { cout << "ERROR: Nodes check failed. Pid: " << pid << endl; }
  
  //Elements check
  pGraphGlobalManip<GEOELEMENT,pMapItemShare,NODEMAP> elementsChecking(commDev);
  bool flagE = elementsChecking.check(Grid2d.getElements());
  
  if(!flagE)
  { cout << "ERROR: Elements check failed. Pid: " << pid << endl; }
  
  //Edges check
  pGraphGlobalManip<GEOELEMENT1D,pMapItemShare,NODEMAP> edgesChecking(commDev);
  bool flagD = edgesChecking.check(Grid2d.getEdges());
  
  if(!flagD)
  { cout << "ERROR: Edges check failed. Pid: " << pid << endl; }
  
  flag = flagN & flagE & flagD;
  
  //Cross map checking
  pMapManip<NODEMAP> crossMapChecker;
  
  //Checking against elements col map
  sVect<UInt> allFlags(size);
  UInt myFlag = crossMapChecker.isEqual(Grid2d.getNodes().getMapRef(),Grid2d.getElements().getColMap());
  
  if(myFlag != 1)
  { cout << "ERROR: Elements col map differs from Nodes col map. Pid: " << pid << endl; }
  
  all_gather(*commDev, myFlag, allFlags);
  
  for(UInt i=1; i <= size; ++i)
  { myFlag = myFlag & allFlags(i); }
  
  flag = flag & myFlag;
  
  //Checking against edges col map
  if(getNumGlobalEdges(Grid2d) != 0)
  { myFlag = crossMapChecker.isEqual(Grid2d.getNodes().getMapRef(),Grid2d.getEdges().getColMap()); }
  else
  { myFlag = true; }
    
  all_gather(*commDev, myFlag, allFlags);
  
  for(UInt i=1; i <= size; ++i)
  { myFlag = myFlag & allFlags(i); }
  
  flag = flag & myFlag;
    
  if(myFlag != 1)
  { cout << "ERROR: Edges col map differs from Nodes col map. Pid: " << pid << endl; }
  
  
  //Cross checking for ownership - elements
  auxiliaryGraphCheck2d<GEOSHAPE,ELMAP,NODEMAP> elOwCheker;
  bool flagWE = elOwCheker.check(Grid2d.getElements());
  
  flag = flag && flagWE;
  
  if(!flagWE)
  { cout << "ERROR: Elements cross ownership failed. Pid: " << pid << endl; }  
  
  //Cross checking for ownership - faces
  if(Grid2d.getEdges().size() != 0)
  {
    auxiliaryGraphCheck2d<GEOSHAPE1D,ELMAP,NODEMAP> edOwCheker;
    bool flagWD = edOwCheker.check(Grid2d.getEdges());
     
    flag = flag && flagWD;
  
    if(!flagWD)
    { cout << "Edges: Faces cross ownership failed. Pid: " << pid << endl; }
  }
  
  
  //Cross check among processors
  sVect<UInt> flags;
  all_gather(*commDev,UInt(flag),flags);
  
  for(UInt i=1; i <= size; ++i)
  { flag = flag & bool(flags(i));}
  
  return(flag);
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalVertices(const Teuchos::RCP<const MESH2D> & Grid2d) const
{
  pVectGlobalManip<point3d,NODEMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid2d->getNodes()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalVertices(const MESH2D & Grid2d) const
{
  pVectGlobalManip<point3d,NODEMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid2d.getNodes()));
}    

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalEdges(const Teuchos::RCP<const MESH2D> & Grid2d) const
{
  pVectGlobalManip<GEOELEMENT1D,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid2d->getEdges()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalEdges(const MESH2D & Grid2d) const
{
  pVectGlobalManip<GEOELEMENT1D,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid2d.getEdges()));
}
    
template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalElements(const Teuchos::RCP<const MESH2D> & Grid2d) const
{
  pVectGlobalManip<GEOELEMENT,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid2d->getElements()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalElements(const MESH2D & Grid2d) const
{
  pVectGlobalManip<GEOELEMENT,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid2d.getElements()));
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshPartition(const Teuchos::RCP<MESH2D> & Grid2d) const
{
  assert(commDevLoaded);
  meshPartition(*Grid2d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshPartition(MESH2D & Grid2d) const
{
  //Asserts
  assert(commDevLoaded);
  assert(Grid2d.getElements().colIsLocal() == false);
  assert(Grid2d.getMeshStandard() == STDL);
  
  //Clearing
  Grid2d.clearEdges();
  Grid2d.transferMap();
  
  //Sharing-owningMap______________________________________________________________________________
  typedef typename MESH2D::GRAPH2D GRAPH2D;
  
  GRAPH2D tempElements = Grid2d.getElements();
  
  for(UInt i=1; i <= tempElements.getRowMap().size(); ++i)
  {
    tempElements.getRowMap().get(i).setShared(false);
    tempElements.getRowMap().get(i).setOwned(true);
  }
  
  Grid2d.setElements(tempElements);
  Grid2d.transferMap();
  
  //PARMETIS_______________________________________________________________________________________
  UInt commSize = commDev->size();
  
  //Sizes communication
  UInt mySize = Grid2d.getNumElements();
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
      node    = Grid2d.getElementL(i).getCid(j) - 1; //Change numbering      
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
  idx_t ncommonnodes = GEOSHAPE1D::numPoints;
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
  
  GRAPH elGraph = Grid2d.getElements();
    
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
  Grid2d.setElements(elGraph);

  //Change the map to the nodes
  pVect<point3d,NODEMAP> nodesVect = Grid2d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Insert the nodes
  Grid2d.setNodes(nodesVect);
  
  //Internal methods
  Grid2d.computeNumVertices();
  Grid2d.transferMap();
  Grid2d.setMeshStandard(STDB);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshBalancing(const Teuchos::RCP<MESH2D> & Grid2d) const
{
  assert(commDevLoaded);
  meshBalancing(*Grid2d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshBalancing(MESH2D & Grid2d) const
{
  assert(commDevLoaded);
  assert(Grid2d.getMeshStandard() == STDB);
  
  typedef typename MESH2D::GRAPH2D      ELGRAPH;
  typedef typename MESH2D::GEOELEMENT2D ELEMENT;
  typedef typename MESH2D::NODESVECT    NODESVECT;
  
  //Clearing
  Grid2d.clearEdges();
  Grid2d.transferMap();
  
  //Grid copy
  ELGRAPH   elementsCopy = Grid2d.getElements();
  ELGRAPH   elList;
  NODESVECT nodesCopy = Grid2d.getNodes();;
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
  MESH2D tempGrid;
  
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
  Grid2d = tempGrid;
  Grid2d.transferMap();
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshBalancing(const Teuchos::RCP<MESH2D> & Grid2d, const sVect<UInt> & elWeights) const
{
  assert(commDevLoaded);
  meshBalancing(*Grid2d,elWeights);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshBalancing(MESH2D & Grid2d, const sVect<UInt> & elWeights) const
{
  assert(commDevLoaded);
  assert(Grid2d.getElements().colIsLocal());
  assert(Grid2d.getMeshStandard() == STDB);
  
  //Clearing
  Grid2d.clearEdges();
  Grid2d.transferMap();
  
  //PARMETIS_______________________________________________________________________________________
  UInt commSize = commDev->size();
  
  //Sizes communication
  UInt mySize = Grid2d.getNumElements();
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
      node    = Grid2d.getElements().getCid_LG(i,j) - 1; //Change numbering      
      eind[k] = node;
      k++;
    }
    
    assert(k == UInt(eptr[i]));
  }
  
  //Weights flags
  idx_t * elmwgt = new idx_t[mySize];
  
  for(UInt i=0; i < mySize; ++i)
  { elmwgt[i] = elWeights(i+1); }
  
  //Some standard flags
  idx_t wgtflag = 2;
  idx_t numflag = 0;
  idx_t ncon    = 1;
  
  //Other flags
  idx_t ncommonnodes = GEOSHAPE1D::numPoints;
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
  
  GRAPH elGraph = Grid2d.getElements();
    
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
  Grid2d.clearEdges();  
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
  Grid2d.setElements(elGraph);

  //Change the map to the nodes
  pVect<point3d,NODEMAP> nodesVect = Grid2d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Insert the nodes
  Grid2d.setNodes(nodesVect);
  
  //Internal methods
  Grid2d.computeNumVertices();
  Grid2d.transferMap();
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshOverlapping(const Teuchos::RCP<MESH2D> & Grid2d) const
{
  assert(commDevLoaded);
  meshOverlapping(*Grid2d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshOverlapping(MESH2D & Grid2d) const
{
  assert(commDevLoaded);
  assert(Grid2d.getElements().colIsLocal());
  assert(Grid2d.getMeshStandard() == STDB);
  
  //Clearing
  Grid2d.clearEdges();
  Grid2d.transferMap();
  
  //Elements overlapping
  typedef pGraph<GEOELEMENT,ELMAP,NODEMAP>  GRAPH;
  
  GRAPH elGraph = Grid2d.getElements();
  
  pGraphGlobalManip<GEOELEMENT,pMapItemShare,NODEMAP> elGraphManip(commDev);
  elGraphManip.overlap(elGraph);
  
  //Nodes overlapping
  pMap<NODEMAP>          colMap    = elGraph.getColMap();
  pVect<point3d,NODEMAP> nodesVect = Grid2d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Grid fill
  Grid2d.setElements(elGraph);
  Grid2d.setNodes(nodesVect);
  Grid2d.transferMap();
  Grid2d.setMeshStandard(STDA);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
gather(const UInt & gatherGid, const Teuchos::RCP<MESH2D> & targetGrid2d, const Teuchos::RCP<const MESH2D> & sourceGrid2d)
{
  assert(commDevLoaded);
  gather(gatherGid,*targetGrid2d,*sourceGrid2d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
gather(const UInt & gatherGid, MESH2D & targetGrid2d, const MESH2D & sourceGrid2d)
{
  typedef typename MESH2D::NODESVECT    NODESVECT;
  typedef typename MESH2D::GRAPH2D      GRAPH2D;
  typedef typename MESH2D::GEOELEMENT2D GEOELEMENT2D;
  
  //Mesh cheking
  assert(commDevLoaded);
  assert(sourceGrid2d.getElements().colIsLocal());
  
  //Nodes and graph extraction
  NODESVECT nodes    = sourceGrid2d.getNodes();
  GRAPH2D   elements = sourceGrid2d.getElements();
  
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

  pVectComm<GEOELEMENT2D,ELMAP> elementsComm(commDev);
  elementsComm.vectorPid(elements);
  
  
  //Nodes - Elements exclusive merging
  pVectManip<point3d,NODEMAP> nodesManip;
  nodesManip.unionExclusive(nodes);  
  
  pVectManip<GEOELEMENT2D,ELMAP> elementsManip;
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
  
  targetGrid2d.clear();
  targetGrid2d.setNodes(nodes);
  targetGrid2d.setElements(elements);
  targetGrid2d.transferMap();
  targetGrid2d.setMeshStandard(STDU);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
destroyOverlap(const Teuchos::RCP<MESH2D> & grid2d) const
{
  assert(commDevLoaded);
  destroyOverlap(*grid2d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
destroyOverlap(MESH2D & grid2d) const
{
  assert(grid2d.getMeshStandard() == STDA);
  
  //Mesh cheking-----------------------------------------------------------------------------------
  assert(commDevLoaded);
  assert(grid2d.getElements().colIsLocal());
  
  //Typedefs
  typedef typename MESH2D::GRAPH2D    GRAPH2D;
  typedef typename MESH2D::NODESVECT  NODESVECT;
  
  //Extract data
  GRAPH2D elements = grid2d.getElements();
  NODESVECT  nodes = grid2d.getNodes();
  
  
  //Eliminate the unused elements------------------------------------------------------------------
  UInt tot = 1;
  GRAPH2D       newElements;
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
    for(UInt j=1; j <= grid2d.getNumPoints(); ++j)  //I points sono locali
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
  geoElement<GEOSHAPE> geoElement3d(true);
  
  for(UInt i=1; i <= newElements.size(); ++i)
  {
    geoElement3d = newElements.getItemL(i);
    
    for(UInt j=1; j <= grid2d.getNumPoints(); ++j)
    {
      id = newElements.getCid_LL(i,j);
      geoElement3d.setCid(j,newOrder(id));
    }
    
    newElements(i) = geoElement3d;
  } 
  
  
  //Rebuild the owning-sharing for the nodes-------------------------------------------------------
  nodesMapFixer(newNodes.getMapRef(),*commDev);
  
  //Build new grid---------------------------------------------------------------------------------
  newNodes.updateFinder();
  
  newElements.setColMap(newNodes.getMapRef());
  newElements.updateRowFinder();
  newElements.updateColFinder();
  
  grid2d.clearEdges();
  
  grid2d.setNodes(newNodes);
  grid2d.setElements(newElements);
  grid2d.computeNumVertices();
  grid2d.transferMap();
  grid2d.setMeshStandard(STDB);
  
  //Final check------------------------------------------------------------------------------------
  assert(check(grid2d));
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MESH2D>       & OldMesh,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MESH2D>             & NewMesh)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldMesh,
                    *NewCommDev,
                    *NewMesh);
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
reduceCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MESH2D       & OldMesh,
                   const communicator & NewCommDev,
                         MESH2D       & NewMesh)
{
  //Typedefs---------------------------------------------------------
  typedef typename MESH2D::NODESVECT NODESVECT;
  typedef typename MESH2D::BOOLVECT  BOOLVECT;
  typedef typename MESH2D::GRAPH2D   GRAPH2D;
  typedef typename MESH2D::GRAPH1D   GRAPH1D;
  
  typedef pVectGlobalManip<point3d,NODEMAP>                     MANIP_POINT;
  typedef pVectGlobalManip<UInt,NODEMAP>                        MANIP_ISVERTEX;
  typedef pGraphGlobalManip<GEOELEMENT,pMapItemShare,NODEMAP>   MANIP_ELEMENT;
  typedef pGraphGlobalManip<GEOELEMENT1D,pMapItemShare,NODEMAP> MANIP_EDGE;
  
  //Manipulators-----------------------------------------------------
  NODESVECT newNodes;
  BOOLVECT  newIsVertex;
  GRAPH2D   newElements;
  GRAPH1D   newEdges;
  
  MANIP_POINT    manipPoint;
  MANIP_ISVERTEX manipIsVertex;
  MANIP_ELEMENT  manipElement;
  MANIP_EDGE     manipEdge;
  
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
    
    manipEdge.reduceCommunicator(isActive,
                                 OldCommDev,
                                 OldMesh.getEdges(),
                                 NewCommDev,
                                 newEdges);
    
    //Fix the mesh
    NewMesh.setNodes(newNodes,newIsVertex);
    NewMesh.setElements(newElements);
    NewMesh.setEdges(newEdges);
    
    NewMesh.computeNumVertices();
    NewMesh.transferMap();
  }
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MESH2D>       & OldMesh,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MESH2D>             & NewMesh)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldMesh,
                    *NewCommDev,
                    *NewMesh);
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh2dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
expandCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MESH2D       & OldMesh,
                   const communicator & NewCommDev,
                         MESH2D       & NewMesh)
{
  //Typedefs---------------------------------------------------------
  typedef typename MESH2D::NODESVECT NODESVECT;
  typedef typename MESH2D::BOOLVECT  BOOLVECT;
  typedef typename MESH2D::GRAPH2D   GRAPH2D;
  typedef typename MESH2D::GRAPH1D   GRAPH1D;
  
  typedef pVectGlobalManip<point3d,NODEMAP>                     MANIP_POINT;
  typedef pVectGlobalManip<UInt,NODEMAP>                        MANIP_ISVERTEX;
  typedef pGraphGlobalManip<GEOELEMENT,pMapItemShare,NODEMAP>   MANIP_ELEMENT;
  typedef pGraphGlobalManip<GEOELEMENT1D,pMapItemShare,NODEMAP> MANIP_EDGE;
  
  //Manipulators-----------------------------------------------------
  NODESVECT newNodes;
  BOOLVECT  newIsVertex;
  GRAPH2D   newElements;
  GRAPH1D   newEdges;
  
  MANIP_POINT    manipPoint;
  MANIP_ISVERTEX manipIsVertex;
  MANIP_ELEMENT  manipElement;
  MANIP_EDGE     manipEdge;
  
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
    
    manipEdge.expandCommunicator(isActive,
                                 OldCommDev,
                                 OldMesh.getEdges(),
                                 NewCommDev,
                                 newEdges);
    
    //Fix the mesh
    NewMesh.setNodes(newNodes,newIsVertex);
    NewMesh.setElements(newElements);
    NewMesh.setEdges(newEdges);
  }
  
  NewMesh.computeNumVertices();
  NewMesh.transferMap();
}

#endif
