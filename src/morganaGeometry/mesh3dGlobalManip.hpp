/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MESH3DGLOBALMANIP_HPP
#define MESH3DGLOBALMANIP_HPP

#include <map>
#include "parmetis.h"

#include "traitsGeometry.hpp"
#include "pMapGlobalManip.h"
#include "pVectGlobalManip.hpp"
#include "pGraphGlobalManip.hpp"
#include "mesh3d.hpp"


//_________________________________________________________________________________________________
// UNSPECIALIZED
//-------------------------------------------------------------------------------------------------

/*! Global mesh manipulator 3d - unspecialized, empty */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP> class mesh3dGlobalManip
{ };


//_________________________________________________________________________________________________
// PMAPITEM
//-------------------------------------------------------------------------------------------------

/*! Global mesh manipulator 3d - \c pMapItem */
template<typename GEOSHAPE, typename NODEMAP>
class mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItem ELMAP;
    typedef mesh3d<GEOSHAPE,ELMAP,NODEMAP> MESH3D;
    
    typedef typename GEOSHAPE::GEOBSHAPE   GEOSHAPE2D;
    typedef typename GEOSHAPE2D::GEOBSHAPE GEOSHAPE1D;
    
    typedef geoElement<GEOSHAPE>   GEOELEMENT;
    typedef geoElement<GEOSHAPE2D> GEOELEMENT2D;
    typedef geoElement<GEOSHAPE1D> GEOELEMENT1D;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
  
    /*! @name Constructors */ //@{
  public:
    mesh3dGlobalManip();
    mesh3dGlobalManip(const Teuchos::RCP<const communicator> & CommDev);
    mesh3dGlobalManip(const communicator & CommDev);
    //@}
    
    /*! @name Information functions */ //@{
  public:
    /*! Check the graphs, the equality of the col maps of the elements, faces and edges with the map of the points */
    bool check(const Teuchos::RCP<const MESH3D> & Grid3d) const;
    bool check(const MESH3D & Grid3d) const;
    
    /*! Return the total number of vertices on the processors */
    UInt getNumGlobalVertices(const Teuchos::RCP<const MESH3D> & Grid3d) const;
    UInt getNumGlobalVertices(const MESH3D & Grid3d) const;
    
    /*! Return the total numer of edges on the processors */
    UInt getNumGlobalEdges(const Teuchos::RCP<const MESH3D> & Grid3d) const;
    UInt getNumGlobalEdges(const MESH3D & Grid3d) const;
    
    /*! Return the total numer of faces on the processors */
    UInt getNumGlobalFaces(const Teuchos::RCP<const MESH3D> & Grid3d) const;
    UInt getNumGlobalFaces(const MESH3D & Grid3d) const;
    
    /*! Return the total numer of elements on the processors */
    UInt getNumGlobalElements(const Teuchos::RCP<const MESH3D> & Grid3d) const;
    UInt getNumGlobalElements(const MESH3D & Grid3d) const;
    //@}
    
    /*! @name Operation functions */ //@{
  public:
    /*! The mesh should be in global numbering and a new mesh in local numbering is built,
    no nodes or elements should be repeated across the processors */
    void meshPartition(const Teuchos::RCP<MESH3D> & Grid3d) const;
    void meshPartition(MESH3D & Grid3d) const;
    
    /*! The mesh should be in local numbering, no element repetition is allowed */
    void meshBalancing(const Teuchos::RCP<MESH3D> & Grid3d) const;
    void meshBalancing(MESH3D & grid3d) const;
    
    /*! The mesh should be in local numbering, the not-owned elements are discarded */
    void meshBalancing(const Teuchos::RCP<MESH3D> & Grid3d, const sVect<UInt> & elWeights) const;
    void meshBalancing(MESH3D & grid3d, const sVect<UInt> & elWeights) const;
    
    /*! The mesh should be in local numbering. The appendended elements and nodes are listed below the original ones */
    void meshOverlapping(const Teuchos::RCP<MESH3D> & Grid3d) const;
    void meshOverlapping(MESH3D & grid3d) const;
    
    /*! Reduces the mesh to a single process and pushes the \c cids to the global numbering.
    Duplicated elements or points are eliminated. */
    void gather(const UInt & gatherGid, const Teuchos::RCP<MESH3D> & targetGrid3d, const Teuchos::RCP<const MESH3D> & sourceGrid3d);
    void gather(const UInt & gatherGid, MESH3D & targetGrid3d, const MESH3D & sourceGrid3d);
    
    /*! Eliminates the overlapped regions. ATTENTION: since there is no way history of the way
    the mesh was overlapped the shared elements are put in a particular domain in an
    un-predictable manner */
    void destroyOverlap(const Teuchos::RCP<MESH3D> & grid3d) const;
    void destroyOverlap(MESH3D & grid3d) const;
    //@}
    
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MESH3D>       & OldMesh,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MESH3D>             & NewMesh);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MESH3D       & OldMesh,
                            const communicator & NewCommDev,
                                  MESH3D       & NewMesh);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MESH3D>       & OldMesh,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MESH3D>             & NewMesh);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MESH3D       & OldMesh,
                            const communicator & NewCommDev,
                                  MESH3D       & NewMesh);
    //@}
};


template<typename GEOSHAPE, typename NODEMAP>
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
mesh3dGlobalManip()
{
  assert(GEOSHAPE::nDim == 3);
  commDevLoaded = false;
}

template<typename GEOSHAPE, typename NODEMAP>
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
mesh3dGlobalManip(const Teuchos::RCP<const communicator> & CommDev)
{
  assert(GEOSHAPE::nDim == 3);
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename NODEMAP>
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
mesh3dGlobalManip(const communicator & CommDev)
{
  assert(GEOSHAPE::nDim == 3);
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename NODEMAP>
bool
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
check(const Teuchos::RCP<const MESH3D> & Grid3d) const
{  
  return(check(*Grid3d));
}

template<typename GEOSHAPE, typename NODEMAP>
bool
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
check(const MESH3D & Grid3d) const
{
  assert(Grid3d.getMapTransferred());
  
  //Allocations
  bool flag;
  UInt pid  = commDev->rank();
  UInt size = commDev->size();
  
  //Nodes vect check
  pVectGlobalManip<point3d,NODEMAP> nodesChecking(commDev);
  bool flagN = nodesChecking.check(Grid3d.getNodes());
  
  if(!flagN)
  { cout << "ERROR: Nodes check failed. Pid: " << pid << endl; }
  
  //Elements check
  pGraphGlobalManip<GEOELEMENT,pMapItem,NODEMAP> elementsChecking(commDev);
  bool flagE = elementsChecking.check(Grid3d.getElements());
  
  if(!flagE)
  { cout << "ERROR: Elements check failed. Pid: " << pid << endl; }
  
  //Faces check
  pGraphGlobalManip<GEOELEMENT2D,pMapItem,NODEMAP> facesChecking(commDev);
  bool flagF = facesChecking.check(Grid3d.getFaces());
  
  if(!flagF)
  { cout << "ERROR: Faces check failed. Pid: " << pid << endl; }
  
  //Edges check
  pGraphGlobalManip<GEOELEMENT1D,pMapItem,NODEMAP> edgesChecking(commDev);
  bool flagD = edgesChecking.check(Grid3d.getEdges());
  
  if(!flagD)
  { cout << "ERROR: Edges check failed. Pid: " << pid << endl; }
  
  flag = flagN & flagE & flagF & flagD;
  
  //Cross map checking
  pMapManip<NODEMAP> crossMapChecker;
  
  //Checking against elements col map
  sVect<UInt> allFlags(size);
  UInt myFlag = crossMapChecker.isEqual(Grid3d.getNodes().getMapRef(),Grid3d.getElements().getColMap());
  
  if(myFlag != 1)
  { cout << "ERROR: Elements col map differs from Nodes col map. Pid: " << pid << endl; }
  
  all_gather(*commDev, myFlag, allFlags);
  
  for(UInt i=1; i <= size; ++i)
  { myFlag = myFlag & allFlags(i); }
  
  flag = flag & myFlag;
  
  //Checking against faces col map
  if(getNumGlobalFaces(Grid3d) != 0)
  { myFlag = crossMapChecker.isEqual(Grid3d.getNodes().getMapRef(),Grid3d.getFaces().getColMap()); }
  else
  { myFlag = true; }
    
  all_gather(*commDev, myFlag, allFlags);
  
  for(UInt i=1; i <= size; ++i)
  { myFlag = myFlag & allFlags(i); }
  
  flag = flag & myFlag;
  
  if(myFlag != 1)
  { cout << "ERROR: Faces col map differs from Nodes col map. Pid: " << pid << endl; }
  
  //Checking against edges col map
  if(getNumGlobalEdges(Grid3d) != 0)
  { myFlag = crossMapChecker.isEqual(Grid3d.getNodes().getMapRef(),Grid3d.getEdges().getColMap()); }
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
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalVertices(const Teuchos::RCP<const MESH3D> & Grid3d) const
{
  pVectGlobalManip<point3d,NODEMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d->getNodes()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalVertices(const MESH3D & Grid3d) const
{
  pVectGlobalManip<point3d,NODEMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d));
}    

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalEdges(const Teuchos::RCP<const MESH3D> & Grid3d) const
{
  pVectGlobalManip<GEOELEMENT1D,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d->getEdges()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalEdges(const MESH3D & Grid3d) const
{
  pVectGlobalManip<GEOELEMENT1D,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d.getEdges()));
}    

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalFaces(const Teuchos::RCP<const MESH3D> & Grid3d) const
{
  pVectGlobalManip<GEOELEMENT2D,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d->getFaces()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalFaces(const MESH3D & Grid3d) const
{
  pVectGlobalManip<GEOELEMENT2D,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d.getFaces()));
}
    
template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalElements(const Teuchos::RCP<const MESH3D> & Grid3d) const
{
  pVectGlobalManip<GEOELEMENT,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d->getElements()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
getNumGlobalElements(const MESH3D & Grid3d) const
{
  pVectGlobalManip<GEOELEMENT,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d.getElements()));
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshPartition(const Teuchos::RCP<MESH3D> & Grid3d) const
{
  assert(commDevLoaded);
  meshPartition(*Grid3d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshPartition(mesh3d<GEOSHAPE,ELMAP,NODEMAP> & Grid3d) const
{
  //Asserts
  assert(commDevLoaded);
  assert(Grid3d.getElements().colIsLocal() == false);
  assert(Grid3d.getMeshStandard() == STDL);
  
  //Clearing
  Grid3d.clearFaces();
  Grid3d.clearEdges();
  Grid3d.transferMap();
  
  //PARMETIS_______________________________________________________________________________________
  UInt commSize = commDev->size();
  
  //Sizes communication
  UInt mySize = Grid3d.getNumElements();
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
      node    = Grid3d.getElementL(i).getCid(j) - 1; //Change numbering      
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
  idx_t ncommonnodes = GEOSHAPE2D::numPoints;
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
  
  GRAPH elGraph = Grid3d.getElements();
    
  //Downloading the data on the graph map
  for(UInt i=1; i <= mySize; ++i)
  {
    elGraph.getRowMapL(i).setPid(part[i-1]);
  }
  
  //Clearing
  delete[] elmdist;
  delete[] eptr;
  delete[] eind;
  delete[] tpwgts;
  delete[] options;
  delete[] part;
  
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
  Grid3d.setElements(elGraph);

  //Change the map to the nodes
  pVect<point3d,NODEMAP> nodesVect = Grid3d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Insert the nodes
  Grid3d.setNodes(nodesVect);
  
  //Internal methods
  Grid3d.computeNumVertices();
  Grid3d.transferMap();
  Grid3d.setMeshStandard(STDB);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshBalancing(const Teuchos::RCP<MESH3D> & Grid3d)  const
{
  assert(commDevLoaded);
  meshBalancing(*Grid3d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshBalancing(mesh3d<GEOSHAPE,ELMAP,NODEMAP> & Grid3d) const
{
  assert(commDevLoaded);
  assert(Grid3d.getElements().colIsLocal());
  assert(Grid3d.getMeshStandard() == STDB);
  
  //Clearing
  Grid3d.clearFaces();
  Grid3d.clearEdges();
  Grid3d.transferMap();
  
  //PARMETIS_______________________________________________________________________________________
  UInt commSize = commDev->size();
  
  //Sizes communication
  UInt mySize = Grid3d.getNumElements();
  sVect<UInt> allSizes(commSize);
  
  all_gather(*commDev, mySize, allSizes);
  
  //Elements distribution
  idx_t * elmdist = new idx_t[commSize + 1];
  
  elmdist[0] = 0;
  
  for(UInt i=1 ; i <= commSize; ++i)
  { elmdist[i] = elmdist[i-1] + allSizes(i); }  
  
  //Element pointer
  idx_t * eptr = new idx_t[mySize + 1];
  
  eptr[0] = 0;
  
  for(UInt i=1 ; i <= mySize; ++i)
  { eptr[i] = eptr[i-1] + GEOSHAPE::numPoints; }
  
  //Element indices
  UInt  node, k = 0;
  idx_t * eind = new idx_t[mySize * GEOSHAPE::numPoints];
  
  for(UInt i=1; i <= mySize; ++i)
  {
    for(UInt j=1; j <= GEOSHAPE::numPoints; ++j)
    {
      node    = Grid3d.getElements().getCid_LG(i,j) - 1; //Change numbering      
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
  idx_t ncommonnodes = GEOSHAPE2D::numPoints;
  idx_t nparts       = commSize;
  
  //Weights
  real_t * tpwgts = new real_t[ncon * nparts];
  
  for(int i=0; i < nparts; ++i)
  { tpwgts[i] = Real(1.0) / Real(nparts); }
  
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
  
  GRAPH elGraph = Grid3d.getElements();
    
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
  Grid3d.clearFaces();
  Grid3d.clearEdges();  
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
  Grid3d.setElements(elGraph);

  //Change the map to the nodes
  pVect<point3d,NODEMAP> nodesVect = Grid3d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Insert the nodes
  Grid3d.setNodes(nodesVect);
  
  //Internal methods
  Grid3d.computeNumVertices();
  Grid3d.transferMap();
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshBalancing(const Teuchos::RCP<MESH3D> & Grid3d, const sVect<UInt> & elWeights) const
{
  assert(Grid3d->getMeshStandard() == STDB);
  meshBalancing(*Grid3d,elWeights);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshBalancing(MESH3D & Grid3d, const sVect<UInt> & elWeights) const
{
  assert(commDevLoaded);
  assert(Grid3d.getElements().colIsLocal());
  assert(Grid3d.getMeshStandard() == STDB);
  assert(elWeights.size() == Grid3d.getNumElements());
  
  //Clearing
  Grid3d.clearFaces();
  Grid3d.clearEdges();
  Grid3d.transferMap();
  
  //PARMETIS_______________________________________________________________________________________
  UInt commSize = commDev->size();
  
  //Sizes communication
  UInt mySize = Grid3d.getNumElements();
  sVect<UInt> allSizes(commSize);
  
  all_gather(*commDev, mySize, allSizes);
  
  //Elements distribution
  idx_t * elmdist = new idx_t[commSize + 1];
  
  elmdist[0] = 0;
  
  for(UInt i=1 ; i <= commSize; ++i)
  { elmdist[i] = elmdist[i-1] + allSizes(i); }  
  
  //Element pointer
  idx_t * eptr = new idx_t[mySize + 1];
  
  eptr[0] = 0;
  
  for(UInt i=1 ; i <= mySize; ++i)
  { eptr[i] = eptr[i-1] + GEOSHAPE::numPoints; }
  
  //Element indices
  UInt  node, k = 0;
  idx_t * eind = new idx_t[mySize * GEOSHAPE::numPoints];
  
  for(UInt i=1; i <= mySize; ++i)
  {
    for(UInt j=1; j <= GEOSHAPE::numPoints; ++j)
    {
      node    = Grid3d.getElements().getCid_LG(i,j) - 1; //Change numbering      
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
  idx_t ncommonnodes = GEOSHAPE2D::numPoints;
  idx_t nparts       = commSize;
  
  //Weights
  real_t * tpwgts = new real_t[ncon * nparts];
  
  for(int i=0; i < nparts; ++i)
  { tpwgts[i] = Real(1.0) / Real(nparts); }
  
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
  
  GRAPH elGraph = Grid3d.getElements();
    
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
  delete elmwgt;
  
  
  //COMMUNICATION ACROSS THE COMUNICATOR___________________________________________________________
  Grid3d.clearFaces();
  Grid3d.clearEdges();  
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
  Grid3d.setElements(elGraph);

  //Change the map to the nodes
  pVect<point3d,NODEMAP> nodesVect = Grid3d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Insert the nodes
  Grid3d.setNodes(nodesVect);
  
  //Internal methods
  Grid3d.computeNumVertices();
  Grid3d.transferMap();
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshOverlapping(const Teuchos::RCP<MESH3D> & Grid3d) const
{
  assert(commDevLoaded);
  meshOverlapping(*Grid3d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
meshOverlapping(mesh3d<GEOSHAPE,ELMAP,NODEMAP> & Grid3d) const
{
  assert(commDevLoaded);
  assert(Grid3d.getElements().colIsLocal());
  assert(Grid3d.getMeshStandard() == STDB);
  
  //Clearing
  Grid3d.clearFaces();
  Grid3d.clearEdges();
  Grid3d.transferMap();
  
  //Elements overlapping
  typedef pGraph<GEOELEMENT,ELMAP,NODEMAP>  GRAPH;
  
  GRAPH elGraph = Grid3d.getElements();
  
  pGraphGlobalManip<GEOELEMENT,pMapItem,NODEMAP> elGraphManip(commDev);
  elGraphManip.overlap(elGraph);
  
  //Nodes overlapping
  pMap<NODEMAP>          colMap    = elGraph.getColMap();
  pVect<point3d,NODEMAP> nodesVect = Grid3d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Grid fill
  Grid3d.setElements(elGraph);
  Grid3d.setNodes(nodesVect);
  Grid3d.transferMap();
  Grid3d.setMeshStandard(STDA);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
gather(const UInt & gatherGid, const Teuchos::RCP<MESH3D> & targetGrid3d, const Teuchos::RCP<const MESH3D> & sourceGrid3d)
{
  assert(commDevLoaded);
  gather(gatherGid,*targetGrid3d,*sourceGrid3d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
gather(const UInt & gatherGid, MESH3D & targetGrid3d, const MESH3D & sourceGrid3d)
{
  assert(commDevLoaded);
  
  typedef typename MESH3D::NODESVECT    NODESVECT;
  typedef typename MESH3D::GRAPH3D      GRAPH3D;
  typedef typename MESH3D::GEOELEMENT3D GEOELEMENT3D;
  
  //Mesh cheking
  assert(commDevLoaded);
  assert(sourceGrid3d.getElements().colIsLocal());
  
  //Nodes and graph extraction
  NODESVECT nodes    = sourceGrid3d.getNodes();
  GRAPH3D   elements = sourceGrid3d.getElements();
  
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

  pVectComm<GEOELEMENT3D,ELMAP> elementsComm(commDev);
  elementsComm.vectorPid(elements);
  
  
  //Nodes - Elements exclusive merging
  pVectManip<point3d,NODEMAP> nodesManip;
  nodesManip.unionExclusive(nodes);  
  
  pVectManip<GEOELEMENT3D,ELMAP> elementsManip;
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
  
  targetGrid3d.clear();
  targetGrid3d.setNodes(nodes);
  targetGrid3d.setElements(elements);
  targetGrid3d.transferMap();
  targetGrid3d.setMeshStandard(STDU);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
destroyOverlap(const Teuchos::RCP<MESH3D> & grid3d) const
{
  assert(commDevLoaded);
  destroyOverlap(*grid3d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
destroyOverlap(MESH3D & grid3d) const
{
  assert(commDevLoaded);
  assert(grid3d.getMeshStandard() == STDA);
  
  //Mesh cheking-----------------------------------------------------------------------------------
  assert(commDevLoaded);
  assert(grid3d.getElements().colIsLocal());
  
  //Typedefs
  typedef typename MESH3D::GRAPH3D    GRAPH3D;
  typedef typename MESH3D::NODESVECT  NODESVECT;
  
  //Extract data
  GRAPH3D elements = grid3d.getElements();
  NODESVECT  nodes = grid3d.getNodes();
  
  
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
  GRAPH3D  newElements;
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
    for(UInt j=1; j <= grid3d.getNumPoints(); ++j)  //I points sono locali
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
    
    for(UInt j=1; j <= grid3d.getNumPoints(); ++j)
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
  
  grid3d.clearFaces();
  grid3d.clearEdges();
  
  grid3d.setNodes(newNodes);
  grid3d.setElements(newElements);
  grid3d.computeNumVertices();
  grid3d.transferMap();
  grid3d.setMeshStandard(STDB);
  
  //Final check------------------------------------------------------------------------------------
  assert(check(grid3d));
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MESH3D>       & OldMesh,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MESH3D>             & NewMesh)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldMesh,
                    *NewCommDev,
                    *NewMesh);
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
reduceCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MESH3D       & OldMesh,
                   const communicator & NewCommDev,
                         MESH3D       & NewMesh)
{
  //Typedefs---------------------------------------------------------
  typedef typename MESH3D::NODESVECT NODESVECT;
  typedef typename MESH3D::BOOLVECT  BOOLVECT;
  typedef typename MESH3D::GRAPH3D   GRAPH3D;
  typedef typename MESH3D::GRAPH2D   GRAPH2D;
  typedef typename MESH3D::GRAPH1D   GRAPH1D;
  
  typedef pVectGlobalManip<point3d,NODEMAP>                MANIP_POINT;
  typedef pVectGlobalManip<UInt,NODEMAP>                   MANIP_ISVERTEX;
  typedef pGraphGlobalManip<GEOELEMENT,pMapItem,NODEMAP>   MANIP_ELEMENT;
  typedef pGraphGlobalManip<GEOELEMENT2D,pMapItem,NODEMAP> MANIP_FACE;
  typedef pGraphGlobalManip<GEOELEMENT1D,pMapItem,NODEMAP> MANIP_EDGE;
  
  //Manipulators-----------------------------------------------------
  NODESVECT newNodes;
  BOOLVECT  newIsVertex;
  GRAPH3D   newElements;
  GRAPH2D   newFaces;
  GRAPH1D   newEdges;
  
  MANIP_POINT    manipPoint;
  MANIP_ISVERTEX manipIsVertex;
  MANIP_ELEMENT  manipElement;
  MANIP_FACE     manipFace;
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
    
    manipFace.reduceCommunicator(isActive,
                                 OldCommDev,
                                 OldMesh.getFaces(),
                                 NewCommDev,
                                 newFaces);
    
    manipEdge.reduceCommunicator(isActive,
                                 OldCommDev,
                                 OldMesh.getEdges(),
                                 NewCommDev,
                                 newEdges);
    
    //Fix the mesh
    NewMesh.setNodes(newNodes,newIsVertex);
    NewMesh.setElements(newElements);
    NewMesh.setFaces(newFaces);
    NewMesh.setEdges(newEdges);
    
    NewMesh.computeNumVertices();
    NewMesh.transferMap();
  }
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MESH3D>       & OldMesh,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MESH3D>             & NewMesh)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldMesh,
                    *NewCommDev,
                    *NewMesh);
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItem,NODEMAP>::
expandCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MESH3D       & OldMesh,
                   const communicator & NewCommDev,
                         MESH3D       & NewMesh)
{
  //Typedefs---------------------------------------------------------
  typedef typename MESH3D::NODESVECT NODESVECT;
  typedef typename MESH3D::BOOLVECT  BOOLVECT;
  typedef typename MESH3D::GRAPH3D   GRAPH3D;
  typedef typename MESH3D::GRAPH2D   GRAPH2D;
  typedef typename MESH3D::GRAPH1D   GRAPH1D;
  
  typedef pVectGlobalManip<point3d,NODEMAP>                MANIP_POINT;
  typedef pVectGlobalManip<UInt,NODEMAP>                   MANIP_ISVERTEX;
  typedef pGraphGlobalManip<GEOELEMENT,pMapItem,NODEMAP>   MANIP_ELEMENT;
  typedef pGraphGlobalManip<GEOELEMENT2D,pMapItem,NODEMAP> MANIP_FACE;
  typedef pGraphGlobalManip<GEOELEMENT1D,pMapItem,NODEMAP> MANIP_EDGE;
  
  //Manipulators-----------------------------------------------------
  NODESVECT newNodes;
  BOOLVECT  newIsVertex;
  GRAPH3D   newElements;
  GRAPH2D   newFaces;
  GRAPH1D   newEdges;
  
  MANIP_POINT    manipPoint;
  MANIP_ISVERTEX manipIsVertex;
  MANIP_ELEMENT  manipElement;
  MANIP_FACE     manipFace;
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
    
    manipFace.expandCommunicator(isActive,
                                  OldCommDev,
                                  OldMesh.getFaces(),
                                  NewCommDev,
                                  newFaces);
    
    manipEdge.expandCommunicator(isActive,
                                 OldCommDev,
                                 OldMesh.getEdges(),
                                 NewCommDev,
                                 newEdges);
    
    //Fix the mesh
    NewMesh.setNodes(newNodes,newIsVertex);
    NewMesh.setElements(newElements);
    NewMesh.setFaces(newFaces);
    NewMesh.setEdges(newEdges);
  }
  
  NewMesh.computeNumVertices();
  NewMesh.transferMap();
}



//_________________________________________________________________________________________________
// PMAPITEMSHARE
//-------------------------------------------------------------------------------------------------

/*! Global mesh manipulator 3d - \c pMapItemShare */
template<typename GEOSHAPE, typename NODEMAP>
class mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItemShare ELMAP;
    typedef mesh3d<GEOSHAPE,ELMAP,NODEMAP> MESH3D;
    
    typedef typename GEOSHAPE::GEOBSHAPE   GEOSHAPE2D;
    typedef typename GEOSHAPE2D::GEOBSHAPE GEOSHAPE1D;
    
    typedef geoElement<GEOSHAPE>   GEOELEMENT;
    typedef geoElement<GEOSHAPE2D> GEOELEMENT2D;
    typedef geoElement<GEOSHAPE1D> GEOELEMENT1D;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
  
    /*! @name Constructors */ //@{
  public:
    mesh3dGlobalManip();
    mesh3dGlobalManip(const Teuchos::RCP<const communicator> & CommDev);
    mesh3dGlobalManip(const communicator & CommDev);
    //@}
    
    /*! @name Information functions */ //@{
  public:
    /*! Check the graphs, the equality of the col maps of the elements, faces and edges with the map of the points */
    bool check(const Teuchos::RCP<const MESH3D> & Grid3d) const;
    bool check(const MESH3D & Grid3d) const;
    
    /*! Return the total number of vertices on the processors */
    UInt getNumGlobalVertices(const Teuchos::RCP<const MESH3D> & Grid3d) const;
    UInt getNumGlobalVertices(const MESH3D & Grid3d) const;
    
    /*! Return the total numer of edges on the processors */
    UInt getNumGlobalEdges(const Teuchos::RCP<const MESH3D> & Grid3d) const;
    UInt getNumGlobalEdges(const MESH3D & Grid3d) const;
    
    /*! Return the total numer of faces on the processors */
    UInt getNumGlobalFaces(const Teuchos::RCP<const MESH3D> & Grid3d) const;
    UInt getNumGlobalFaces(const MESH3D & Grid3d) const;
    
    /*! Return the total numer of elements on the processors */
    UInt getNumGlobalElements(const Teuchos::RCP<const MESH3D> & Grid3d) const;
    UInt getNumGlobalElements(const MESH3D & Grid3d) const;
    //@}
    
    /*! @name Operation functions */ //@{
  public:
    /*! The mesh should be in global numbering and a new mesh in local numbering is built */
    void meshPartition(const Teuchos::RCP<MESH3D> & Grid3d) const;
    void meshPartition(MESH3D & Grid3d) const;
    
    /*! The mesh should be in local numbering, the not-owned elements are discarded */
    void meshBalancing(const Teuchos::RCP<MESH3D> & Grid3d) const;
    void meshBalancing(MESH3D & Grid3d) const;
    
    /*! The mesh should be in local numbering, the not-owned elements are discarded */
    void meshBalancing(const Teuchos::RCP<MESH3D> & Grid3d, const sVect<UInt> & elWeights) const;
    void meshBalancing(MESH3D & Grid3d, const sVect<UInt> & elWeights) const;
    
    /*! The mesh should be in local numbering. The appendended elements and nodes are listed below the original ones */
    void meshOverlapping(const Teuchos::RCP<MESH3D> & Grid3d) const;
    void meshOverlapping(MESH3D & Grid3d) const;
    
    /*! Reduces the mesh to a single process and pushes the \c cids to the global numbering.
    Duplicated elements or points are eliminated. */
    void gather(const UInt & gatherGid, const Teuchos::RCP<MESH3D> & targetGrid3d, const Teuchos::RCP<const MESH3D> & sourceGrid3d);
    void gather(const UInt & gatherGid, MESH3D & targetGrid3d, const MESH3D & sourceGrid3d);
    
    /*! Eliminates the overlapped regions. The shared elements are eliminated */
    void destroyOverlap(const Teuchos::RCP<MESH3D> & grid3d) const;
    void destroyOverlap(MESH3D & grid3d) const;
    //@}
    
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MESH3D>       & OldMesh,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MESH3D>             & NewMesh);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MESH3D       & OldMesh,
                            const communicator & NewCommDev,
                                  MESH3D       & NewMesh);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MESH3D>       & OldMesh,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MESH3D>             & NewMesh);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MESH3D       & OldMesh,
                            const communicator & NewCommDev,
                                  MESH3D       & NewMesh);
    //@}
};


template<typename GEOSHAPE, typename NODEMAP>
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
mesh3dGlobalManip()
{
  assert(GEOSHAPE::nDim == 3);
  commDevLoaded = false;
}

template<typename GEOSHAPE, typename NODEMAP>
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
mesh3dGlobalManip(const Teuchos::RCP<const communicator> & CommDev)
{
  assert(GEOSHAPE::nDim == 3);
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename NODEMAP>
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
mesh3dGlobalManip(const communicator & CommDev)
{
  assert(GEOSHAPE::nDim == 3);
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename NODEMAP>
bool
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
check(const Teuchos::RCP<const MESH3D> & Grid3d) const
{
  return(check(*Grid3d));
}

template<typename GEOSHAPE, typename NODEMAP>
bool
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
check(const MESH3D & Grid3d) const
{
  assert(Grid3d.getMapTransferred());
  
  //Allocations
  bool flag;
  UInt pid  = commDev->rank();
  UInt size = commDev->size();
  
  //Nodes vect check
  pVectGlobalManip<point3d,NODEMAP> nodesChecking(commDev);
  bool flagN = nodesChecking.check(Grid3d.getNodes());
  
  if(!flagN)
  { cout << "ERROR: Nodes check failed. Pid: " << pid << endl; }
  
  //Elements check
  pGraphGlobalManip<GEOELEMENT,pMapItemShare,NODEMAP> elementsChecking(commDev);
  bool flagE = elementsChecking.check(Grid3d.getElements());
  
  if(!flagE)
  { cout << "ERROR: Elements check failed. Pid: " << pid << endl; }
  
  //Faces check
  pGraphGlobalManip<GEOELEMENT2D,pMapItemShare,NODEMAP> facesChecking(commDev);
  bool flagF = facesChecking.check(Grid3d.getFaces());
  
  if(!flagF)
  { cout << "ERROR: Faces check failed. Pid: " << pid << endl; }
  
  //Edges check
  pGraphGlobalManip<GEOELEMENT1D,pMapItemShare,NODEMAP> edgesChecking(commDev);
  bool flagD = edgesChecking.check(Grid3d.getEdges());
  
  if(!flagD)
  { cout << "ERROR: Edges check failed. Pid: " << pid << endl; }
  
  flag = flagN & flagE & flagF & flagD;
  
  //Cross map checking
  pMapManip<NODEMAP> crossMapChecker;
  
  //Checking against elements col map
  sVect<UInt> allFlags(size);
  UInt myFlag = crossMapChecker.isEqual(Grid3d.getNodes().getMapRef(),Grid3d.getElements().getColMap());
  
  if(myFlag != 1)
  { cout << "ERROR: Elements col map differs from Nodes col map. Pid: " << pid << endl; }
  
  all_gather(*commDev, myFlag, allFlags);
  
  for(UInt i=1; i <= size; ++i)
  { myFlag = myFlag & allFlags(i); }
  
  flag = flag & myFlag;
  
  //Checking against faces col map
  if(getNumGlobalFaces(Grid3d) != 0)
  {
    myFlag = crossMapChecker.isEqual(Grid3d.getNodes().getMapRef(),Grid3d.getFaces().getColMap());
    
    all_gather(*commDev, myFlag, allFlags);
  
    for(UInt i=1; i <= size; ++i)
    { myFlag = myFlag & allFlags(i); }
  
    flag = flag & myFlag;
    
    if(myFlag != 1)
    { cout << "ERROR: Faces col map differs from Nodes col map. Pid: " << pid << endl; }
  }
  
  //Checking against edges col map
  if(getNumGlobalEdges(Grid3d) != 0)
  {
    myFlag = crossMapChecker.isEqual(Grid3d.getNodes().getMapRef(),Grid3d.getEdges().getColMap());
    
    all_gather(*commDev, myFlag, allFlags);
  
    for(UInt i=1; i <= size; ++i)
    { myFlag = myFlag & allFlags(i); }
  
    flag = flag & myFlag;
    
    if(myFlag != 1)
    { cout << "ERROR: Edges col map differs from Nodes col map. Pid: " << pid << endl; }
  }
  
  
  //Cross checking for ownership - elements
  auxiliaryGraphCheck3d<GEOSHAPE,ELMAP,NODEMAP> elOwCheker;
  bool flagWE = elOwCheker.check(Grid3d.getElements());
  
  flag = flag && flagWE;
  
  if(!flagWE)
  { cout << "ERROR: Elements cross ownership failed. Pid: " << pid << endl; }
  
  
  //Cross checking for ownership - faces
  if(Grid3d.getFaces().size() != 0)
  {
    auxiliaryGraphCheck3d<GEOSHAPE2D,ELMAP,NODEMAP> faOwCheker;
    bool flagWF = faOwCheker.check(Grid3d.getFaces());
     
    flag = flag && flagWF;
  
    if(!flagWF)
    { cout << "ERROR: Faces cross ownership failed. Pid: " << pid << endl; }
  }
  
  
  //Cross checking for ownership - faces
  if(Grid3d.getEdges().size() != 0)
  {
    auxiliaryGraphCheck3d<GEOSHAPE1D,ELMAP,NODEMAP> edOwCheker;
    bool flagWD = edOwCheker.check(Grid3d.getEdges());
     
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
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalVertices(const Teuchos::RCP<const MESH3D> & Grid3d) const
{
  pVectGlobalManip<point3d,NODEMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d->getNodes()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalVertices(const MESH3D & Grid3d) const
{
  pVectGlobalManip<point3d,NODEMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d.getNodes()));
}    

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalEdges(const Teuchos::RCP<const MESH3D> & Grid3d) const
{
  pVectGlobalManip<GEOELEMENT1D,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d->getEdges()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalEdges(const MESH3D & Grid3d) const
{
  pVectGlobalManip<GEOELEMENT1D,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d.getEdges()));
}    

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalFaces(const Teuchos::RCP<const MESH3D> & Grid3d) const
{
  pVectGlobalManip<GEOELEMENT2D,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d->getFaces()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalFaces(const MESH3D & Grid3d) const
{
  pVectGlobalManip<GEOELEMENT2D,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d.getFaces()));
}
    
template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalElements(const Teuchos::RCP<const MESH3D> & Grid3d) const
{
  pVectGlobalManip<GEOELEMENT,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d->getElements()));
}

template<typename GEOSHAPE, typename NODEMAP>
UInt
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
getNumGlobalElements(const MESH3D & Grid3d) const
{
  pVectGlobalManip<GEOELEMENT,ELMAP> manipulator(commDev);
  return(manipulator.sizeG(Grid3d.getElements()));
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshPartition(const Teuchos::RCP<MESH3D> & Grid3d) const
{
  //Asserts
  assert(commDevLoaded);
  meshPartition(*Grid3d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshPartition(MESH3D & Grid3d) const
{
  //Asserts
  assert(commDevLoaded);
  assert(Grid3d.getElements().colIsLocal() == false);
  assert(Grid3d.getMeshStandard() == STDL);
  
  //Clearing
  Grid3d.clearFaces();
  Grid3d.clearEdges();
  Grid3d.transferMap();
  
  //Sharing-owningMap______________________________________________________________________________
  typedef typename MESH3D::GRAPH3D GRAPH3D;
  
  GRAPH3D tempElements = Grid3d.getElements();
  
  for(UInt i=1; i <= tempElements.getRowMap().size(); ++i)
  {
    tempElements.getRowMap().get(i).setShared(false);
    tempElements.getRowMap().get(i).setOwned(true);
  }
  
  Grid3d.setElements(tempElements);
  Grid3d.transferMap();
  
  //PARMETIS_______________________________________________________________________________________
  UInt commSize = commDev->size();
  
  //Sizes communication
  UInt mySize = Grid3d.getNumElements();
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
      node    = Grid3d.getElementL(i).getCid(j) - 1; //Change numbering      
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
  idx_t ncommonnodes = GEOSHAPE2D::numPoints;
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
  
  GRAPH elGraph = Grid3d.getElements();
    
  //Downloading the data on the graph map
  for(UInt i=1; i <= mySize; ++i)
  {
    elGraph.getRowMapL(i).setPid(part[i-1]);
  }
  
  //Clearing
  delete[] elmdist;
  delete[] eptr;
  delete[] eind;
  delete[] tpwgts;
  delete[] options;
  delete[] part;
  
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
  Grid3d.setElements(elGraph);

  //Change the map to the nodes
  pVect<point3d,NODEMAP> nodesVect = Grid3d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Insert the nodes
  Grid3d.setNodes(nodesVect);
  
  //Internal methods
  Grid3d.computeNumVertices();
  Grid3d.transferMap();
  Grid3d.setMeshStandard(STDB);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshBalancing(const Teuchos::RCP<MESH3D> & Grid3d) const
{
  assert(Grid3d->getMeshStandard() == STDB);
  meshBalancing(*Grid3d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshBalancing(MESH3D & Grid3d) const
{
  assert(Grid3d.getMeshStandard() == STDB);
  
  typedef typename MESH3D::GRAPH3D      ELGRAPH;
  typedef typename MESH3D::GEOELEMENT3D ELEMENT;
  typedef typename MESH3D::NODESVECT    NODESVECT;
  
  //Clearing
  Grid3d.clearFaces();
  Grid3d.clearEdges();
  Grid3d.transferMap();
  
  //Grid copy
  ELGRAPH   elementsCopy = Grid3d.getElements();
  ELGRAPH   elList;
  NODESVECT nodesCopy = Grid3d.getNodes();;
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
  MESH3D tempGrid;
  
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
  Grid3d = tempGrid;
  Grid3d.transferMap();
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshBalancing(const Teuchos::RCP<MESH3D> & Grid3d, const sVect<UInt> & elWeights) const
{
  assert(Grid3d->getMeshStandard() == STDB);
  meshBalancing(*Grid3d,elWeights);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshBalancing(MESH3D & Grid3d, const sVect<UInt> & elWeights) const
{
  assert(commDevLoaded);
  assert(Grid3d.getElements().colIsLocal());
  assert(Grid3d.getMeshStandard() == STDB);
  assert(elWeights.size() == Grid3d.getNumElements());
  
  //Clearing
  Grid3d.clearFaces();
  Grid3d.clearEdges();
  Grid3d.transferMap();
  
  //PARMETIS_______________________________________________________________________________________
  UInt commSize = commDev->size();
  
  //Sizes communication
  UInt mySize = Grid3d.getNumElements();
  sVect<UInt> allSizes(commSize);
  
  all_gather(*commDev, mySize, allSizes);
  
  //Elements distribution
  idx_t * elmdist = new idx_t[commSize + 1];
  
  elmdist[0] = 0;
  
  for(UInt i=1 ; i <= commSize; ++i)
  { elmdist[i] = elmdist[i-1] + allSizes(i); }  
  
  //Element pointer
  idx_t * eptr = new idx_t[mySize + 1];
  
  eptr[0] = 0;
  
  for(UInt i=1 ; i <= mySize; ++i)
  { eptr[i] = eptr[i-1] + GEOSHAPE::numPoints; }
  
  //Element indices
  UInt  node, k = 0;
  idx_t * eind = new idx_t[mySize * GEOSHAPE::numPoints];
  
  for(UInt i=1; i <= mySize; ++i)
  {
    for(UInt j=1; j <= GEOSHAPE::numPoints; ++j)
    {
      node    = Grid3d.getElements().getCid_LG(i,j) - 1; //Change numbering      
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
  idx_t ncommonnodes = GEOSHAPE2D::numPoints;
  idx_t nparts       = commSize;
  
  //Weights
  real_t * tpwgts = new real_t[ncon * nparts];
  
  for(int i=0; i < nparts; ++i)
  { tpwgts[i] = Real(1.0) / Real(nparts); }
  
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
  
  GRAPH elGraph = Grid3d.getElements();
    
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
  delete elmwgt;
  
  
  //COMMUNICATION ACROSS THE COMUNICATOR___________________________________________________________
  Grid3d.clearFaces();
  Grid3d.clearEdges();  
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
  Grid3d.setElements(elGraph);

  //Change the map to the nodes
  pVect<point3d,NODEMAP> nodesVect = Grid3d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Insert the nodes
  Grid3d.setNodes(nodesVect);
  
  //Internal methods
  Grid3d.computeNumVertices();
  Grid3d.transferMap();
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshOverlapping(const Teuchos::RCP<MESH3D> & Grid3d) const
{
  assert(commDevLoaded);
  meshOverlapping(*Grid3d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
meshOverlapping(MESH3D & Grid3d) const
{
  assert(commDevLoaded);
  assert(Grid3d.getElements().colIsLocal());
  assert(Grid3d.getMeshStandard() == STDB);
  
  //Clearing
  Grid3d.clearFaces();
  Grid3d.clearEdges();
  Grid3d.transferMap();
  
  //Elements overlapping
  typedef pGraph<GEOELEMENT,ELMAP,NODEMAP>  GRAPH;
  
  GRAPH elGraph = Grid3d.getElements();
  
  pGraphGlobalManip<GEOELEMENT,pMapItemShare,NODEMAP> elGraphManip(commDev);
  elGraphManip.overlap(elGraph);
  
  //Nodes overlapping
  pMap<NODEMAP>          colMap    = elGraph.getColMap();
  pVect<point3d,NODEMAP> nodesVect = Grid3d.getNodes();
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(nodesVect,colMap);
  
  //Grid fill
  Grid3d.setElements(elGraph);
  Grid3d.setNodes(nodesVect);
  Grid3d.transferMap();
  Grid3d.setMeshStandard(STDA);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
gather(const UInt & gatherGid, const Teuchos::RCP<MESH3D> & targetGrid3d, const Teuchos::RCP<const MESH3D> & sourceGrid3d)
{
  assert(commDevLoaded);
  gather(gatherGid,*targetGrid3d,*sourceGrid3d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
gather(const UInt & gatherGid, MESH3D & targetGrid3d, const MESH3D & sourceGrid3d)
{
  assert(commDevLoaded);
  
  typedef typename MESH3D::NODESVECT    NODESVECT;
  typedef typename MESH3D::GRAPH3D      GRAPH3D;
  typedef typename MESH3D::GEOELEMENT3D GEOELEMENT3D;
  
  //Mesh cheking
  assert(commDevLoaded);
  assert(sourceGrid3d.getElements().colIsLocal());
  
  //Nodes and graph extraction
  NODESVECT nodes    = sourceGrid3d.getNodes();
  GRAPH3D   elements = sourceGrid3d.getElements();
  
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

  pVectComm<GEOELEMENT3D,ELMAP> elementsComm(commDev);
  elementsComm.vectorPid(elements);
  
  
  //Nodes - Elements exclusive merging
  pVectManip<point3d,NODEMAP> nodesManip;
  nodesManip.unionExclusive(nodes);  
  
  pVectManip<GEOELEMENT3D,ELMAP> elementsManip;
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
  
  targetGrid3d.clear();
  targetGrid3d.setNodes(nodes);
  targetGrid3d.setElements(elements);
  targetGrid3d.transferMap();
  targetGrid3d.setMeshStandard(STDU);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
destroyOverlap(const Teuchos::RCP<MESH3D> & grid3d) const
{
  assert(commDevLoaded);
  destroyOverlap(*grid3d);
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
destroyOverlap(MESH3D & grid3d) const
{
  assert(commDevLoaded);
  assert(grid3d.getMeshStandard() == STDA);
  
  //Mesh cheking-----------------------------------------------------------------------------------
  assert(commDevLoaded);
  assert(grid3d.getElements().colIsLocal());
  
  //Typedefs
  typedef typename MESH3D::GRAPH3D    GRAPH3D;
  typedef typename MESH3D::NODESVECT  NODESVECT;
  
  //Extract data
  GRAPH3D elements = grid3d.getElements();
  NODESVECT  nodes = grid3d.getNodes();
  
  
  //Eliminate the unused elements------------------------------------------------------------------
  UInt tot = 1;
  GRAPH3D       newElements;
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
    for(UInt j=1; j <= grid3d.getNumPoints(); ++j)  //I points sono locali
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
    
    for(UInt j=1; j <= grid3d.getNumPoints(); ++j)
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
  
  grid3d.clearFaces();
  grid3d.clearEdges();
  
  grid3d.setNodes(newNodes);
  grid3d.setElements(newElements);
  grid3d.computeNumVertices();
  grid3d.transferMap();
  grid3d.setMeshStandard(STDB);
  
  //Final check------------------------------------------------------------------------------------
  assert(check(grid3d));
}

template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MESH3D>       & OldMesh,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MESH3D>             & NewMesh)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldMesh,
                    *NewCommDev,
                    *NewMesh);
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
reduceCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MESH3D       & OldMesh,
                   const communicator & NewCommDev,
                         MESH3D       & NewMesh)
{
  //Typedefs---------------------------------------------------------
  typedef typename MESH3D::NODESVECT NODESVECT;
  typedef typename MESH3D::BOOLVECT  BOOLVECT;
  typedef typename MESH3D::GRAPH3D   GRAPH3D;
  typedef typename MESH3D::GRAPH2D   GRAPH2D;
  typedef typename MESH3D::GRAPH1D   GRAPH1D;
  
  typedef pVectGlobalManip<point3d,NODEMAP>                     MANIP_POINT;
  typedef pVectGlobalManip<UInt,NODEMAP>                        MANIP_ISVERTEX;
  typedef pGraphGlobalManip<GEOELEMENT,pMapItemShare,NODEMAP>   MANIP_ELEMENT;
  typedef pGraphGlobalManip<GEOELEMENT2D,pMapItemShare,NODEMAP> MANIP_FACE;
  typedef pGraphGlobalManip<GEOELEMENT1D,pMapItemShare,NODEMAP> MANIP_EDGE;
  
  //Manipulators-----------------------------------------------------
  NODESVECT newNodes;
  BOOLVECT  newIsVertex;
  GRAPH3D   newElements;
  GRAPH2D   newFaces;
  GRAPH1D   newEdges;
  
  MANIP_POINT    manipPoint;
  MANIP_ISVERTEX manipIsVertex;
  MANIP_ELEMENT  manipElement;
  MANIP_FACE     manipFace;
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
    
    manipFace.reduceCommunicator(isActive,
                                  OldCommDev,
                                  OldMesh.getFaces(),
                                  NewCommDev,
                                  newFaces);
    
    manipEdge.reduceCommunicator(isActive,
                                 OldCommDev,
                                 OldMesh.getEdges(),
                                 NewCommDev,
                                 newEdges);
    
    //Fix the mesh
    NewMesh.setNodes(newNodes,newIsVertex);
    NewMesh.setElements(newElements);
    NewMesh.setFaces(newFaces);
    NewMesh.setEdges(newEdges);
    
    NewMesh.computeNumVertices();
    NewMesh.transferMap();
  }
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MESH3D>       & OldMesh,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MESH3D>             & NewMesh)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldMesh,
                    *NewCommDev,
                    *NewMesh);
}
    
template<typename GEOSHAPE, typename NODEMAP>
void
mesh3dGlobalManip<GEOSHAPE,pMapItemShare,NODEMAP>::
expandCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MESH3D       & OldMesh,
                   const communicator & NewCommDev,
                         MESH3D       & NewMesh)
{
  //Typedefs---------------------------------------------------------
  typedef typename MESH3D::NODESVECT NODESVECT;
  typedef typename MESH3D::BOOLVECT  BOOLVECT;
  typedef typename MESH3D::GRAPH3D   GRAPH3D;
  typedef typename MESH3D::GRAPH2D   GRAPH2D;
  typedef typename MESH3D::GRAPH1D   GRAPH1D;
  
  typedef pVectGlobalManip<point3d,NODEMAP>                     MANIP_POINT;
  typedef pVectGlobalManip<UInt,NODEMAP>                        MANIP_ISVERTEX;
  typedef pGraphGlobalManip<GEOELEMENT,pMapItemShare,NODEMAP>   MANIP_ELEMENT;
  typedef pGraphGlobalManip<GEOELEMENT2D,pMapItemShare,NODEMAP> MANIP_FACE;
  typedef pGraphGlobalManip<GEOELEMENT1D,pMapItemShare,NODEMAP> MANIP_EDGE;
  
  //Manipulators-----------------------------------------------------
  NODESVECT newNodes;
  BOOLVECT  newIsVertex;
  GRAPH3D   newElements;
  GRAPH2D   newFaces;
  GRAPH1D   newEdges;
  
  MANIP_POINT    manipPoint;
  MANIP_ISVERTEX manipIsVertex;
  MANIP_ELEMENT  manipElement;
  MANIP_FACE     manipFace;
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
    
    manipFace.expandCommunicator(isActive,
                                  OldCommDev,
                                  OldMesh.getFaces(),
                                  NewCommDev,
                                  newFaces);
    
    manipEdge.expandCommunicator(isActive,
                                 OldCommDev,
                                 OldMesh.getEdges(),
                                 NewCommDev,
                                 newEdges);
    
    //Fix the mesh
    NewMesh.setNodes(newNodes,newIsVertex);
    NewMesh.setElements(newElements);
    NewMesh.setFaces(newFaces);
    NewMesh.setEdges(newEdges);
  }
  
  NewMesh.computeNumVertices();
  NewMesh.transferMap();
}

#endif
