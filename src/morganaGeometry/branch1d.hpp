/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef BRANCH1D_HPP
#define BRANCH1D_HPP

#include "branch1d_card.h"
#include "branch1d_searchCard.h"
#include "branch1d_bif.h"
#include "mesh1d.hpp"
#include "connect1d.hpp"


/*! The branch 1d */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class branch1d : public mesh1d<GEOSHAPE,ELMAP,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef branch1d_card          ENDCARD;
    typedef pVect<ENDCARD,NODEMAP> ENDNODES;
    typedef sVect<ENDCARD>         ENDNODEDATA;
    
    typedef branch1d_bif                   BIFURCATIONCARD;
    typedef pVect<BIFURCATIONCARD,NODEMAP> BIFURCATIONNODES;
    typedef sVect<BIFURCATIONCARD>         BIFURCATIONDATA;
    
    typedef pVect<point3d,NODEMAP> COORDVECT;
    
    typedef mesh1d<GEOSHAPE,ELMAP,NODEMAP>    MESH1D;
    typedef connect1d<GEOSHAPE,ELMAP,NODEMAP> CONNECT1D;
    typedef typename MESH1D::NODESVECT     NODESVECT;
    typedef typename MESH1D::BOOLVECT      BOOLVECT;
    typedef typename MESH1D::GRAPH1D       GRAPH1D;
    typedef typename MESH1D::GEOELEMENT1D  GEOELEMENT1D;
    //@}
    
    /*! @name Internal variables and links */ //@{
  public:
    ENDNODES         endNodes;
    BIFURCATIONNODES bifurcationNodes;
    
    Teuchos::RCP<const communicator> commDev;
    Teuchos::RCP<CONNECT1D>          connectGrid1d;
    //@}
    
    /*! @name Internal flags and structures */ //@{
  public:
    bool branchOk;
    bool commDevOk;
    bool connectOk;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    branch1d();
    branch1d(const NODESVECT & Nodes, const BOOLVECT & IsVertex, const GRAPH1D & Elements);
    branch1d(const NODESVECT & Nodes, const GRAPH1D & Elements);
    branch1d(const branch1d & Mesh);
    branch1d(const MESH1D   & Mesh);
    branch1d & operator=(const branch1d & Mesh);
    branch1d & operator=(const MESH1D   & Mesh);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    void setConnect1d(const Teuchos::RCP<CONNECT1D> & ConnectGrid1d);
    //@}
    
    /*! @name Build functions */ //@{
  public:
    void buildTriplePoints();
    void clearTriplePoints();
    void overlapGrid();
    //@}
    
    /*! @name Get functions - end nodes */ //@{
  public:
    UInt             getNumEndNodes() const;
    ENDNODES       & getEndNodes();
    const ENDNODES & getEndNodes() const;
    const point3d  & getEndNode(const UInt & i);
    UInt             getEndNodeToNodeL(const UInt & i) const;
    UInt             getEndNodeToNodeG(const UInt & i) const;
    UInt             getEndNodeToIntElementL(const UInt & i) const;
    UInt             getEndNodeToIntElementG(const UInt & i) const;
    point3d          getEndNodeToIntNormal(const UInt & i) const;
    UInt             getEndNodeToExtNodeL(const UInt & i) const;
    UInt             getEndNodeToExtNodeG(const UInt & i) const;
    UInt             getNumEndNodeToExtElements(const UInt & i) const;
    UInt             getEndNodeToExtElementsL(const UInt & i, const UInt & j) const;
    UInt             getEndNodeToExtElementsG(const UInt & i, const UInt & j) const;
    UInt             getNumEndNodeToExtNormals(const UInt & i) const;
    point3d          getEndNodeToExtNormals(const UInt & i, const UInt & j) const;
    bool             isGidEndNode(const UInt & gid) const;
    //@}
    
    /*! @name Get functions - bifurcation nodes */ //@{
  public:
    UInt                     getNumBifurcationNodes() const;
    BIFURCATIONNODES       & getBifurcationNodes();
    const BIFURCATIONNODES & getBifurcationNodes() const;
    const UInt             & getBifurcationNodeL(const UInt & i) const;
    const UInt             & getBifurcationNodeG(const UInt & i) const;
    UInt                     getBifurcationIntToNumElements(const UInt & i) const;
    const UInt             & getBifurcationToIntElementsL(const UInt & i, const UInt & j) const;
    const UInt             & getBifurcationToIntElementsG(const UInt & i, const UInt & j) const;
    UInt                     getBifurcationToNumIntNormals(const UInt & i) const;
    const point3d          & getBifurcationToIntNormals(const UInt & i, const UInt & j) const;
    const UInt             & getBifurcationToExtNodeL(const UInt & i) const;
    const UInt             & getBifurcationToExtNodeG(const UInt & i) const;
    const UInt             & getBifurcationToExtElementL(const UInt & i) const;
    const UInt             & getBifurcationToExtElementG(const UInt & i) const;
    const point3d          & getBifurcationToExtNormal(const UInt & i) const;
    bool                     isGidBifurcation(const UInt & gid) const;
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
branch1d() : mesh1d<GEOSHAPE,ELMAP,NODEMAP>()
{
  branchOk  = false;
  commDevOk = false;
  connectOk = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
branch1d(const NODESVECT & Nodes, const BOOLVECT & IsVertex, const GRAPH1D & Elements) : mesh1d<GEOSHAPE,ELMAP,NODEMAP>(Nodes,IsVertex,Elements)
{
  branchOk  = false;
  commDevOk = false;
  connectOk = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
branch1d(const NODESVECT & Nodes, const GRAPH1D & Elements) : mesh1d<GEOSHAPE,ELMAP,NODEMAP>(Nodes,Elements)
{
  branchOk  = false;
  commDevOk = false;
  connectOk = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
branch1d(const branch1d & Mesh) : mesh1d<GEOSHAPE,ELMAP,NODEMAP>(Mesh)
{
  branchOk  = Mesh.branchOk;
  commDevOk = Mesh.commDevOk;
  connectOk = Mesh.connectOk;
  
  commDev       = Mesh.commDev;
  connectGrid1d = Mesh.connectGrid1d;
  
  endNodes = Mesh.endNodes;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
branch1d(const MESH1D & Mesh) : mesh1d<GEOSHAPE,ELMAP,NODEMAP>(Mesh)
{
  branchOk  = false;
  commDevOk = false;
  connectOk = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
branch1d<GEOSHAPE,ELMAP,NODEMAP> &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
operator=(const branch1d & Mesh)
{
  MESH1D::operator=(Mesh);
    
  branchOk  = Mesh.branchOk;
  commDevOk = Mesh.commDevOk;
  connectOk = Mesh.connectOk;
  
  commDev       = Mesh.commDev;
  connectGrid1d = Mesh.connectGrid1d;
  
  endNodes = Mesh.endNodes;
  
  return *this;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
branch1d<GEOSHAPE,ELMAP,NODEMAP> &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
operator=(const MESH1D & Mesh)
{
  MESH1D::operator=(Mesh);
    
  branchOk  = false;
  commDevOk = false;
  connectOk = false;
  
  return *this;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevOk = true;
  commDev   = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const communicator & CommDev)
{
  commDevOk = true;
  commDev   = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
setConnect1d(const Teuchos::RCP<CONNECT1D> & ConnectGrid1d)
{
  connectOk = true;
  connectGrid1d = ConnectGrid1d;
}


//_________________________________________________________________________________________________
// BUILD FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
buildTriplePoints()
{  
  //Flags------------------------------------------------------------
  assert(commDevOk);
  assert(connectOk);
  
  branchOk = true;
  
  //Typedefs---------------------------------------------------------
  typedef std::multimap<point3d,UInt>          NODESTL;
  typedef std::pair<point3d,UInt>              NODESTL_DATA;
  typedef typename NODESTL::iterator           NODESTL_ITER;
  typedef std::pair<NODESTL_ITER,NODESTL_ITER> NODESTL_RANGE;
  
  typedef std::multimap<UInt,UInt>   IDSSTL;
  typedef std::pair<UInt,UInt>       IDSSTL_DATA;
  typedef typename IDSSTL::iterator  IDSSTL_ITER;
  
  typedef pVect<branch1d_searchCard,pMapItem>     PVECT_BRANCH;
  typedef pVectComm<branch1d_searchCard,pMapItem> COMM_BRANCH;
  typedef sVect<branch1d_searchCard>              DATA_BRANCH;
  typedef sVect<pMapItem>                         PMAP_BRANCH;
  
  typedef geoElement<GEOSHAPE>                          GEOELEMENT1D;
  typedef pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> ELMANIP;
  typedef typename ELMANIP::PINVGRAPH  PINVGRAPH;
  
  //Clearing vectors-------------------------------------------------
  endNodes.clear();
  
  ENDNODEDATA     endNodesData;
  BIFURCATIONDATA bifurcationNodesData;
  COORDVECT       tempCoordVect;
  
  //Alloc------------------------------------------------------------
  UInt nodeGid, nodeLid, elGid, elLid, count;
  UInt lid1, lid2, gid1, gid2;
  UInt nodeLid1, nodeLid2, nodeGid1, nodeGid2;
  UInt srcNodeGid, tgtNodeGid, tgtNodeLid;
  point3d P, P1, P2, N;
  
  NODESTL       nodeStl;
  NODESTL_DATA  nodeData;
  NODESTL_RANGE nodeRange;
  
  ENDCARD     card;
  pMapItem mapItem;
  branch1d_searchCard searchCard;
  branch1d_bif        bifurcationCard;
  
  PVECT_BRANCH pVectBranch;
  DATA_BRANCH  endNodeSearchData;
  PMAP_BRANCH  endNodeSearchPmap;
  sVect<UInt>    gidVect(commDev->size());
  sVect<point3d> intNormals(2);
  
  //Create the list of end points------------------------------------
  PINVGRAPH invGraph;
  ELMANIP   elementInverter(commDev);
  elementInverter.inversion(MESH1D::elements,invGraph);
  
  assert(MESH1D::getNumNodes() == invGraph.rowSize());
  
  for(UInt p = 1; p <= MESH1D::getNumNodes(); ++p)
  {
    if(invGraph.getItemL(p).size() == 1)
    {
      nodeGid = MESH1D::nodes.getMapL(p).getGid();
      elLid   = connectGrid1d->getVertexToElement(p,1);
      elGid   = MESH1D::elements.getRowMapL(elLid).getGid();
      
      lid1 = MESH1D::elements.getItemL(elLid).getCid(1);
      lid2 = MESH1D::elements.getItemL(elLid).getCid(2);
      
      if(lid1 == p)
      {
        P1 = MESH1D::nodes(lid2);
        P2 = MESH1D::nodes(lid1);
      }
      else
      {
        P1 = MESH1D::nodes(lid1);
        P2 = MESH1D::nodes(lid2);
      }
      
      N  = P2 - P1;
      N /= point3d::norm2(N);
      
      card.setNodeGid(nodeGid);
      card.setIntElementGid(elGid);
      card.setIntNormal(N);
      
      endNodesData.push_back(card);
    }
  }  
  
  //Create the list of all the nodes---------------------------------
  for(UInt p = 1; p <= MESH1D::getNumNodes(); ++p)
  {
    nodeData.first  = MESH1D::nodes(p);
    nodeData.second = MESH1D::nodes.getMapL(p).getGid();
    nodeStl.insert(nodeData);
  }
  
  //Search the triple points-----------------------------------------
  for(UInt i=1; i <= endNodesData.size(); ++i)
  {
    nodeGid = endNodesData(i).getNodeGid();
    P       = MESH1D::getNodeG(nodeGid);
    count   = nodeStl.count(P);
    
    if(count != 1)
    {
      nodeRange = nodeStl.equal_range(P);
      
      for(NODESTL_ITER iter = nodeRange.first; iter != nodeRange.second; ++iter)
      {
        if( (nodeGid != iter->second) && (invGraph.getItemG(iter->second).size() == 2) )
        { endNodesData(i).setExtNodeGid(iter->second); }
      }
    }
  }
  
  //Get elements-----------------------------------------------------
  for(UInt i=1; i <= endNodesData.size(); ++i)
  {
    if(endNodesData(i).getExtNodeGid() != 0)
    {
      nodeGid = endNodesData(i).getExtNodeGid();
      nodeLid = MESH1D::nodes.getMapG(nodeGid).getLid();
      
      assert(connectGrid1d->getNumVertexToElement(nodeLid) <= 2);
      
      for(UInt j=1; j <= connectGrid1d->getNumVertexToElement(nodeLid); ++j)
      {
        elLid = connectGrid1d->getVertexToElement(nodeLid,j);
        elGid = MESH1D::elements.getRowMapL(elLid).getGid();
        
        endNodesData(i).addExtElementsGid(elGid);
      }
    }
  }
  
  //Get normals------------------------------------------------------
  for(UInt i=1; i <= endNodesData.size(); ++i)
  {
    for(UInt j=1; j <= endNodesData(i).getExtElementsNum(); ++j)
    {
      elGid = endNodesData(i).getExtElementsGid(j);
      
      lid1 = MESH1D::getElementG(elGid).getCid(1);
      lid2 = MESH1D::getElementG(elGid).getCid(2);
      
      gid1 = MESH1D::nodes.getMapL(lid1).getGid();
      gid2 = MESH1D::nodes.getMapL(lid2).getGid();
      
      if(endNodesData(i).getExtNodeGid() == gid1)
      {
        P1 = MESH1D::getNodeL(lid2);
        P2 = MESH1D::getNodeL(lid1);
      }
      
      if(endNodesData(i).getExtNodeGid() == gid2)
      {
        P1 = MESH1D::getNodeL(lid1);
        P2 = MESH1D::getNodeL(lid2);
      }
      
      N  = P2 - P1;
      N /= point3d::norm2(N);
      
      endNodesData(i).addExtNormals(N);
    }
  }
  
  //Prepare parallel card--------------------------------------------
  mapItem.setPid(commDev->rank());
  endNodeSearchData.resize(endNodesData.size());
  endNodeSearchPmap.resize(endNodesData.size());

  for(UInt i=1; i <= endNodesData.size(); ++i)
  {
    mapItem.setLid(i);
    
    nodeGid = endNodesData(i).getNodeGid();
    searchCard.setTgtNodeGid(nodeGid);
    searchCard.setTgtNode(MESH1D::nodes.getDataG(nodeGid));
    searchCard.setTgtElGid(endNodesData(i).getIntElementGid());
    searchCard.setTgtNormal(endNodesData(i).getIntNormal());
    
    endNodeSearchData(i) = searchCard;
    endNodeSearchPmap(i) = mapItem;
  }
  
  //Communicate
  sVect<DATA_BRANCH> tempDataBranch(commDev->size());
  sVect<PMAP_BRANCH> tempPmapBranch(commDev->size());
  
  all_gather(*commDev, endNodeSearchData, tempDataBranch);
  all_gather(*commDev, endNodeSearchPmap, tempPmapBranch);

  //Assemble the vector
  UInt localSize   = 0;
  UInt globalStray = 0;
  
  for(UInt i=1; i <= UInt(commDev->size()); ++i)
  {
    if( int(i-1) != commDev->rank())
    { localSize += tempDataBranch(i).size(); }
  }
  
  pVectBranch.reserve(localSize);
  
  for(UInt i=1; i <= UInt(commDev->size()); ++i)
  {
    if( int(i-1) != commDev->rank())
    {
      for(UInt j=1; j <= tempDataBranch(i).size(); ++j)
      { pVectBranch.push_back(tempPmapBranch(i)(j),tempDataBranch(i)(j)); }
    }
  }
  
  //Fill the parallel vector-----------------------------------------
  PVECT_BRANCH newBranch;
  localSize = 0;
  
  for(UInt i=1; i <= pVectBranch.size(); ++i)
  {
    mapItem    = pVectBranch.getMapL(i);
    searchCard = pVectBranch.getDataL(i);
      
    tgtNodeGid = searchCard.getTgtNodeGid();
    P          = searchCard.getTgtNode();
    
    if(nodeStl.count(P) != 0)
    {
      nodeRange = nodeStl.equal_range(P);
      
      for(NODESTL_ITER iter = nodeRange.first; iter != nodeRange.second; ++iter)
      {
        if( (tgtNodeGid != iter->second) && (invGraph.getItemG(iter->second).size() == 2) )
        {
          localSize++;
          mapItem.setLid(localSize);
          
          srcNodeGid = iter->second;
          searchCard.setSrcNodeGid(srcNodeGid);
          
          newBranch.push_back(mapItem,searchCard);
        }
      }
    }
  }
  
  pVectBranch = newBranch;
  
  //Global numbering-------------------------------------------------
  sVect<UInt> allSizes(commDev->rank());
  all_gather(*commDev, UInt(pVectBranch.size()), allSizes);
  
  for(UInt i=1; i <= allSizes.size(); ++i)
  {
    if( int(i-1) < commDev->rank())
    { globalStray += allSizes(i); }
  }
  
  for(UInt i=1; i <= pVectBranch.size(); ++i)
  { pVectBranch.getMapL(i).setGid(i + globalStray); }
  
  //Parallel elements------------------------------------------------
  for(UInt i=1; i <= pVectBranch.size(); ++i)
  {
    nodeGid = pVectBranch(i).getSrcNodeGid();
    nodeLid = MESH1D::nodes.getMapG(nodeGid).getLid();
    
    assert(connectGrid1d->getNumVertexToElement(nodeLid) <= 2);
    
    for(UInt j=1; j <= connectGrid1d->getNumVertexToElement(nodeLid); ++j)
    {
      elLid = connectGrid1d->getVertexToElement(nodeLid,j);
      elGid = MESH1D::elements.getRowMapL(elLid).getGid();
      pVectBranch(i).getSrcElementsGids().push_back(elGid);
    }
  }
  
  //Parallel normals-------------------------------------------------
  for(UInt i=1; i <= pVectBranch.size(); ++i)
  {
    for(UInt j=1; j <= pVectBranch(i).getSrcElementsGids().size(); ++j)
    {
      elGid = pVectBranch(i).getSrcElementsGids()(j);
      
      lid1 = MESH1D::getElementG(elGid).getCid(1);
      lid2 = MESH1D::getElementG(elGid).getCid(2);
      
      gid1 = MESH1D::nodes.getMapL(lid1).getGid();
      gid2 = MESH1D::nodes.getMapL(lid2).getGid();
      
      if(pVectBranch(i).getSrcNodeGid() == gid1)
      {
        P1 = MESH1D::getNodeL(lid2);
        P2 = MESH1D::getNodeL(lid1);
      }
      
      if(pVectBranch(i).getSrcNodeGid() == gid2)
      {
        P1 = MESH1D::getNodeL(lid1);
        P2 = MESH1D::getNodeL(lid2);
      }
      
      N  = P2 - P1;
      N /= point3d::norm2(N);
      
      pVectBranch(i).getSrcNormals().push_back(N);
    }
  }
  
  //Communicate back-------------------------------------------------
  COMM_BRANCH commBranch(commDev);
  
  pVectBranch.updateFinder();
  commBranch.vectorPid(pVectBranch);
  
  //Merge in the local structure-------------------------------------
  IDSSTL      endNodeStl;
  IDSSTL_DATA endNodeStl_data;
  IDSSTL_ITER endNodeStl_iter;
  
  UInt index;
  std::set<UInt> refSet;
  
  for(UInt i=1; i <= endNodesData.size(); ++i)
  {
    endNodeStl_data.first  = endNodesData(i).getNodeGid();
    endNodeStl_data.second = i;
    
    endNodeStl.insert(endNodeStl_data);
  }
  
  for(UInt i=1; i <= pVectBranch.size(); ++i)
  {
    assert(endNodeStl.count(pVectBranch(i).getTgtNodeGid()) == 1);
    
    endNodeStl_iter = endNodeStl.find(pVectBranch(i).getTgtNodeGid());
    index           = endNodeStl_iter->second;
    
    if(endNodesData(index).getExtNodeGid() == 0)
    { endNodesData(index).setExtNodeGid(pVectBranch(i).getSrcNodeGid()); }
    else
    { assert(endNodesData(index).getExtNodeGid() == pVectBranch(i).getSrcNodeGid()); } //If fails means that two pieces of line crosses in a non-end point
    
    refSet.clear();
    for(UInt j=1; j <= endNodesData(index).getExtElementsNum(); ++j)
    { refSet.insert(endNodesData(index).getExtElementsGid(j)); }
    
    for(UInt j=1; j <= pVectBranch(i).getSrcElementsGids().size(); ++j)
    {
      if(refSet.count(pVectBranch(i).getSrcElementsGids()(j)) == 0)
      {
        endNodesData(index).addExtElementsGid(pVectBranch(i).getSrcElementsGids()(j));
        endNodesData(index).addExtNormals(pVectBranch(i).getSrcNormals()(j));
      }
    }
  }
  
  //Build end node vector--------------------------------------------
  tempCoordVect.resize(endNodesData.size());
  
  UInt pid = commDev->rank();
  NODEMAP nodeMapItem;
  nodeMapItem.setPid(pid);
  
  for(UInt i=1; i <= endNodesData.size(); ++i)
  {
    gid1 = endNodesData(i).getNodeGid();
    P1   = MESH1D::getNodeG(gid1);
    
    nodeMapItem.setLid(i);
    nodeMapItem.setGid(i);
    tempCoordVect.getMapL(i)  = nodeMapItem;
    tempCoordVect.getDataL(i) = P1;
  }
  
  pVectGlobalManip<point3d,NODEMAP> endNodeManip(commDev);
  endNodeManip.buildGlobalNumbering(tempCoordVect);
  
  //Build the final end-node vector----------------------------------
  endNodes.setData(tempCoordVect.getMapRef(),endNodesData);
  
  //Build bifurcation vect-------------------------------------------
  tempCoordVect.clear();
  gidVect.resize(2);
  index = 0;
  
  for(UInt i=1; i <= tempDataBranch.size(); ++i)
  {
    for(UInt j=1; j <= tempDataBranch(i).size(); ++j)
    {
      P         = tempDataBranch(i)(j).getTgtNode();
      nodeRange = nodeStl.equal_range(P);
      
      for(NODESTL_ITER iter = nodeRange.first; iter != nodeRange.second; ++iter)
      {
        tgtNodeGid = iter->second;
        tgtNodeLid = MESH1D::nodes.getMapG(tgtNodeGid).getLid();
        
        if(connectGrid1d->getNumVertexToElement(tgtNodeLid) >= 2)
        {
          assert(connectGrid1d->getNumVertexToElement(tgtNodeLid) == 2);
        
          bifurcationCard.setNodeGid(tgtNodeGid);
          bifurcationCard.setExtNodeGid(tempDataBranch(i)(j).getTgtNodeGid());
          bifurcationCard.setExtElementGid(tempDataBranch(i)(j).getTgtElGid());
          bifurcationCard.setExtNormal(tempDataBranch(i)(j).getTgtNormal());
          
          lid1 = connectGrid1d->getVertexToElement(tgtNodeLid,1);
          lid2 = connectGrid1d->getVertexToElement(tgtNodeLid,2);
          
          gid1 = MESH1D::elements.getMapL(lid1).getGid();
          gid2 = MESH1D::elements.getMapL(lid2).getGid();
          
          gidVect(1) = gid1;
          gidVect(2) = gid2;
          bifurcationCard.setIntElementsGid(gidVect);
          
          //Norm1
          elGid = gidVect(1);
          
          nodeLid1 = MESH1D::getElementG(elGid).getCid(1);
          nodeLid2 = MESH1D::getElementG(elGid).getCid(2);
          
          nodeGid1 = MESH1D::nodes.getMapL(nodeLid1).getGid();
          nodeGid2 = MESH1D::nodes.getMapL(nodeLid2).getGid();
          
          if(tgtNodeGid == nodeGid1)
          {
            P1 = MESH1D::getNodeL(nodeLid2);
            P2 = MESH1D::getNodeL(nodeLid1);
          }
          
          if(tgtNodeGid == nodeGid2)
          {
            P1 = MESH1D::getNodeL(nodeLid1);
            P2 = MESH1D::getNodeL(nodeLid2);
          }
          
          N  = P2 - P1;
          N /= point3d::norm2(N);
          
          intNormals(1) = N;
          
          //Norm2
          elGid = gidVect(2);
          
          nodeLid1 = MESH1D::getElementG(elGid).getCid(1);
          nodeLid2 = MESH1D::getElementG(elGid).getCid(2);
          
          nodeGid1 = MESH1D::nodes.getMapL(nodeLid1).getGid();
          nodeGid2 = MESH1D::nodes.getMapL(nodeLid2).getGid();
          
          if(tgtNodeGid == nodeGid1)
          {
            P1 = MESH1D::getNodeL(nodeLid2);
            P2 = MESH1D::getNodeL(nodeLid1);
          }
          
          if(tgtNodeGid == nodeGid2)
          {
            P1 = MESH1D::getNodeL(nodeLid1);
            P2 = MESH1D::getNodeL(nodeLid2);
          }
          
          N  = P2 - P1;
          N /= point3d::norm2(N);
          
          intNormals(2) = N;
          
          bifurcationCard.setIntNormal(intNormals);
          
          //Upload data
          bifurcationNodesData.push_back(bifurcationCard);
          
          //Upload node vect
          index++;
          nodeMapItem.setLid(index);
          nodeMapItem.setGid(index);
          tempCoordVect.push_back(nodeMapItem,P);          
        }
      }
    }
  }
  
  endNodeManip.buildGlobalNumbering(tempCoordVect);
  bifurcationNodes.setData(tempCoordVect.getMapRef(),bifurcationNodesData);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
clearTriplePoints()
{
  branchOk = false;
  endNodes.clear();
}


template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
overlapGrid()
{
  //Assert-----------------------------------------------------------
  assert(branchOk);
  
  //Alloc------------------------------------------------------------
  UInt elLid, elGid;
  set<UInt> elmSet;
  
  //Mesh overlap-----------------------------------------------------
  ELMAP elMapItem; elMapItem.setPid(commDev->rank());
  pMap<ELMAP> newElMap = MESH1D::elements.getRowMap();
  
  for(UInt i=1; i <= endNodes.size(); ++i)
  {
    for(UInt j=1; j <= endNodes(i).getExtElementsNum(); ++j)
    {
      elLid = newElMap.size() + 1;
      elGid = endNodes(i).getExtElementsGid(j);
      
      if( (!MESH1D::elements.isRowG(elGid)) && (elmSet.count(elGid) == 0) )
      {
        elmSet.insert(elGid);
        elMapItem.setLid(elLid);
        elMapItem.setGid(elGid);
        newElMap.push_back(elMapItem);
      }
    }
  }
  
  for(UInt i=1; i <= bifurcationNodes.size(); ++i)
  {
    elLid = newElMap.size() + 1;
    elGid = bifurcationNodes(i).getExtElementGid();
    
    if( (!MESH1D::elements.isRowG(elGid)) && (elmSet.count(elGid) == 0) )
    {
      elmSet.insert(elGid);
      elMapItem.setLid(elLid);
      elMapItem.setGid(elGid);
      newElMap.push_back(elMapItem);
    }
  }
  
  colMapFixer_changeMap(newElMap,*commDev);
  
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> elementsManip(commDev);
  GRAPH1D newElements = MESH1D::elements;
  elementsManip.changeRowMap(newElements,newElMap);  
  
  //Change node map--------------------------------------------------
  pMap<NODEMAP>        newNodeMap = newElements.getColMap();
  pVect<point3d,NODEMAP> newNodes = MESH1D::nodes;
  
  pVectGlobalManip<point3d,NODEMAP> nodesManip(commDev);
  nodesManip.changeMap(newNodes,newNodeMap);
  
  //Upload-----------------------------------------------------------
  MESH1D::setElements(newElements);
  MESH1D::setNodes(newNodes);
}


//_________________________________________________________________________________________________
// GET FUNCTIONS - ENDNODES
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumEndNodes() const
{
  assert(branchOk);
  return(endNodes.size());
}
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
typename branch1d<GEOSHAPE,ELMAP,NODEMAP>::ENDNODES &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getEndNodes()
{
  return(endNodes);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename branch1d<GEOSHAPE,ELMAP,NODEMAP>::ENDNODES &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getEndNodes() const
{
  return(endNodes);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getEndNode(const UInt & i)
{
  assert(i <= endNodes.size());
  
  UInt gid = endNodes(i).getNodeGid();
  
  return(MESH1D::getNodeG(gid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getEndNodeToNodeL(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= endNodes.size());
  assert(MESH1D::nodes.isG(endNodes(i).getNodeGid()));
  
  return(MESH1D::nodes.getMapG(endNodes(i).getNodeGid()).getLid());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getEndNodeToNodeG(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= endNodes.size());
  
  return(endNodes(i).getNodeGid());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getEndNodeToIntElementL(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= endNodes.size());
  assert(MESH1D::elements.isG(endNodes(i).getIntElementGid()));
  
  return(MESH1D::elements.getMapG(endNodes(i).getIntElementGid()).getLid());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getEndNodeToIntElementG(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= endNodes.size());
  
  return(endNodes(i).getIntElementGid());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
point3d
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getEndNodeToIntNormal(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= endNodes.size());
  
  return(endNodes(i).getIntNormal());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getEndNodeToExtNodeL(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= endNodes.size());
  
  return(MESH1D::nodes.getMapG(endNodes(i).getExtNodeGid()).getLid());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getEndNodeToExtNodeG(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= endNodes.size());
  
  return(endNodes(i).getExtNodeGid());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumEndNodeToExtElements(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= endNodes.size());
  
  return(endNodes(i).getExtElementsNum());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getEndNodeToExtElementsL(const UInt & i, const UInt & j) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= endNodes.size());
  assert(j >= 1);
  assert(j <= endNodes(i).getExtElementsNum());
  assert(MESH1D::elements.isG(endNodes(i).getExtElementsGid(j)));
  
  return(MESH1D::elements.getMapG(endNodes(i).getExtElementsGid(j)).getLid());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getEndNodeToExtElementsG(const UInt & i, const UInt & j) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= endNodes.size());
  assert(j >= 1);
  assert(j <= endNodes(i).getExtElementsNum());
  
  return(endNodes(i).getExtElementsGid(j));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumEndNodeToExtNormals(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= endNodes.size());
  
  return(endNodes(i).getExtNormalsNum());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
point3d
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getEndNodeToExtNormals(const UInt & i, const UInt & j) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= endNodes.size());
  assert(j >= 1);
  assert(j <= endNodes(i).getExtNormalsNum());
  
  return(endNodes(i).getExtNormals(j));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
bool
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
isGidEndNode(const UInt & gid) const
{
  return(endNodes.isG(gid));
}


//_________________________________________________________________________________________________
// GET FUNCTIONS - BIFURACTION NODES
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumBifurcationNodes() const
{
  assert(branchOk);
  return(bifurcationNodes.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
typename branch1d<GEOSHAPE,ELMAP,NODEMAP>::BIFURCATIONNODES &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getBifurcationNodes()
{
  assert(branchOk);
  return(bifurcationNodes);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename branch1d<GEOSHAPE,ELMAP,NODEMAP>::BIFURCATIONNODES &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getBifurcationNodes() const
{
  assert(branchOk);
  return(bifurcationNodes);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getBifurcationNodeL(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= bifurcationNodes.size());
  assert(MESH1D::nodes.isG(bifurcationNodes(i).getNodeGid()));
  
  return(MESH1D::nodes.getMapG(bifurcationNodes(i).getNodeGid()).getLid());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getBifurcationNodeG(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= bifurcationNodes.size());
  
  return(bifurcationNodes(i).getNodeGid());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getBifurcationIntToNumElements(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= bifurcationNodes.size());
    
  return(bifurcationNodes(i).getIntElementsGidNum());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getBifurcationToIntElementsL(const UInt & i, const UInt & j) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= bifurcationNodes.size());
  assert(j <= bifurcationNodes(i).getIntElementsGidNum());
  assert(MESH1D::elements.isG(bifurcationNodes(i).getIntElementsGid()(j)));
  
  return(MESH1D::elements.getMapG(bifurcationNodes(i).getIntElementsGid()(j)).getLid());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getBifurcationToIntElementsG(const UInt & i, const UInt & j) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= bifurcationNodes.size());
  assert(j <= bifurcationNodes(i).getIntElementsGidNum());
  
  return(bifurcationNodes(i).getIntElementsGid()(j));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getBifurcationToNumIntNormals(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= bifurcationNodes.size());
  
  return(bifurcationNodes(i).getIntNormalsNum());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d  &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getBifurcationToIntNormals(const UInt & i, const UInt & j) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= bifurcationNodes.size());
  
  return(bifurcationNodes(i).getIntNormals()(j));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getBifurcationToExtNodeL(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= bifurcationNodes.size());
  assert(MESH1D::nodes.isG(bifurcationNodes(i).getExtNodeGid()));
  
  return(MESH1D::nodes.getMapG(bifurcationNodes(i).getExtNodeGid()).getLid());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getBifurcationToExtNodeG(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= bifurcationNodes.size());
  
  return(bifurcationNodes(i).getExtNodeGid());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getBifurcationToExtElementL(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= bifurcationNodes.size());
  
  assert(MESH1D::elements.isG(bifurcationNodes(i).getExtElementGid()));
  
  return(MESH1D::elements.getMapG(bifurcationNodes(i).getExtElementGid()).getLid());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getBifurcationToExtElementG(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= bifurcationNodes.size());
  
  return(bifurcationNodes(i).getExtElementGid());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
getBifurcationToExtNormal(const UInt & i) const
{
  assert(branchOk);
  assert(i >= 1);
  assert(i <= bifurcationNodes.size());
  
  return(bifurcationNodes(i).getExtNormal());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>    
bool
branch1d<GEOSHAPE,ELMAP,NODEMAP>::
isGidBifurcation(const UInt & gid) const
{
  assert(branchOk);
  return(bifurcationNodes.isG(gid));
}

#endif
