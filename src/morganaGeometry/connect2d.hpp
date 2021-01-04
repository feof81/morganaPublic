/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef CONNECT2D_HPP
#define CONNECT2D_HPP

#include <map>
#include <set>

#include "pMapItem.h"
#include "pMapGlobalManip.h"

#include "pGraph.hpp"
#include "pGraphItem.h"
#include "pGraphItemOriented.h"
#include "pGraphItemSubLoc.h"
#include "pGraphGlobalManip.hpp"
#include "traitsConnectOptimization.hpp"

#include "mesh1d.hpp"
#include "mesh2d.hpp"

using namespace std;


/*! Topological connection information of a 2d mesh */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class connect2d
{
  /*! @name Typedefs */ //@{
  public:
    typedef ELMAP    GRID_ELMAP;
    typedef NODEMAP  GRID_NODEMAP;
    typedef GEOSHAPE GRID_GEOSHAPE;
    
    typedef GEOSHAPE                            GEOSHAPE2D;
    typedef typename   GEOSHAPE::GEOBSHAPE      GEOSHAPE1D;
    
    typedef geoElement<GEOSHAPE2D>              GEOELEMENT2D;
    typedef geoElement<GEOSHAPE1D>              GEOELEMENT1D;
    
    typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>    MESH2D;
    typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>    MESH1D;
    
    typedef pGraph<pGraphItem,NODEMAP,NODEMAP>       VERTEX_TO_VERTEX;
    typedef pGraph<pGraphItem,NODEMAP,ELMAP>         VERTEX_TO_ELEMENT;
    typedef pGraph<pGraphItemOriented,NODEMAP,ELMAP> VERTEX_TO_EDGE;
    typedef pGraph<pGraphItemSubLoc,ELMAP,ELMAP>     EDGE_TO_ELEMENT;
    typedef pGraph<pGraphItem,ELMAP,ELMAP>           ELEMENT_TO_EDGE;
    typedef pGraph<pGraphItem,ELMAP,ELMAP>           ELEMENT_TO_ELEMENT;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 2;
    //@}
    
    /*! @name Internal control flags */ //@{
  public:
    bool connectCreated;
    bool boundaryConnectCreated;
    bool commDevLoaded;
    bool grid2dLoaded;
    bool grid1dLoaded;
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<const communicator> commDev;
    Teuchos::RCP<MESH2D> grid2d;
    Teuchos::RCP<MESH1D> grid1d;
    //@}
    
    /*! @name Internal graphs */ //@{
  public:
    VERTEX_TO_VERTEX   vertexToVertex;
    VERTEX_TO_ELEMENT  vertexToElement;
    VERTEX_TO_EDGE     vertexToEdge;
    EDGE_TO_ELEMENT    edgeToElement;
    ELEMENT_TO_EDGE    elementToEdge;
    ELEMENT_TO_ELEMENT elementToElement;
    //@}
    
    /*! @name Internal vectors */ //@{
  public:
    pVect<UInt,NODEMAP> vertexIsBoundary;
    pVect<UInt,ELMAP>   elementIsBoundary;
    pVect<UInt,ELMAP>   edgeIsBoundary;
    
    pVect<UInt,NODEMAP> vertexBVertex;
    pVect<UInt,ELMAP>   edgeBEdge;
    //@}
    
    /*! @name Reference shapes */ //@{
  public:
    geoMapInterface<GEOSHAPE2D> refShape2d;
    geoMapInterface<GEOSHAPE1D> refShape1d;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    connect2d();
    connect2d(const Teuchos::RCP<const communicator> & CommDev);
    connect2d(const communicator & CommDev);
    connect2d(const connect2d & C);
    connect2d & operator=(const connect2d & C);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    void setMesh2d(const Teuchos::RCP<MESH2D> & Grid2d);
    void setMesh1d(const Teuchos::RCP<MESH1D> & Grid1d);
    //@}
    
    /*! @name Computing funcions */ //@{
  public:
    void buildConnectivity();
    void buildBoundaryConnectivity();
    void clearConnectivity();
    void clearBoundaryConnectivity();
    void clear();
    //@}
    
    /*! @name Get num functions */ //@{
  public:
    UInt getNumVertexToVertex(const UInt & i) const;
    UInt getNumVertexToElement(const UInt & i) const;
    UInt getNumVertexToEdge(const UInt & i) const;
    UInt getNumEdgeToElement(const UInt & i) const;    
    UInt getNumElementToEdge(const UInt & i) const;
    UInt getNumElementToElement(const UInt & i) const;
    //@}
    
    /*! @name Get boundary functions */ //@{
  public:
    bool getVertexIsBoundary(const UInt & i) const;
    bool getElementIsBoundary(const UInt & i) const;
    bool getEdgeIsBoundary(const UInt & i) const;
    const UInt & getVertexBVertex(const UInt & i) const;
    const UInt & getEdgeBEdge(const UInt & i) const;
    //@}
    
    /*! @name Get volume functions */ //@{
  public:
    const UInt & getVertexToVertex(const UInt & i, const UInt & j) const;
    const UInt & getVertexToElement(const UInt & i, const UInt & j) const;
    const UInt & getVertexToEdge(const UInt & i, const UInt & j) const;
    const UInt & getEdgeToElement(const UInt & i, const UInt & j) const;
    const UInt & getElementToEdge(const UInt & i, const UInt & j) const;
    const UInt & getElementToElement(const UInt & i, const UInt & j) const;
    //@}
    
    /*! @name Dump graphs functions */ //@{
  public:
    const VERTEX_TO_VERTEX   & getVertexToVertex() const;
    const VERTEX_TO_ELEMENT  & getVertexToElement() const;
    const VERTEX_TO_EDGE     & getVertexToEdge() const;
    const EDGE_TO_ELEMENT    & getEdgeToElement() const;
    const ELEMENT_TO_EDGE    & getElementToEdge() const;
    const ELEMENT_TO_ELEMENT & getElementToElement() const;
    //@}
    
    /*! @name Dump vects functions */ //@{
  public:
    const pVect<UInt,NODEMAP>  & getVertexIsBoundary() const;
    const pVect<UInt,ELMAP>    & getElementIsBoundary() const;
    const pVect<UInt,ELMAP>    & getEdgeIsBoundary() const;
    const pVect<UInt,NODEMAP>  & getVertexBVertex() const;
    const pVect<UInt,ELMAP>    & getEdgeBEdge() const;
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
connect2d()
{
  assert(GEOSHAPE::nDim == 2);
  
  connectCreated         = false;
  boundaryConnectCreated = false;
  commDevLoaded          = false;
  grid2dLoaded           = false;
  grid1dLoaded           = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
connect2d(const Teuchos::RCP<const communicator> & CommDev)
{
  assert(GEOSHAPE::nDim == 2);
  
  connectCreated         = false;
  boundaryConnectCreated = false;
  commDevLoaded          = true;
  grid2dLoaded           = false;
  grid1dLoaded           = false;
  
  commDev = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
connect2d(const communicator & CommDev)
{
  assert(GEOSHAPE::nDim == 2);
  
  connectCreated         = false;
  boundaryConnectCreated = false;
  commDevLoaded          = true;
  grid2dLoaded           = false;
  grid1dLoaded           = false;
  
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
connect2d(const connect2d & C)
{
  assert(GEOSHAPE::nDim == 2);
  
  connectCreated         = C.connectCreated;
  boundaryConnectCreated = C.boundaryConnectCreated;
  commDevLoaded          = C.commDevLoaded;
  grid2dLoaded           = C.grid2dLoaded;
  grid1dLoaded           = C.grid1dLoaded;
  
  commDev = C.commDev;
  grid2d  = C.grid2d;
  grid1d  = C.grid1d;
  
  vertexToVertex   = C.vertexToVertex;
  vertexToElement  = C.vertexToElement;
  vertexToEdge     = C.vertexToEdge;
  elementToEdge    = C.elementToEdge;
  elementToElement = C.elementToElement;

  vertexIsBoundary  = C.vertexIsBoundary;
  elementIsBoundary = C.elementIsBoundary;
  edgeIsBoundary    = C.edgeIsBoundary;
    
  vertexBVertex = C.vertexBVertex;
  edgeBEdge     = C.edgeBEdge;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect2d<GEOSHAPE,ELMAP,NODEMAP> &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
operator=(const connect2d & C)
{
  assert(GEOSHAPE::nDim == 2);
  
  connectCreated         = C.connectCreated;
  boundaryConnectCreated = C.boundaryConnectCreated;
  commDevLoaded          = C.commDevLoaded;
  grid2dLoaded           = C.grid2dLoaded;
  grid1dLoaded           = C.grid1dLoaded;
  
  commDev = C.commDev;
  grid2d  = C.grid2d;
  grid1d  = C.grid1d;
  
  vertexToVertex   = C.vertexToVertex;
  vertexToElement  = C.vertexToElement;
  vertexToEdge     = C.vertexToEdge;
  elementToEdge    = C.elementToEdge;
  elementToElement = C.elementToElement;

  vertexIsBoundary  = C.vertexIsBoundary;
  elementIsBoundary = C.elementIsBoundary;
  edgeIsBoundary    = C.edgeIsBoundary;
    
  vertexBVertex = C.vertexBVertex;
  edgeBEdge     = C.edgeBEdge;
  
  return *this;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
setMesh2d(const Teuchos::RCP<MESH2D> & Grid2d)
{
  grid2dLoaded = true;
  grid2d       = Grid2d;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
setMesh1d(const Teuchos::RCP<MESH1D> & Grid1d)
{
  grid1dLoaded = true;
  grid1d       = Grid1d;
}



//_________________________________________________________________________________________________
// COMPUTING FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
buildConnectivity()
{
  //Control logic----------------------------------------------------
  assert(grid2dLoaded);
  connectCreated = true;
  
  //Clearing---------------------------------------------------------
  vertexToVertex.clear();
  vertexToElement.clear();
  vertexToEdge.clear();
  edgeToElement.clear();
  elementToEdge.clear();
  elementToElement.clear();
  
  //Allocations Edges------------------------------------------------
  typedef map<GEOELEMENT1D,ELMAP>      SETEDGES;
  
  typedef typename SETEDGES::iterator  TEMPSPIGOLIITERATOR;
  pair<GEOELEMENT1D,ELMAP>             valueSpigoli;
  SETEDGES                             tempSpigoli;
  pair<TEMPSPIGOLIITERATOR,bool>       paioSpigoli;
  pGraph<GEOELEMENT1D,ELMAP,NODEMAP>   tempGraphEdges;
  GEOELEMENT1D                         spigolo(true);
  ELMAP                                spigoloMap;
  
  UInt numEdges = 1;
  ELMAP      elementEdgeMap;
  pGraphItem elementEdgeItem;
  elementEdgeItem.resize(refShape2d.getNumEdges());
  
  //Other connecting allocations-------------------------------------
  ELMAP              edgeToElementMap;
  pGraphItemSubLoc   edgeToElementItem;
  pGraphItem         vertexToVertexItem;
  pGraphItem         elementToElementItem;
  pGraphItem         vertexToElementItem;
  pGraphItemOriented vertexToEdgeItem;
  
  //Other allocations------------------------------------------------ 
  UInt id1, id2;
  UInt lid, gid;
  UInt pid = commDev->rank();
  
  //Edges build + elementToEdge + edgeToElement----------------------
  for(UInt i=1; i <= grid2d->getNumElements(); ++i)
  {
    //Loop onto local edges
    for(UInt j=1; j <= refShape2d.getNumEdges(); ++j)
    {
      //Local edge preparation
      for(UInt k=1; k <= refShape1d.getNumPoints(); ++k)
      { 
	lid = refShape2d.edgeToPoint(j,k);
	gid = grid2d->getElementL(i).getCid(lid);
	spigolo.setCid(k,gid);
      }
      
      spigolo.setGeoId(0);
      
      spigoloMap = grid2d->getElements().getRowMapL(i);
      spigoloMap.setPid(pid);
      spigoloMap.setLid(numEdges);
      
      valueSpigoli.first  = spigolo;
      valueSpigoli.second = spigoloMap;
      
      //Insert edge
      paioSpigoli = tempSpigoli.insert(valueSpigoli);
      
      
      if(!paioSpigoli.second)
      {
	//already existing
	spigoloMap = (paioSpigoli.first)->second;
	
	//edge to element
	edgeToElementItem = edgeToElement(spigoloMap.getLid());
	edgeToElementItem.resize(2);
	edgeToElementItem.setCid(2,i);
	edgeToElementItem.setSubLocId(2,j);
	edgeToElement.getItemL(spigoloMap.getLid()) = edgeToElementItem;
	
	//elements to face
	elementEdgeItem.setCid(j,spigoloMap.getLid());
      }
      else
      {
	//should be created - edge to element
	edgeToElementMap.setLid(numEdges);
	edgeToElementItem.resize(1);
	edgeToElementItem.setCid(1,i);
	edgeToElementItem.setSubLocId(1,j);
	edgeToElement.push_back(edgeToElementItem,edgeToElementMap);
	
	//element to edge
	elementEdgeMap.setLid(i);
	elementEdgeItem.setCid(j,numEdges);
	
	//face creation
	tempGraphEdges.push_back(spigolo,spigoloMap); 
	++numEdges;
      }
    }
    
    elementToEdge.push_back(elementEdgeItem,elementEdgeMap);
  }
  
  //Global numbering of the edges
  tempGraphEdges.setColMap(grid2d->getNodes().getMapRef());
  tempGraphEdges.updateColFinder();
  
  connectOptimization<GEOELEMENT1D,ELMAP,NODEMAP> edgeOptimizer;
  sVect<bool> edgeIsLocal = edgeOptimizer.getIsLocal(tempGraphEdges);  
  
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> edgeGlobalMapper(commDev);
  edgeGlobalMapper.buildGlobalNumbering(tempGraphEdges,edgeIsLocal);
  
  tempGraphEdges.updateRowFinder();
  tempGraphEdges.updateColFinder();
  
  //Element to Edge filling
  elementToEdge.setRowMap(grid2d->getElements().getRowMap());
  elementToEdge.setColMap(tempGraphEdges.getRowMap());
  
  elementToEdge.updateRowFinder();
  elementToEdge.updateColFinder();
  
  //Face to Element filling
  edgeToElement.setRowMap(tempGraphEdges.getRowMap());
  edgeToElement.setColMap(grid2d->getElements().getRowMap());
  
  edgeToElement.updateRowFinder();
  edgeToElement.updateColFinder();
  
  
  //Edges loading----------------------------------------------------
  grid2d->setEdges(tempGraphEdges);
  
  
  //Vertex to Vertex connectivity------------------------------------
  vertexToVertex.resize(grid2d->getNumNodes());
  
  for(UInt i=1; i <= grid2d->getNumEdges(); ++i)
  {    
    id1  = grid2d->getEdgeL(i).getCid(1);
    id2  = grid2d->getEdgeL(i).getCid(2);
    
    vertexToVertexItem = vertexToVertex.getItemL(id1);
    vertexToVertexItem.push_back(id2);
    vertexToVertex.getItemL(id1) = vertexToVertexItem;
    
    vertexToVertexItem = vertexToVertex.getItemL(id2);
    vertexToVertexItem.push_back(id1);
    vertexToVertex.getItemL(id2) = vertexToVertexItem;
  }
  
  vertexToVertex.setRowMap(grid2d->getNodes().getMapRef());
  vertexToVertex.setColMap(grid2d->getNodes().getMapRef());
  
  vertexToVertex.updateRowFinder();
  vertexToVertex.updateColFinder();
  
  
  //Element to element connectivity----------------------------------
  elementToElement.resize(grid2d->getNumElements());
  
  for(UInt i=1; i <= grid2d->getNumEdges(); ++i)
  {
    edgeToElementItem = edgeToElement.getItemL(i);
    
    if(edgeToElementItem.size() == 2)
    {
      id1 = edgeToElementItem.getCid(1);
      id2 = edgeToElementItem.getCid(2);
      
      
      elementToElementItem = elementToElement.getItemL(id1);
      elementToElementItem.push_back(id2);
      elementToElement.getItemL(id1) = elementToElementItem;
      
      elementToElementItem = elementToElement.getItemL(id2);
      elementToElementItem.push_back(id1);
      elementToElement.getItemL(id2) = elementToElementItem;
    }
  }
  
  elementToElement.setRowMap(grid2d->getElements().getMapRef());
  elementToElement.setColMap(grid2d->getElements().getMapRef());
  
  elementToElement.updateRowFinder();
  elementToElement.updateColFinder();
  
  
  //Vertex to element connectivity-----------------------------------
  vertexToElement.resize(grid2d->getNumNodes());
  
  for(UInt i=1; i <= grid2d->getNumElements(); ++i)
  {
    for(UInt j=1; j <= refShape2d.getNumVertices(); j++)
    {
      id1 = grid2d->getElementL(i).getCid(j);
      
      vertexToElementItem = vertexToElement.getItemL(id1);
      vertexToElementItem.push_back(i);
      vertexToElement.getItemL(id1) = vertexToElementItem;
    }
  }
  
  vertexToElement.setRowMap(grid2d->getNodes().getMapRef());
  vertexToElement.setColMap(grid2d->getElements().getMapRef());
  
  vertexToElement.updateRowFinder();
  vertexToElement.updateColFinder();
  
  
  //Vertex edge connectivity-----------------------------------------
  vertexToEdge.resize(grid2d->getNumNodes());
  
  for(UInt i=1; i <= grid2d->getNumEdges(); ++i)
  {
    id1 = grid2d->getEdgeL(i).getCid(1);
    id2 = grid2d->getEdgeL(i).getCid(2);
    
    vertexToEdgeItem = vertexToEdge.getItemL(id1);
    vertexToEdgeItem.push_back(i,true);
    vertexToEdge.getItemL(id1) = vertexToEdgeItem;
    
    vertexToEdgeItem = vertexToEdge(id2);
    vertexToEdgeItem.push_back(i,false);
    vertexToEdge.getItemL(id2) = vertexToEdgeItem;
  }
  
  vertexToEdge.setRowMap(grid2d->getNodes().getMapRef());
  vertexToEdge.setColMap(grid2d->getEdges().getMapRef());
  
  vertexToEdge.updateRowFinder();
  vertexToEdge.updateColFinder();
  
  
  //Transferring map-------------------------------------------------
  grid2d->transferMap();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
buildBoundaryConnectivity()
{
  //Control logic----------------------------------------------------
  assert(grid2dLoaded);
  assert(grid1dLoaded);
  
  boundaryConnectCreated = true;
  
  
  //Allocations------------------------------------------------------
  GEOELEMENT2D element;
  UInt elemId;
  UInt cid;
  
  
  //Resizing---------------------------------------------------------
  vertexIsBoundary.resize(grid2d->getNumNodes());
  vertexBVertex.resize(grid2d->getNumNodes()); 
  for(UInt i=1; i <= vertexIsBoundary.size(); ++i)
  {
    vertexIsBoundary(i) = false;
  }
  
  elementIsBoundary.resize(grid2d->getNumElements());
  for(UInt i=1; i <= elementIsBoundary.size(); ++i)
  {
    elementIsBoundary(i) = false;
  }
  
  edgeIsBoundary.resize(grid2d->getNumEdges());
  edgeBEdge.resize(grid2d->getNumEdges());
  for(UInt i=1; i <= edgeIsBoundary.size(); ++i)
  {
    edgeIsBoundary(i) = false;
  }
  
  
  //Boundary nodes set construction----------------------------------
  typedef          map<point3d,NODEMAP>     POINTSET;      //Loading nodes
  typedef typename POINTSET::iterator       POINTITERATOR;
  typedef          pair<point3d,NODEMAP>    POINTVALUE;
  
  POINTSET       pointSet;
  POINTVALUE     pointValue;
  POINTITERATOR  pointIterator;
  point3d        P;

  for(UInt i=1; i <= grid1d->getNumNodes(); ++i) 
  {
    pointValue.first  = grid1d->getNodes().getDataL(i);
    pointValue.second = grid1d->getNodes().getMapL(i);
    
    pointSet.insert(pointValue);
  }
  
  //Vertex to Bvertex association------------------------------------
  for(UInt i=1; i <= vertexIsBoundary.size(); ++i)
  {
    P = grid2d->getNodeL(i);
    
    if(pointSet.count(P) == 1)
    {
      pointIterator       = pointSet.find(P);
      vertexBVertex(i)    = pointIterator->second.getLid();
      vertexIsBoundary(i) = true;
    }
  }
  
  //Rowmap fixing
  vertexBVertex.setMap(grid2d->getNodes().getMapRef());
  vertexBVertex.updateFinder();
  
  vertexIsBoundary.setMap(grid2d->getNodes().getMapRef());
  vertexIsBoundary.updateFinder();
  
  
  //Boundary edges connectivity--------------------------------------
  typedef          map<GEOELEMENT1D,ELMAP>   BEDGESET;      //Loading faces
  typedef typename BEDGESET::iterator        BEDGEITERATOR;
  typedef          pair<GEOELEMENT1D,ELMAP>  BEDGEVALUE;
  
  BEDGESET       bEdgeSet;
  BEDGEVALUE     bEdgeValue;
  BEDGEITERATOR  bEdgeIterator;
  
  for(UInt i=1; i <= grid1d->getNumElements(); ++i) 
  {
    bEdgeValue.first  = grid1d->getElements().getDataL(i);
    bEdgeValue.second = grid1d->getElements().getMapL(i);
    
    bEdgeSet.insert(bEdgeValue);
  }
  
  
  //BEdges edges remapping-------------------------------------------
  typedef typename MESH1D::GRAPH1D GRAPH1D;
  
  GEOELEMENT1D D(true);
  
  for(UInt i=1; i <= grid2d->getNumEdges(); ++i)
  {
    D = grid2d->getEdges().getItemL(i);
    
    //Change nodes numbering to match mesh3d
    for(UInt j=1; j <= D.size(); ++j)
    {
      cid = D.getCid(j);
      cid = vertexBVertex(cid);
      D.setCid(j,cid);
    }
    
    assert(bEdgeSet.count(D) <= 1);
    
    if(bEdgeSet.count(D) == 1)
    {     
      bEdgeIterator     = bEdgeSet.find(D);
      edgeBEdge(i)      = bEdgeIterator->second.getLid();
      edgeIsBoundary(i) = true;
    }
  }
  
  edgeBEdge.setMap(grid2d->getEdges().getRowMap());
  edgeBEdge.updateFinder();
  
  edgeIsBoundary.setMap(grid2d->getEdges().getRowMap());
  edgeIsBoundary.updateFinder();
  
  
  //Element is boundary----------------------------------------------
  for(UInt i=1; i <= edgeIsBoundary.size(); ++i)
  {
    if(edgeIsBoundary(i))
    {
      assert((edgeToElement(i).size() == 1) || (edgeToElement(i).size() == 2));
      elemId = edgeToElement(i).getCid(1);
      elementIsBoundary(elemId) = true;
    }
  }
  
  elementIsBoundary.setMap(grid2d->getElements().getRowMap());
  elementIsBoundary.updateFinder();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
clearConnectivity()
{
  connectCreated = false;
  
  vertexToVertex.clear();
  vertexToElement.clear();
  vertexToEdge.clear();
  edgeToElement.clear();
  elementToEdge.clear();
  elementToElement.clear();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
clearBoundaryConnectivity()
{
  boundaryConnectCreated = false;
  
  vertexIsBoundary.clear();
  elementIsBoundary.clear();
  edgeIsBoundary.clear();
  edgeIsBoundary.clear();
  
  vertexBVertex.clear();
  edgeBEdge.clear();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
clear()
{
  connectCreated         = false;
  boundaryConnectCreated = false;
  
  vertexToVertex.clear();
  vertexToElement.clear();
  vertexToEdge.clear();
  edgeToElement.clear();
  elementToEdge.clear();
  elementToElement.clear();
  
  vertexIsBoundary.clear();
  elementIsBoundary.clear();
  edgeIsBoundary.clear();
  
  vertexBVertex.clear();
  edgeBEdge.clear();
}



//_________________________________________________________________________________________________
// GET NUM FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumVertexToVertex(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= vertexToVertex.rowSize());
  return(vertexToVertex.getItemL(i).size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumVertexToElement(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= vertexToElement.rowSize());
  return(vertexToElement.getItemL(i).size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumVertexToEdge(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= vertexToEdge.rowSize());
  return(vertexToEdge.getItemL(i).size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumEdgeToElement(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= edgeToElement.rowSize());
  return(edgeToElement.getItemL(i).size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumElementToEdge(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= elementToEdge.rowSize());
  return(elementToEdge.getItemL(i).size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumElementToElement(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= elementToElement.rowSize());
  return(elementToElement.getItemL(i).size());
}



//_________________________________________________________________________________________________
// GET BOUNDARY FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
bool
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexIsBoundary(const UInt & i) const
{
  assert(boundaryConnectCreated);
  assert(i >= 1);
  assert(i <= vertexIsBoundary.size());
  
  return(vertexIsBoundary(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
bool
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getElementIsBoundary(const UInt & i) const
{
  assert(boundaryConnectCreated);
  assert(i >= 1);
  assert(i <= elementIsBoundary.size());
  
  return(elementIsBoundary(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
bool
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeIsBoundary(const UInt & i) const
{
  assert(boundaryConnectCreated);
  assert(i >= 1);
  assert(i <= edgeIsBoundary.size());
  
  return(edgeIsBoundary(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexBVertex(const UInt & i) const
{
  assert(boundaryConnectCreated);
  assert(i >= 1);
  assert(i <= vertexBVertex.size());
  
  return(vertexBVertex(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeBEdge(const UInt & i) const
{
  assert(boundaryConnectCreated);
  assert(i >= 1);
  assert(i <= edgeBEdge.size());
  
  return(edgeBEdge(i));
}



//_________________________________________________________________________________________________
// GET VOLUME FUINCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexToVertex(const UInt & i, const UInt & j) const
{
  assert(connectCreated);
  assert(i >=1); assert(j >= 1);
  assert(i <= vertexToVertex.rowSizeL());
  assert(j <= vertexToVertex.getItemL(i).size());
  return(vertexToVertex.getCid_LL(i,j));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexToElement(const UInt & i, const UInt & j) const
{
  assert(connectCreated);
  assert(i >=1); assert(j >= 1);
  assert(i <= vertexToElement.rowSize());
  assert(j <= vertexToElement.getItemL(i).size());
  return(vertexToElement.getCid_LL(i,j));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexToEdge(const UInt & i, const UInt & j) const
{
  assert(connectCreated);
  assert(i >=1); assert(j >= 1);
  assert(i <= vertexToEdge.rowSize());
  assert(j <= vertexToEdge.getItemL(i).size());
  return(vertexToEdge.getCid_LL(i,j));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeToElement(const UInt & i, const UInt & j) const
{
  assert(connectCreated);
  assert(i >=1); assert(j >= 1);
  assert(i <= edgeToElement.rowSize());
  assert(j <= edgeToElement.getItemL(i).size());
  return(edgeToElement.getCid_LL(i,j));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getElementToEdge(const UInt & i, const UInt & j) const
{
  assert(connectCreated);
  assert(i >=1); assert(j >= 1);
  assert(i <= elementToEdge.rowSize());
  assert(j <= elementToEdge.getItemL(i).size());
  return(elementToEdge.getCid_LL(i,j));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getElementToElement(const UInt & i, const UInt & j) const
{
  assert(connectCreated);
  assert(i >=1); assert(j >= 1);
  assert(i <= elementToElement.rowSize());
  assert(j <= elementToElement.getItemL(i).size());
  return(elementToElement.getCid_LL(i,j));
}



//_________________________________________________________________________________________________
// DUMP GRAPH FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect2d<GEOSHAPE,ELMAP,NODEMAP>::VERTEX_TO_VERTEX &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexToVertex() const
{
  assert(connectCreated);
  return(vertexToVertex);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect2d<GEOSHAPE,ELMAP,NODEMAP>::VERTEX_TO_ELEMENT &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexToElement() const
{
  assert(connectCreated);
  return(vertexToElement);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect2d<GEOSHAPE,ELMAP,NODEMAP>::VERTEX_TO_EDGE &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexToEdge() const
{
  assert(connectCreated);
  return(vertexToEdge);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect2d<GEOSHAPE,ELMAP,NODEMAP>::EDGE_TO_ELEMENT &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeToElement() const
{
  assert(connectCreated);
  return(edgeToElement);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect2d<GEOSHAPE,ELMAP,NODEMAP>::ELEMENT_TO_EDGE &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getElementToEdge() const
{
  assert(connectCreated);
  return(elementToEdge);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect2d<GEOSHAPE,ELMAP,NODEMAP>::ELEMENT_TO_ELEMENT &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getElementToElement() const
{
  assert(connectCreated);
  return(elementToElement);
}



//_________________________________________________________________________________________________
// DUNMP VECT FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const pVect<UInt,NODEMAP> &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexIsBoundary() const
{
  assert(boundaryConnectCreated);
  return(vertexIsBoundary);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const pVect<UInt,ELMAP> &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getElementIsBoundary() const
{
  assert(boundaryConnectCreated);
  return(elementIsBoundary);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const pVect<UInt,ELMAP> &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeIsBoundary() const
{
  assert(boundaryConnectCreated);
  return(edgeIsBoundary);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const pVect<UInt,NODEMAP> &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexBVertex() const
{
  assert(boundaryConnectCreated);
  return(vertexBVertex);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const pVect<UInt,ELMAP> &
connect2d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeBEdge() const
{
  assert(boundaryConnectCreated);
  return(edgeBEdge);
}

#endif
