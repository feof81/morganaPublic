/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef CONNECT3D_HPP
#define CONNECT3D_HPP

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

#include "mesh2d.hpp"
#include "mesh3d.hpp"


using namespace std;


/*! Topological connection information of a 3d mesh */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class connect3d
{
    /*! @name Typedefs */ //@{
  public:
    typedef ELMAP    GRID_ELMAP;
    typedef NODEMAP  GRID_NODEMAP;
    typedef GEOSHAPE GRID_GEOSHAPE;
    
    typedef GEOSHAPE                            GEOSHAPE3D;
    typedef typename   GEOSHAPE::GEOBSHAPE      GEOSHAPE2D;
    typedef typename GEOSHAPE2D::GEOBSHAPE      GEOSHAPE1D;
    
    typedef geoElement<GEOSHAPE3D>              GEOELEMENT3D;
    typedef geoElement<GEOSHAPE2D>              GEOELEMENT2D;
    typedef geoElement<GEOSHAPE1D>              GEOELEMENT1D;
    
    typedef mesh3d<GEOSHAPE3D,ELMAP,NODEMAP>    MESH3D;
    typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>    MESH2D;
    
    typedef pGraph<pGraphItem,NODEMAP,NODEMAP>       VERTEX_TO_VERTEX;
    typedef pGraph<pGraphItem,NODEMAP,ELMAP>         VERTEX_TO_ELEMENT;
    typedef pGraph<pGraphItemOriented,NODEMAP,ELMAP> VERTEX_TO_EDGE;
    typedef pGraph<pGraphItemSubLoc,ELMAP,ELMAP>     FACE_TO_ELEMENT;
    typedef pGraph<pGraphItem,ELMAP,ELMAP>           ELEMENT_TO_EDGE;
    typedef pGraph<pGraphItem,ELMAP,ELMAP>           ELEMENT_TO_FACE;
    typedef pGraph<pGraphItem,ELMAP,ELMAP>           ELEMENT_TO_ELEMENT;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 3;
    //@}
    
    /*! @name Internal control flags */ //@{
  public:
    bool connectCreated;
    bool boundaryConnectCreated;
    bool commDevLoaded;
    bool grid3dLoaded;
    bool grid2dLoaded;
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<const communicator> commDev;
    Teuchos::RCP<MESH3D> grid3d;
    Teuchos::RCP<MESH2D> grid2d;
    //@}
    
    /*! @name Internal graphs */ //@{
  public:
    VERTEX_TO_VERTEX   vertexToVertex;
    VERTEX_TO_ELEMENT  vertexToElement;
    VERTEX_TO_EDGE     vertexToEdge;
    FACE_TO_ELEMENT    faceToElement;
    ELEMENT_TO_EDGE    elementToEdge;
    ELEMENT_TO_FACE    elementToFace;
    ELEMENT_TO_ELEMENT elementToElement;
    //@}
    
    /*! @name Internal vectors */ //@{
  public:
    pVect<UInt,NODEMAP> vertexIsBoundary;
    pVect<UInt,ELMAP>   elementIsBoundary;
    pVect<UInt,ELMAP>   faceIsBoundary;
    pVect<UInt,ELMAP>   edgeIsBoundary;
    
    pVect<UInt,NODEMAP> vertexBVertex;
    pVect<UInt,ELMAP>   faceBFace;
    pVect<UInt,ELMAP>   edgeBEdge;
    //@}
    
    /*! @name Reference shapes */ //@{
  public:
    geoMapInterface<GEOSHAPE3D> refShape3d;
    geoMapInterface<GEOSHAPE2D> refShape2d;
    geoMapInterface<GEOSHAPE1D> refShape1d;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    connect3d();
    connect3d(const Teuchos::RCP<const communicator> & CommDev);
    connect3d(const communicator & CommDev);
    connect3d(const connect3d & C);
    connect3d & operator=(const connect3d & C);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    void setMesh3d(const Teuchos::RCP<MESH3D> & Grid3d);
    void setMesh2d(const Teuchos::RCP<MESH2D> & Grid2d);
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
    UInt getNumFaceToElement(const UInt & i) const;
    UInt getNumElementToEdge(const UInt & i) const;
    UInt getNumElementToFace(const UInt & i) const;
    UInt getNumElementToElement(const UInt & i) const;
    //@}
    
    /*! @name Get boundary functions */ //@{
  public:
    bool getVertexIsBoundary(const UInt & i) const;
    bool getElementIsBoundary(const UInt & i) const;
    bool getFaceIsBoundary(const UInt & i) const;
    bool getEdgeIsBoundary(const UInt & i) const;
    const UInt & getVertexBVertex(const UInt & i) const;
    const UInt & getFaceBFace(const UInt & i) const;
    const UInt & getEdgeBEdge(const UInt & i) const;
    //@}
    
    /*! @name Get volume functions */ //@{
  public:
    const UInt & getVertexToVertex(const UInt & i, const UInt & j) const;
    const UInt & getVertexToElement(const UInt & i, const UInt & j) const;
    const UInt & getVertexToEdge(const UInt & i, const UInt & j) const;
    const UInt & getFaceToElement(const UInt & i, const UInt & j) const;
    const UInt & getElementToEdge(const UInt & i, const UInt & j) const;
    const UInt & getElementToFace(const UInt & i, const UInt & j) const;
    const UInt & getElementToElement(const UInt & i, const UInt & j) const;
    //@}
    
    /*! @name Dump graphs functions */ //@{
  public:
    const VERTEX_TO_VERTEX   & getVertexToVertex() const;
    const VERTEX_TO_ELEMENT  & getVertexToElement() const;
    const VERTEX_TO_EDGE     & getVertexToEdge() const;
    const FACE_TO_ELEMENT    & getFaceToElement() const;
    const ELEMENT_TO_EDGE    & getElementToEdge() const;
    const ELEMENT_TO_FACE    & getElementToFace() const;
    const ELEMENT_TO_ELEMENT & getElementToElement() const;
    //@}
    
    /*! @name Dump vects functions */ //@{
  public:
    const pVect<UInt,NODEMAP>  & getVertexIsBoundary() const;
    const pVect<UInt,ELMAP>    & getElementIsBoundary() const;
    const pVect<UInt,ELMAP>    & getFaceIsBoundary() const;
    const pVect<UInt,ELMAP>    & getEdgeIsBoundary() const;
    const pVect<UInt,NODEMAP>  & getVertexBVertex() const;
    const pVect<UInt,ELMAP>    & getFaceBFace() const;
    const pVect<UInt,ELMAP>    & getEdgeBEdge() const;
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
connect3d()
{
  assert(GEOSHAPE::nDim == 3);
  
  connectCreated         = false;
  boundaryConnectCreated = false;
  commDevLoaded          = false;
  grid3dLoaded           = false;
  grid2dLoaded           = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
connect3d(const Teuchos::RCP<const communicator> & CommDev)
{
  assert(GEOSHAPE::nDim == 3);
  
  connectCreated         = false;
  boundaryConnectCreated = false;
  commDevLoaded          = true;
  grid3dLoaded           = false;
  grid2dLoaded           = false;
  
  commDev = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
connect3d(const communicator & CommDev)
{
  assert(GEOSHAPE::nDim == 3);
  
  connectCreated         = false;
  boundaryConnectCreated = false;
  commDevLoaded          = true;
  grid3dLoaded           = false;
  grid2dLoaded           = false;
  
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
connect3d(const connect3d & C)
{
  assert(GEOSHAPE::nDim == 3);
  
  connectCreated         = C.connectCreated;
  boundaryConnectCreated = C.boundaryConnectCreated;
  commDevLoaded          = C.commDevLoaded;
  grid3dLoaded           = C.grid3dLoaded;
  grid2dLoaded           = C.grid2dLoaded;
  
  commDev = C.commDev;
  grid3d  = C.grid3d;
  grid2d  = C.grid2d;
  
  vertexToVertex   = C.vertexToVertex;
  vertexToElement  = C.vertexToElement;
  vertexToEdge     = C.vertexToEdge;
  faceToElement    = C.faceToElement;
  elementToEdge    = C.elementToEdge;
  elementToFace    = C.elementToFace;
  elementToElement = C.elementToElement;

  vertexIsBoundary  = C.vertexIsBoundary;
  elementIsBoundary = C.elementIsBoundary;
  faceIsBoundary    = C.faceIsBoundary;
  edgeIsBoundary    = C.edgeIsBoundary;
    
  vertexBVertex = C.vertexBVertex;
  faceBFace     = C.faceBFace;
  edgeBEdge     = C.edgeBEdge;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect3d<GEOSHAPE,ELMAP,NODEMAP> &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
operator=(const connect3d & C)
{
  assert(GEOSHAPE::nDim == 3);
  
  connectCreated         = C.connectCreated;
  boundaryConnectCreated = C.boundaryConnectCreated;
  commDevLoaded          = C.commDevLoaded;
  grid3dLoaded           = C.grid3dLoaded;
  grid2dLoaded           = C.grid2dLoaded;
  
  commDev = C.commDev;
  grid3d  = C.grid3d;
  grid2d  = C.grid2d;
  
  vertexToVertex   = C.vertexToVertex;
  vertexToElement  = C.vertexToElement;
  vertexToEdge     = C.vertexToEdge;
  faceToElement    = C.faceToElement;
  elementToEdge    = C.elementToEdge;
  elementToFace    = C.elementToFace;
  elementToElement = C.elementToElement;

  vertexIsBoundary  = C.vertexIsBoundary;
  elementIsBoundary = C.elementIsBoundary;
  faceIsBoundary    = C.faceIsBoundary;
  edgeIsBoundary    = C.edgeIsBoundary;
    
  vertexBVertex = C.vertexBVertex;
  faceBFace     = C.faceBFace;
  edgeBEdge     = C.edgeBEdge;
  
  return *this;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
setMesh3d(const Teuchos::RCP<MESH3D> & Grid3d)
{
  grid3dLoaded = true;
  grid3d       = Grid3d;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
setMesh2d(const Teuchos::RCP<MESH2D> & Grid2d)
{
  grid2dLoaded = true;
  grid2d       = Grid2d;
}


//_________________________________________________________________________________________________
// COMPUTING FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
buildConnectivity()
{  
  //Control logic----------------------------------------------------
  assert(grid3dLoaded);
  connectCreated = true;
  
  //Clearing---------------------------------------------------------
  vertexToVertex.clear();
  vertexToElement.clear();
  vertexToEdge.clear();
  faceToElement.clear();
  elementToEdge.clear();
  elementToFace.clear();
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
  elementEdgeItem.resize(refShape3d.getNumEdges());
  
  //Allocations Faces------------------------------------------------
  typedef map<GEOELEMENT2D,ELMAP>      SETFACES;
  
  typedef typename SETFACES::iterator  TEMPFACESITERATOR;
  pair<GEOELEMENT2D,ELMAP>             valueFaces;
  SETFACES                             tempFaces; 
  pair<TEMPFACESITERATOR,bool>         paioFaces;
  pGraph<GEOELEMENT2D,ELMAP,NODEMAP>   tempGraphFaces;
  GEOELEMENT2D                         face(true);
  ELMAP                                faceMap;
  
  UInt numFaces = 1;
  ELMAP      elementFaceMap;
  pGraphItem elementFaceItem;
  elementFaceItem.resize(refShape3d.getNumFaces());
  
  //Other connecting allocations-------------------------------------
  ELMAP              faceToElementMap;
  pGraphItemSubLoc   faceToElementItem;
  pGraphItem         vertexToVertexItem;
  pGraphItem         elementToElementItem;
  pGraphItem         vertexToElementItem;
  pGraphItemOriented vertexToEdgeItem;
  
  //Other allocations------------------------------------------------ 
  UInt lId, gId;
  UInt id1, id2;
  UInt pid = commDev->rank();  

  //Edges build + elementToedge--------------------------------------
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    //Loop onto local edges
    for(UInt j=1; j <= refShape3d.getNumEdges(); ++j)
    {
      //Local edge preparation
      for(UInt k=1; k <= refShape1d.getNumPoints(); ++k)
      { 
	lId = refShape3d.edgeToPoint(j,k);
	gId = grid3d->getElementL(i).getCid(lId);
	spigolo.getCid(k) = gId;
      }
      
      spigolo.updateSorting();
      spigolo.setGeoId(0);
      
      spigoloMap = grid3d->getElements().getRowMapL(i);
      spigoloMap.setPid(pid);
      spigoloMap.setLid(numEdges);
      
      valueSpigoli.first  = spigolo;
      valueSpigoli.second = spigoloMap;
      
      //Insert edge
      paioSpigoli = tempSpigoli.insert(valueSpigoli);
      
      if(!paioSpigoli.second)
      {
	//already existing
	spigoloMap = (paioSpigoli.first)->second; //Map extraction
	elementEdgeItem.getCid(j) = spigoloMap.getLid();
      }
      else
      {
	//should be created
	elementEdgeItem.getCid(j) = numEdges;
	tempGraphEdges.push_back(spigolo,spigoloMap);
	++numEdges;
      }
    } //End loop onto local edges
    
    elementEdgeItem.updateSorting();
    
    elementEdgeMap = grid3d->getElements().getRowMapL(i);
    elementEdgeMap.setPid(pid);
    
    elementToEdge.push_back(elementEdgeItem,elementEdgeMap);
  }
  
  
  //Global numbering of the edges
  tempGraphEdges.setColMap(grid3d->getNodes().getMapRef());
  tempGraphEdges.updateColFinder();
  
  connectOptimization<GEOELEMENT1D,ELMAP,NODEMAP> edgeOptimizer;
  sVect<bool> edgeIsLocal = edgeOptimizer.getIsLocal(tempGraphEdges);
    
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> edgeGlobalMapper(commDev);
  edgeGlobalMapper.buildGlobalNumbering(tempGraphEdges,edgeIsLocal);  
  
  tempGraphEdges.updateRowFinder();
  tempGraphEdges.updateColFinder();
  
  
  //Element to Edge filling
  elementToEdge.setRowMap(grid3d->getElements().getRowMap());
  elementToEdge.setColMap(tempGraphEdges.getRowMap());
  
  elementToEdge.updateRowFinder();
  elementToEdge.updateColFinder();
  
  
  //Edges loading----------------------------------------------------
  grid3d->setEdges(tempGraphEdges);
  
  
  //Faces build + elementToFace + faceToElement----------------------
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    //Loop onto local faces
    for(UInt j=1; j <= refShape3d.getNumFaces(); ++j)
    {
      //Local face preparation
      for(UInt k=1; k <= refShape2d.getNumPoints(); ++k)
      { 
	lId = refShape3d.faceToPoint(j,k);
	gId = grid3d->getElementL(i).getCid(lId);
	face.getCid(k) = gId;
      }
      
      face.updateSorting();
      face.setGeoId(0);
      
      faceMap = grid3d->getElements().getRowMapL(i);
      faceMap.setPid(pid);
      faceMap.setLid(numFaces);
      
      valueFaces.first  = face;
      valueFaces.second = faceMap;
      
      //Insert face
      paioFaces = tempFaces.insert(valueFaces);
      
      if(!paioFaces.second)
      {
	//already existing
	faceMap = (paioFaces.first)->second;
	
	//face to element
	faceToElementItem = faceToElement.getItemL(faceMap.getLid());
	faceToElementItem.resize(2);
	faceToElementItem.setCid(2,i);
	faceToElementItem.setSubLocId(2,j);
	faceToElement.getItemL(faceMap.getLid()) = faceToElementItem;
	
	//elements to face
	elementFaceItem.getCid(j) = faceMap.getLid();
	
      }
      else
      {
	//should be created - face to element
	faceToElementMap.setLid(numFaces);
	faceToElementItem.resize(1);
	faceToElementItem.setCid(1,i);
	faceToElementItem.setSubLocId(1,j);
	faceToElement.push_back(faceToElementItem,faceToElementMap);
	
	//element to face
	elementFaceMap.setLid(i);
	elementFaceItem.getCid(j) = numFaces;
	
	//face creation
	tempGraphFaces.push_back(face,faceMap); 
	++numFaces;
      }
    } //End loop onto local faces
    
    elementFaceItem.updateSorting();
    elementToFace.push_back(elementFaceItem,elementFaceMap);
  }
  
  //Global numbering of the faces
  tempGraphFaces.setColMap(grid3d->getNodes().getMapRef());
  tempGraphFaces.updateColFinder();
  
  connectOptimization<GEOELEMENT2D,ELMAP,NODEMAP> faceOptimizer;
  sVect<bool> faceIsLocal = faceOptimizer.getIsLocal(tempGraphFaces);
   
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> faceGlobalMapper(commDev);
  faceGlobalMapper.buildGlobalNumbering(tempGraphFaces,faceIsLocal);
  
  tempGraphFaces.updateRowFinder();
  tempGraphFaces.updateColFinder();
  
  
  //Element to Face filling
  elementToFace.setRowMap(grid3d->getElements().getRowMap());
  elementToFace.setColMap(tempGraphFaces.getRowMap());
  
  elementToFace.updateRowFinder();
  elementToFace.updateColFinder();
  
  //Face to Element filling
  faceToElement.setRowMap(tempGraphFaces.getRowMap());
  faceToElement.setColMap(grid3d->getElements().getRowMap());
  
  faceToElement.updateRowFinder();
  faceToElement.updateColFinder();
  
  
  //Faces loading----------------------------------------------------
  grid3d->setFaces(tempGraphFaces);
  
  
  //Vertex to Vertex connectivity------------------------------------
  vertexToVertex.resize(grid3d->getNumNodes());
  
  for(UInt i=1; i <= grid3d->getNumEdges(); ++i)
  {    
    id1  = grid3d->getEdgeL(i).getCid(1);
    id2  = grid3d->getEdgeL(i).getCid(2);
    
    vertexToVertexItem = vertexToVertex.getItemL(id1);
    vertexToVertexItem.push_back(id2);
    vertexToVertex.getItemL(id1) = vertexToVertexItem;
    
    vertexToVertexItem = vertexToVertex.getItemL(id2);
    vertexToVertexItem.push_back(id1);
    vertexToVertex.getItemL(id2) = vertexToVertexItem;
  }
  
  vertexToVertex.setRowMap(grid3d->getNodes().getMapRef());
  vertexToVertex.setColMap(grid3d->getNodes().getMapRef());
  
  vertexToVertex.updateRowFinder();
  vertexToVertex.updateColFinder();
  
  
  //Element to element connectivity----------------------------------
  elementToElement.resize(grid3d->getNumElements());
  
  for(UInt i=1; i <= grid3d->getNumFaces(); ++i)
  {
    faceToElementItem = faceToElement.getItemL(i);
    
    if(faceToElementItem.size() == 2)
    {
      id1 = faceToElementItem.getCid(1);
      id2 = faceToElementItem.getCid(2);
      
      
      elementToElementItem = elementToElement.getItemL(id1);
      elementToElementItem.push_back(id2);
      elementToElement.getItemL(id1) = elementToElementItem;
      
      elementToElementItem = elementToElement.getItemL(id2);
      elementToElementItem.push_back(id1);
      elementToElement.getItemL(id2) = elementToElementItem;
    }
  }
  
  elementToElement.setRowMap(grid3d->getElements().getMapRef());
  elementToElement.setColMap(grid3d->getElements().getMapRef());
  
  elementToElement.updateRowFinder();
  elementToElement.updateColFinder();
  
  
  //Vertex to element connectivity-----------------------------------
  vertexToElement.resize(grid3d->getNumNodes());
  
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    for(UInt j=1; j <= refShape3d.getNumVertices(); j++)
    {
      id1 = grid3d->getElementL(i).getCid(j);
      
      vertexToElementItem = vertexToElement.getItemL(id1);
      vertexToElementItem.push_back(i);
      vertexToElement.getItemL(id1) = vertexToElementItem;
    }
  }
  
  vertexToElement.setRowMap(grid3d->getNodes().getMapRef());
  vertexToElement.setColMap(grid3d->getElements().getMapRef());
  
  vertexToElement.updateRowFinder();
  vertexToElement.updateColFinder();
  
  
  //Vertex edge connectivity-----------------------------------------
  vertexToEdge.resize(grid3d->getNumNodes());
  
  for(UInt i=1; i <= grid3d->getNumEdges(); ++i)
  {
    id1 = grid3d->getEdgeL(i).getCid(1);
    id2 = grid3d->getEdgeL(i).getCid(2);
    
    vertexToEdgeItem = vertexToEdge.getItemL(id1);
    vertexToEdgeItem.push_back(i,true);
    vertexToEdge.getItemL(id1) = vertexToEdgeItem;
    
    vertexToEdgeItem = vertexToEdge(id2);
    vertexToEdgeItem.push_back(i,false);
    vertexToEdge.getItemL(id2) = vertexToEdgeItem;
  }
  
  vertexToEdge.setRowMap(grid3d->getNodes().getMapRef());
  vertexToEdge.setColMap(grid3d->getEdges().getMapRef());
  
  vertexToEdge.updateRowFinder();
  vertexToEdge.updateColFinder();
  
  
  //Transferring map-------------------------------------------------
  grid3d->transferMap();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
buildBoundaryConnectivity()
{
  //Control logic----------------------------------------------------
  assert(grid3dLoaded);
  assert(grid2dLoaded);
  
  boundaryConnectCreated = true;
  
  
  //Allocations------------------------------------------------------
  GEOELEMENT3D element;
  UInt elemId;
  UInt cid;
  
  
  //Resizing---------------------------------------------------------
  vertexIsBoundary.resize(grid3d->getNumNodes());
  vertexBVertex.resize(grid3d->getNumNodes()); 
  for(UInt i=1; i <= vertexIsBoundary.size(); ++i)
  {
    vertexIsBoundary(i) = false;
  }
  
  elementIsBoundary.resize(grid3d->getNumElements());
  for(UInt i=1; i <= elementIsBoundary.size(); ++i)
  {
    elementIsBoundary(i) = false;
  }
  
  faceIsBoundary.resize(grid3d->getNumFaces());
  faceBFace.resize(grid3d->getNumFaces());
  for(UInt i=1; i <= faceIsBoundary.size(); ++i)
  {
    faceIsBoundary(i) = false;
  }
  
  edgeIsBoundary.resize(grid3d->getNumEdges());
  edgeBEdge.resize(grid3d->getNumEdges());
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

  for(UInt i=1; i <= grid2d->getNumNodes(); ++i) 
  {
    pointValue.first  = grid2d->getNodes().getDataL(i);
    pointValue.second = grid2d->getNodes().getMapL(i);
    
    pointSet.insert(pointValue);
  }  
  
  
  //Vertex to Bvertex association------------------------------------
  for(UInt i=1; i <= vertexIsBoundary.size(); ++i)
  {
    P = grid3d->getNodeL(i);
    
    if(pointSet.count(P) == 1)
    {
      pointIterator       = pointSet.find(P);
      vertexBVertex(i)    = pointIterator->second.getLid();
      vertexIsBoundary(i) = true;
    }
  }
  
  //Rowmap fixing
  vertexBVertex.setMap(grid3d->getNodes().getMapRef());
  vertexBVertex.updateFinder();
  
  vertexIsBoundary.setMap(grid3d->getNodes().getMapRef());
  vertexIsBoundary.updateFinder();
  
  
  //Boundary faces set construction----------------------------------
  typedef          map<GEOELEMENT2D,ELMAP>   BFACESET;      //Loading faces
  typedef typename BFACESET::iterator        BFACEITERATOR;
  typedef          pair<GEOELEMENT2D,ELMAP>  BFACEVALUE;
  
  BFACESET       bFaceSet;
  BFACEVALUE     bFaceValue;
  BFACEITERATOR  bFaceIterator;
  
  for(UInt i=1; i <= grid2d->getNumElements(); ++i) 
  {
    bFaceValue.first  = grid2d->getElements().getDataL(i);
    bFaceValue.second = grid2d->getElements().getMapL(i);
    
    bFaceSet.insert(bFaceValue);
  }
  
  
  //BFaces nodes remapping-------------------------------------------
  typedef typename MESH3D::GRAPH2D GRAPH2D;
  
  GEOELEMENT2D F(true);
  
  for(UInt i=1; i <= grid3d->getNumFaces(); ++i)
  {
    F = grid3d->getFaces().getItemL(i);
    
    //Change nodes numbering to match mesh3d
    for(UInt j=1; j <= F.size(); ++j)
    {
      cid = F.getCid(j);
      cid = vertexBVertex(cid);
      F.getCid(j) = cid;
    }
    
    F.updateSorting();

    assert(bFaceSet.count(F) <= 1);
    
    if(bFaceSet.count(F) == 1)
    {     
      bFaceIterator     = bFaceSet.find(F);
      faceBFace(i)      = bFaceIterator->second.getLid();
      faceIsBoundary(i) = true;
    }
  }
  
  faceBFace.setMap(grid3d->getFaces().getRowMap());
  faceBFace.updateFinder();
  
  faceIsBoundary.setMap(grid3d->getFaces().getRowMap());
  faceIsBoundary.updateFinder();
  
  
  //Boundary edges connectivity--------------------------------------
  typedef          map<GEOELEMENT1D,ELMAP>   BEDGESET;      //Loading faces
  typedef typename BEDGESET::iterator        BEDGEITERATOR;
  typedef          pair<GEOELEMENT1D,ELMAP>  BEDGEVALUE;
  
  BEDGESET       bEdgeSet;
  BEDGEVALUE     bEdgeValue;
  BEDGEITERATOR  bEdgeIterator;
  
  for(UInt i=1; i <= grid2d->getNumEdges(); ++i) 
  {
    bEdgeValue.first  = grid2d->getEdges().getDataL(i);
    bEdgeValue.second = grid2d->getEdges().getMapL(i);
    
    bEdgeSet.insert(bEdgeValue);
  }
  
  
  //BEdges edges remapping-------------------------------------------
  typedef typename MESH3D::GRAPH1D GRAPH1D;
  
  GEOELEMENT1D D(true);
  
  for(UInt i=1; i <= grid3d->getNumEdges(); ++i)
  {
    D = grid3d->getEdges().getItemL(i);
    
    //Change nodes numbering to match mesh3d
    for(UInt j=1; j <= D.size(); ++j)
    {
      cid = D.getCid(j);
      cid = vertexBVertex(cid);
      D.getCid(j) = cid;
    }
    
    D.updateSorting();

    assert(bEdgeSet.count(D) <= 1);
    
    if(bEdgeSet.count(D) == 1)
    {     
      bEdgeIterator     = bEdgeSet.find(D);
      edgeBEdge(i)      = bEdgeIterator->second.getLid();
      edgeIsBoundary(i) = true;
    }
  }
  
  edgeBEdge.setMap(grid3d->getEdges().getRowMap());
  edgeBEdge.updateFinder();
  
  edgeIsBoundary.setMap(grid3d->getEdges().getRowMap());
  edgeIsBoundary.updateFinder();
  
  
  //Element is boundary----------------------------------------------
  for(UInt i=1; i <= faceIsBoundary.size(); ++i)
  {
    if(faceIsBoundary(i))
    {
      assert((faceToElement(i).size() == 1) || (faceToElement(i).size() == 2));
      elemId = faceToElement(i).getCid(1);
      elementIsBoundary(elemId) = true;
    }
  }
  
  elementIsBoundary.setMap(grid3d->getElements().getRowMap());
  elementIsBoundary.updateFinder();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
clearConnectivity()
{
  connectCreated = false;
  
  vertexToVertex.clear();
  vertexToElement.clear();
  vertexToEdge.clear();
  faceToElement.clear();
  elementToEdge.clear();
  elementToFace.clear();
  elementToElement.clear();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
clearBoundaryConnectivity()
{
  boundaryConnectCreated = false;
  
  vertexIsBoundary.clear();
  elementIsBoundary.clear();
  faceIsBoundary.clear();
  edgeIsBoundary.clear();
  
  vertexBVertex.clear();
  faceBFace.clear();
  edgeBEdge.clear();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
clear()
{
  connectCreated         = false;
  boundaryConnectCreated = false;
  
  vertexToVertex.clear();
  vertexToElement.clear();
  vertexToEdge.clear();
  faceToElement.clear();
  elementToEdge.clear();
  elementToFace.clear();
  elementToElement.clear();
  
  vertexIsBoundary.clear();
  elementIsBoundary.clear();
  faceIsBoundary.clear();
  edgeIsBoundary.clear();
  
  vertexBVertex.clear();
  faceBFace.clear();
  edgeBEdge.clear();
}



//_________________________________________________________________________________________________
// GET NUM FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getNumVertexToVertex(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= vertexToVertex.rowSize());
  return(vertexToVertex.getItemL(i).size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getNumVertexToElement(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= vertexToElement.rowSize());
  return(vertexToElement.getItemL(i).size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getNumVertexToEdge(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= vertexToEdge.rowSize());
  return(vertexToEdge.getItemL(i).size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getNumFaceToElement(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= faceToElement.rowSize());
  return(faceToElement.getItemL(i).size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getNumElementToEdge(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= elementToEdge.rowSize());
  return(elementToEdge.getItemL(i).size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getNumElementToFace(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= elementToFace.rowSize());
  return(elementToFace.getItemL(i).size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
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
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexIsBoundary(const UInt & i) const
{
  assert(boundaryConnectCreated);
  assert(i >= 1);
  assert(i <= vertexIsBoundary.size());
  
  return(vertexIsBoundary(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
bool
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getElementIsBoundary(const UInt & i) const
{
  assert(boundaryConnectCreated);
  assert(i >= 1);
  assert(i <= elementIsBoundary.size());
  
  return(elementIsBoundary(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
bool
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getFaceIsBoundary(const UInt & i) const
{
  assert(boundaryConnectCreated);
  assert(i >= 1);
  assert(i <= faceIsBoundary.size());
  
  return(faceIsBoundary(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
bool
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeIsBoundary(const UInt & i) const
{
  assert(boundaryConnectCreated);
  assert(i >= 1);
  assert(i <= edgeIsBoundary.size());
  
  return(edgeIsBoundary(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexBVertex(const UInt & i) const
{
  assert(boundaryConnectCreated);
  assert(i >= 1);
  assert(i <= vertexBVertex.size());
  
  return(vertexBVertex(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getFaceBFace(const UInt & i) const
{
  assert(boundaryConnectCreated);
  assert(i >= 1);
  assert(i <= faceBFace.size());
  
  return(faceBFace(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
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
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
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
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
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
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexToEdge(const UInt & i, const UInt & j) const
{
  assert(connectCreated);
  assert(i >=1); assert(j >= 1);
  assert(i <= vertexToEdge.rowSizeL());
  assert(j <= vertexToEdge.getItemL(i).size());
  return(vertexToEdge.getCid_LL(i,j));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getFaceToElement(const UInt & i, const UInt & j) const
{
  assert(connectCreated);
  assert(i >=1); assert(j >= 1);
  assert(i <= faceToElement.rowSize());
  assert(j <= faceToElement.getItemL(i).size());
  return(faceToElement.getCid_LL(i,j));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
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
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getElementToFace(const UInt & i, const UInt & j) const
{
  assert(connectCreated);
  assert(i >=1); assert(j >= 1);
  assert(i <= elementToFace.rowSize());
  assert(j <= elementToFace.getItemL(i).size());
  return(elementToFace.getCid_LL(i,j));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
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
const typename connect3d<GEOSHAPE,ELMAP,NODEMAP>::VERTEX_TO_VERTEX &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexToVertex() const
{
  assert(connectCreated);
  return(vertexToVertex);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect3d<GEOSHAPE,ELMAP,NODEMAP>::VERTEX_TO_ELEMENT &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexToElement() const
{
  assert(connectCreated);
  return(vertexToElement);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect3d<GEOSHAPE,ELMAP,NODEMAP>::VERTEX_TO_EDGE &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexToEdge() const
{
  assert(connectCreated);
  return(vertexToEdge);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect3d<GEOSHAPE,ELMAP,NODEMAP>::FACE_TO_ELEMENT &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getFaceToElement() const
{
  assert(connectCreated);
  return(faceToElement);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect3d<GEOSHAPE,ELMAP,NODEMAP>::ELEMENT_TO_EDGE &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getElementToEdge() const
{
  assert(connectCreated);
  return(elementToEdge);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect3d<GEOSHAPE,ELMAP,NODEMAP>::ELEMENT_TO_FACE &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getElementToFace() const
{
  assert(connectCreated);
  return(elementToFace);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect3d<GEOSHAPE,ELMAP,NODEMAP>::ELEMENT_TO_ELEMENT &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
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
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexIsBoundary() const
{
  assert(boundaryConnectCreated);
  return(vertexIsBoundary);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const pVect<UInt,ELMAP> &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getElementIsBoundary() const
{
  assert(boundaryConnectCreated);
  return(elementIsBoundary);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const pVect<UInt,ELMAP> &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getFaceIsBoundary() const
{
  assert(boundaryConnectCreated);
  return(faceIsBoundary);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const pVect<UInt,ELMAP> &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeIsBoundary() const
{
  assert(boundaryConnectCreated);
  return(edgeIsBoundary);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const pVect<UInt,NODEMAP> &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexBVertex() const
{
  assert(boundaryConnectCreated);
  return(vertexBVertex);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const pVect<UInt,ELMAP> &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getFaceBFace() const
{
  assert(boundaryConnectCreated);
  return(faceBFace);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const pVect<UInt,ELMAP> &
connect3d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeBEdge() const
{
  assert(boundaryConnectCreated);
  return(edgeBEdge);
}


#endif
