/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef DUAL2D_HPP
#define DUAL2D_HPP

#include "algorithm"
#include "traitsGeometry.hpp"

#include "mesh2d.hpp"
#include "connect2d.hpp"
#include "geoMapSupport1d.hpp"
#include "geoMapSupport2d.hpp"
#include <boost/concept_check.hpp>

using namespace std;



//_______________________________________________________________________________________________________
// DUMMY IMPLEMENTATION
//-------------------------------------------------------------------------------------------------------

/*! Dual mesh 2d topology and geometry - dummy implementation */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class dual2d
{ };



//_______________________________________________________________________________________________________
// LINEAR TRIANGLE SPECILIZATION
//-------------------------------------------------------------------------------------------------------
/*! Dual mesh 2d topology and geometry - \c linearTriangle implementation. Supports class for the implementation of vertex-centered finite volume schemes. Every triangle is decomposed in sub-cells and then aggregated to the mesh vertices.
Also the interfaces between the cells are created.

<b>Data</b>
<ol>
<li> dualPoints = the dual points: they are the union of the original grid vertices and of the barycentric nodes of the elements and of the edges;
<li> dualCells = the grid subcells. Every triangle is divide into 3 x 2 subcells. Every subcell is a right-hand triangle. Every subcell is linked to a vertex of the 
original triangle, therefore each of the original 3 vertices of the triangle has 2 subcells. Every subcell is built joining a vertex of the original 
triangle to and edge midpoint and the barycenter of the triangle. The order of the last to nodes is decided in order to create a right-hand triangle;
<li> interfaces = the interfaces between the subcells. The finite volume method require the definition of interaces between cells. In our case the cell is the
union of the subcells sharing a given vertex of the original mesh. In the triangle case the interfaces are made of lines. Denoting with p1, p2, the nodes
of the interface the node p1 is on ad edge and the node p2 is the barycenter of the element.
</ol>

<b>Connectivity data</b>
<ol>
<li> subCellToElement: the connectivity of the subCells to the father element;
<li> subCellToNode: the connetivity of the subCell to a node;
<li> interfaceToEdge: the connectivity of the interface to the edge;
<li> interfaceToElement: the connectivity of the interface to the element
</ol>
  
<b>Finite volume data</b>
<ol>
<li> cellVolume: the aggregated volume of the cell;
<li> subCellVolume: the volume of each sub cell
<li> subCellN: the orthogonal vector to the plane of the subcell
<li> edgeN: the mean normal vector to the interface. The normal is obtained as a weighted normal vector of the nornals of the triangles that compose the interface.
Every interface is associated to an edge.
<li> edgeSurf: the surface of the interface
</ol>
*/
template<typename ELMAP, typename NODEMAP>
class dual2d<linearTriangle,ELMAP,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef linearTriangle             GEOSHAPE2D;
    typedef linearTriangle::GEOBSHAPE  GEOSHAPE1D;
    typedef geoElement<GEOSHAPE2D>     GEOELEMENT2D;
    typedef geoElement<GEOSHAPE1D>     GEOELEMENT1D;
    
    typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>    MESH2D;
    typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP> CONNECT2D;
    
    typedef pVect<point3d,NODEMAP>              NODESVECT;
    typedef pGraph<GEOELEMENT2D,ELMAP,NODEMAP>  GRAPH2D;
    typedef pGraph<GEOELEMENT1D,ELMAP,NODEMAP>  GRAPH1D;
    
    typedef pVect<Real,ELMAP>       CELLVOLUMESVECT;
    typedef pVect<point3d,ELMAP>    EDGENSVECT;
    typedef pVect<point3d,ELMAP>    SUBCELLNSVECT;
    typedef pVect<Real,ELMAP>       EDGESSURFVECT;
    
    typedef pVect<UInt,ELMAP>       SUBCELL_TO_ELEMENT;
    typedef pVect<UInt,ELMAP>       SUBCELL_TO_NODE;
    typedef pVect<UInt,ELMAP>       INTERFACE_TO_EDGE;
    typedef pVect<UInt,ELMAP>       INTERFACE_TO_ELEMENT;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 2;
    //@}
      
    /*! @name Internal flags */ //@{ 
  public:
    bool commDevLoaded;
    bool grid2dLoaded;
    bool connectGrid2dLoaded;
    bool dualCreated;
    bool connectCreated;
    bool finiteVolumeCreated;
    //@}
      
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<communicator> commDev;
    Teuchos::RCP<MESH2D>       grid2d;
    Teuchos::RCP<CONNECT2D>    connectGrid2d;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    geoMapInterface<GEOSHAPE2D>  refShape2d;
    NODESVECT dualPoints;
    GRAPH2D   subCells;
    GRAPH1D   interfaces;
    //@}
    
    /*! @name Finite Volumes data */ //@{ 
  public:
    CELLVOLUMESVECT cellVolume;
    CELLVOLUMESVECT subCellVolume;
    SUBCELLNSVECT   subCellN;
    EDGENSVECT      edgeN;
    EDGESSURFVECT   edgeSurf;
    //@}
    
    /*! @name Internal connectivity data */ //@{ 
  public:
    SUBCELL_TO_ELEMENT   subCellToElement;
    SUBCELL_TO_NODE      subCellToNode;  
    INTERFACE_TO_EDGE    interfaceToEdge;
    INTERFACE_TO_ELEMENT interfaceToElement;
    //@}
    
    /*! @name Constructors */ //@{ 
  public:
    dual2d();
    dual2d(const Teuchos::RCP<communicator> & CommDev);
    dual2d(communicator & CommDev);
    void setMesh2d(const Teuchos::RCP<MESH2D> & Grid2d);
    void setMesh2d(const MESH2D & Grid2d);
    void setConnect2d(const Teuchos::RCP<CONNECT2D> & Connect2d);
    void setConnect2d(const CONNECT2D & Connect2d);
    //@}
    
    /*! @name Funzioni interne e buildup */ //@{ 
  public:    
    /*! The dual nodes is the union of the original nodes + the elements barycentric nodes + the edges ones.
    Given the id of a geoItem (such as an element, face, etc) this function produces the local id of the dual node associated.
    \param flag = 1 the original vertices, = 2 the element barycenter, = 3 the edge one
    \param lid the local numbering of the geoItem (such as an element, face, etc) */
    UInt getLocalNodeId(const UInt & flag, const UInt & lid) const;
      
    /*! Build dual mesh */
    void buildDualMesh();
      
    /*! Build finite volume data. Requires an explicit construction of the dual mesh */
    void buildFiniteVolumeData();
    
    /*! Build finite volume data without explicitly constructing the dual grid.
    The connections and the finite volume data are build. No global ownership is built*/
    //TODO : this function should be reimplemented and enhanced when the data needed by the 
    // finite volume solvers will become clear - this is due in 2013
    void buildFiniteVolumeDataOnly();
    
    /*! Clear the dual data only */
    void clearDual();
    
     /*! Clear cennection data only */
    void clearConnect();
      
    /*! Clear the finite volume data only */
    void clearFiniteVolume();
      
    /*! Clear all the data */
    void clear();
    //@}
};


template<typename ELMAP, typename NODEMAP>
dual2d<linearTriangle,ELMAP,NODEMAP>::
dual2d()
{
  commDevLoaded       = false;
  grid2dLoaded        = false;
  connectGrid2dLoaded = false;
  dualCreated         = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
}

template<typename ELMAP, typename NODEMAP>
dual2d<linearTriangle,ELMAP,NODEMAP>::
dual2d(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded       = true;
  grid2dLoaded        = false;
  connectGrid2dLoaded = false;
  dualCreated         = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
  
  commDev = CommDev;
}

template<typename ELMAP, typename NODEMAP>
dual2d<linearTriangle,ELMAP,NODEMAP>::
dual2d(communicator & CommDev)
{
  commDevLoaded       = true;
  grid2dLoaded        = false;
  connectGrid2dLoaded = false;
  dualCreated         = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
  
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename ELMAP, typename NODEMAP>
void
dual2d<linearTriangle,ELMAP,NODEMAP>::
setMesh2d(const Teuchos::RCP<MESH2D> & Grid2d)
{
  grid2dLoaded = true;
  grid2d       = Grid2d;
}

template<typename ELMAP, typename NODEMAP>
void
dual2d<linearTriangle,ELMAP,NODEMAP>::
setMesh2d(const MESH2D & Grid2d)
{
  grid2dLoaded = true;
  grid2d       = Teuchos::rcp(new CONNECT2D(Grid2d));
}

template<typename ELMAP, typename NODEMAP>
void
dual2d<linearTriangle,ELMAP,NODEMAP>::
setConnect2d(const Teuchos::RCP<CONNECT2D> & Connect2d)
{
  connectGrid2dLoaded = true;
  connectGrid2d       = Connect2d;
}

template<typename ELMAP, typename NODEMAP>
void
dual2d<linearTriangle,ELMAP,NODEMAP>::
setConnect2d(const CONNECT2D & Connect2d)
{
  connectGrid2dLoaded = true;
  connectGrid2d       = Teuchos::rcp(new CONNECT2D(Connect2d));
}

template<typename ELMAP, typename NODEMAP> 
UInt
dual2d<linearTriangle,ELMAP,NODEMAP>::
getLocalNodeId(const UInt & flag, const UInt & lid) const
{ 
  return( UInt(flag > 1) * grid2d->getNumNodes()
        + UInt(flag > 2) * grid2d->getNumElements()
        + lid );
}

template<typename ELMAP, typename NODEMAP> 
void
dual2d<linearTriangle,ELMAP,NODEMAP>::
buildDualMesh()
{
  //Assert_________________________________________________________________________________________
  assert(commDevLoaded);
  assert(grid2dLoaded);
  assert(connectGrid2dLoaded);
  
  dualCreated    = true;
  connectCreated = true;
  
  
  //Clearing_______________________________________________________________________________________
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
  
  interfaceToElement.clear();
  interfaceToEdge.clear();
  
  subCellToElement.clear();
  subCellToNode.clear();  
  
  
  
  //Num local and global elements__________________________________________________________________
  UInt numLocalElements = grid2d->getNumElements();
  UInt numLocalEdges    = grid2d->getNumEdges();
  
  pMapGlobalManip<NODEMAP> vrMapManip(commDev);
  pMapGlobalManip<ELMAP>   elMapManip(commDev);
  pMapGlobalManip<ELMAP>   edMapManip(commDev);
  
  pMap<NODEMAP> vrMap = grid2d->getNodes().getMapRef();
  pMap<ELMAP>   elMap = grid2d->getElements().getRowMap();
  pMap<ELMAP>   edMap = grid2d->getEdges().getRowMap();
  
  UInt numGlobalNodes    = vrMapManip.sizeG(vrMap);
  UInt numGlobalElements = elMapManip.sizeG(elMap);
  
  
  //Allocations____________________________________________________________________________________
  UInt pid = commDev->rank();
  
  UInt    id1, id2, id3;
  point3d P1, P2, P3, PM;
  
  GEOELEMENT2D subCell(true);
  GEOELEMENT1D interface(true);
  
  UInt numSubCells = 1;
  UInt numInterfaces = 1;
  
  UInt ba, locStray, elGid, geoId, edgeId, numPoints;
  sVect<UInt> p(3), e(3);
  
  ELMAP   elMapItem;
  NODEMAP nodeMapItem;
  
  
  //Nodes construction_____________________________________________________________________________
  dualPoints = grid2d->getNodes();   //Orginal nodes
  numPoints  = dualPoints.size() + 1;
  
  for(UInt i=1; i <= numLocalElements; ++i)    //Nodi baricentrici
  {
    id1 = grid2d->getElementL(i).getCid(1);
    id2 = grid2d->getElementL(i).getCid(2);
    id3 = grid2d->getElementL(i).getCid(3);
    
    P1  = grid2d->getNodeL(id1);
    P2  = grid2d->getNodeL(id2);
    P3  = grid2d->getNodeL(id3);
    
    PM  = (P1 + P2 + P3) / 3.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes +  grid2d->getElements().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  for(UInt i=1; i <= numLocalEdges; ++i)       //Nodi degli spigoli 
  {
    id1 = grid2d->getEdgeL(i).getCid(1);
    id2 = grid2d->getEdgeL(i).getCid(2);
    
    P1  = grid2d->getNodeL(id1);
    P2  = grid2d->getNodeL(id2);
    
    PM = (P1 + P2) / 2.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes + numGlobalElements + grid2d->getEdges().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  //Map ownership fixing
  pMap<NODEMAP> tempNodeMap = dualPoints.getMapRef();
  
  nodesMapFixer(tempNodeMap,*commDev);
  dualPoints.setMap(tempNodeMap);
  dualPoints.updateFinder();
  
  //Checking
  pVectGlobalManip<point3d,NODEMAP> dualNodesChecker(commDev);
  assert(dualNodesChecker.check(dualPoints));
  
  
  //Sub Cells construction_________________________________________________________________________
  for(UInt i=1; i <= numLocalElements; ++i)
  {
    elMapItem = grid2d->getElements().getRowMapL(i);
    elGid     = grid2d->getElements().getRowMapL(i).getGid();
    
    //Triangle nodes
    for(UInt j=1; j <= 3; ++j) 
    {
      p(j) = grid2d->getElementL(i).getCid(j);
    }
    
    //Hexa barycenter
    ba = getLocalNodeId(2,i);
    
    //Nodes on the Hexa edges
    for(UInt j=1; j <= 3; ++j)
    {
      e(j) = connectGrid2d->getElementToEdge(i,j);
      e(j) = getLocalNodeId(3,e(j));
    }
    
    //Subcells creation
    locStray = 1;
    elGid    = grid2d->getElements().getRowMapL(i).getGid();
    geoId    = grid2d->getElements().getItemL(i).getGeoId();
    
    elMapItem.setPid(pid);
    subCell.setGeoId(geoId);
    
    
    //Subcell 1--------------------------------------------
    subCell.setCid(1,p(1));
    subCell.setCid(2,e(1));
    subCell.setCid(3,ba);
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(6 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(1));
    
    
    //Subcell 2--------------------------------------------
    subCell.setCid(1,p(2));
    subCell.setCid(2,ba);
    subCell.setCid(3,e(1));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(6 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(2));
    
    
    //Subcell 3--------------------------------------------
    subCell.setCid(1,p(2));
    subCell.setCid(2,e(2));
    subCell.setCid(3,ba);
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(6 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(2));
    
    
    //Subcell 4--------------------------------------------
    subCell.setCid(1,p(3));
    subCell.setCid(2,ba);
    subCell.setCid(3,e(2));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(6 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(3));
    
    
    //Subcell 5--------------------------------------------
    subCell.setCid(1,p(3));
    subCell.setCid(2,e(3));
    subCell.setCid(3,ba);
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(6 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(3));
    
    
    //Subcell 6--------------------------------------------
    subCell.setCid(1,p(1));
    subCell.setCid(2,ba);
    subCell.setCid(3,e(3));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(6 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(1));
    
    
    //Interfaces creation__________________________________________________________________________
    locStray = 1;
    
    //Interface 1------------------------------------------
    edgeId = 1;
    
    interface.setCid(1,ba);
    interface.setCid(2,e(1));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(3 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid2d->getElementToEdge(i,edgeId));
    interfaceToElement.push_back(elMapItem, i);
    
    
    //Interface 2------------------------------------------
    edgeId = 2;
    
    interface.setCid(1,ba);
    interface.setCid(2,e(2));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(3 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid2d->getElementToEdge(i,edgeId));
    interfaceToElement.push_back(elMapItem, i);
    
    
    //Interface 3------------------------------------------
    edgeId = 3;
    
    interface.setCid(1,ba);
    interface.setCid(2,e(3));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(3 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid2d->getElementToEdge(i,edgeId));
    interfaceToElement.push_back(elMapItem, i);
  }
  
  
  //Map fixing for subCells
  subCells.setColMap(dualPoints.getMapRef());
  subCells.updateRowFinder();
  subCells.updateColFinder();
  
  subCellToElement.setMap(subCells.getRowMap());
  subCellToElement.updateFinder();

  subCellToNode.setMap(subCells.getRowMap());
  subCellToNode.updateFinder();
  
  //Checking for subCells
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> subCellsCheker(commDev);
  assert(subCellsCheker.check(subCells));
  
  
  //MapFixing for interfaces
  interfaces.setColMap(dualPoints.getMapRef());
  interfaces.updateRowFinder();
  interfaces.updateColFinder();
  
  interfaceToEdge.setMap(interfaces.getRowMap());
  interfaceToEdge.updateFinder();
  
  interfaceToElement.setMap(interfaces.getRowMap());
  interfaceToElement.updateFinder();
  
  //Checking interfaces map
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> interfacesCheker(commDev);
  assert(interfacesCheker.check(interfaces));
}

template<typename ELMAP, typename NODEMAP> 
void
dual2d<linearTriangle,ELMAP,NODEMAP>::
buildFiniteVolumeDataOnly()
{
  //Assert_________________________________________________________________________________________
  assert(commDevLoaded);
  assert(grid2dLoaded);
  assert(connectGrid2dLoaded);
  
  dualCreated    = false;
  connectCreated = true;
  
  
  //Clearing_______________________________________________________________________________________
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
  
  interfaceToElement.clear();
  interfaceToEdge.clear();
  
  subCellToElement.clear();
  subCellToNode.clear();  
  
  
  
  //Num local and global elements__________________________________________________________________
  UInt numLocalElements = grid2d->getNumElements();
  UInt numLocalEdges    = grid2d->getNumEdges();
  
  pMapGlobalManip<NODEMAP> vrMapManip(commDev);
  pMapGlobalManip<ELMAP>   elMapManip(commDev);
  pMapGlobalManip<ELMAP>   edMapManip(commDev);
  
  pMap<NODEMAP> vrMap = grid2d->getNodes().getMapRef();
  pMap<ELMAP>   elMap = grid2d->getElements().getRowMap();
  pMap<ELMAP>   edMap = grid2d->getEdges().getRowMap();
  
  UInt numGlobalNodes    = vrMapManip.sizeG(vrMap);
  UInt numGlobalElements = elMapManip.sizeG(elMap);
  
  
  //Allocations____________________________________________________________________________________
  UInt pid = commDev->rank();
  
  UInt    id1, id2, id3;
  point3d P1, P2, P3, PM;
  
  GEOELEMENT2D subCell(true);
  GEOELEMENT1D interface(true);
  
  UInt numSubCells = 1;
  UInt numInterfaces = 1;
  
  UInt ba, locStray, elGid, geoId, edgeId, numPoints;
  sVect<UInt> p(3), e(3);
  
  ELMAP   elMapItem;
  NODEMAP nodeMapItem;
  
  
  //Nodes construction_____________________________________________________________________________
  dualPoints = grid2d->getNodes();   //Orginal nodes
  numPoints  = dualPoints.size() + 1;
  
  for(UInt i=1; i <= numLocalElements; ++i)    //Nodi baricentrici
  {
    id1 = grid2d->getElementL(i).getCid(1);
    id2 = grid2d->getElementL(i).getCid(2);
    id3 = grid2d->getElementL(i).getCid(3);
    
    P1  = grid2d->getNodeL(id1);
    P2  = grid2d->getNodeL(id2);
    P3  = grid2d->getNodeL(id3);
    
    PM  = (P1 + P2 + P3) / 3.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes +  grid2d->getElements().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  for(UInt i=1; i <= numLocalEdges; ++i)       //Nodi degli spigoli 
  {
    id1 = grid2d->getEdgeL(i).getCid(1);
    id2 = grid2d->getEdgeL(i).getCid(2);
    
    P1  = grid2d->getNodeL(id1);
    P2  = grid2d->getNodeL(id2);
    
    PM = (P1 + P2) / 2.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes + numGlobalElements + grid2d->getEdges().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  
  //Sub Cells construction_________________________________________________________________________
  for(UInt i=1; i <= numLocalElements; ++i)
  {
    elGid = grid2d->getElements().getRowMapL(i).getGid();
    
    //Triangle nodes
    for(UInt j=1; j <= 3; ++j) 
    {
      p(j) = grid2d->getElementL(i).getCid(j);
    }
    
    //Hexa barycenter
    ba = getLocalNodeId(2,i);
    
    //Nodes on the Hexa edges
    for(UInt j=1; j <= 3; ++j)
    {
      e(j) = connectGrid2d->getElementToEdge(i,j);
      e(j) = getLocalNodeId(3,e(j));
    }
    
    //Subcells creation
    locStray = 1;
    elGid    = grid2d->getElements().getRowMapL(i).getGid();
    geoId    = grid2d->getElements().getItemL(i).getGeoId();
    
    elMapItem.setPid(pid);
    subCell.setGeoId(geoId);
    
    
    //Subcell 1--------------------------------------------
    subCell.setCid(1,p(1));
    subCell.setCid(2,e(1));
    subCell.setCid(3,ba);
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(6 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(1));
    
    
    //Subcell 2--------------------------------------------
    subCell.setCid(1,p(2));
    subCell.setCid(2,ba);
    subCell.setCid(3,e(1));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(6 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(2));
    
    
    //Subcell 3--------------------------------------------
    subCell.setCid(1,p(2));
    subCell.setCid(2,e(2));
    subCell.setCid(3,ba);
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(6 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(2));
    
    
    //Subcell 4--------------------------------------------
    subCell.setCid(1,p(3));
    subCell.setCid(2,ba);
    subCell.setCid(3,e(2));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(6 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(3));
    
    
    //Subcell 5--------------------------------------------
    subCell.setCid(1,p(3));
    subCell.setCid(2,e(3));
    subCell.setCid(3,ba);
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(6 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(3));
    
    
    //Subcell 6--------------------------------------------
    subCell.setCid(1,p(1));
    subCell.setCid(2,ba);
    subCell.setCid(3,e(3));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(6 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(1));
    
    
    //Interfaces creation__________________________________________________________________________
    locStray = 1;
    
    //Interface 1------------------------------------------
    edgeId = 1;
    
    interface.setCid(1,ba);
    interface.setCid(2,e(1));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(3 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid2d->getElementToEdge(i,edgeId));
    interfaceToElement.push_back(elMapItem, i);
    
    
    //Interface 2------------------------------------------
    edgeId = 2;
    
    interface.setCid(1,ba);
    interface.setCid(2,e(2));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(3 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid2d->getElementToEdge(i,edgeId));
    interfaceToElement.push_back(elMapItem, i);
    
    
    //Interface 3------------------------------------------
    edgeId = 3;
    
    interface.setCid(1,ba);
    interface.setCid(2,e(3));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(3 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid2d->getElementToEdge(i,edgeId));
    interfaceToElement.push_back(elMapItem, i);
  }
  
  
  //Map fixing for subCells
  pMap<ELMAP> tempSubCellsMap = subCells.getRowMap();

  subCells.setColMap(dualPoints.getMapRef());
  subCells.updateRowFinder();
  subCells.updateColFinder();
  
  subCellToElement.setMap(tempSubCellsMap);
  subCellToElement.updateFinder();

  subCellToNode.setMap(tempSubCellsMap);
  subCellToNode.updateFinder();
  
  
  //MapFixing for interfaces
  pMap<ELMAP> tempInterfacesMap = interfaces.getRowMap();

  interfaces.setColMap(dualPoints.getMapRef());
  interfaces.updateRowFinder();
  interfaces.updateColFinder();
  
  interfaceToEdge.setMap(tempInterfacesMap);
  interfaceToEdge.updateFinder();
  
  interfaceToElement.setMap(tempInterfacesMap);
  interfaceToElement.updateFinder();
  
  
  //Clearing_______________________________________________________________________________________
  buildFiniteVolumeData();
  
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
}

template<typename ELMAP, typename NODEMAP> 
void
dual2d<linearTriangle,ELMAP,NODEMAP>::
buildFiniteVolumeData()
{
  //Assert_________________________________________________________________________________________
  assert(commDevLoaded);
  assert(grid2dLoaded);
  assert(connectGrid2dLoaded);
  
  finiteVolumeCreated = true;
  
  //Clearing_______________________________________________________________________________________
  cellVolume.clear();
  edgeN.clear();
  edgeSurf.clear();
  
  subCellVolume.clear();
  subCellN.clear();
  
  
  //Resizing_______________________________________________________________________________________
  edgeN.resize(grid2d->getNumEdges());
  edgeSurf.resize(grid2d->getNumEdges());
  cellVolume.resize(grid2d->getNumNodes());
  
  subCellVolume.resize(subCells.size());
  subCellN.resize(subCells.size());
  
  
  //Load maps______________________________________________________________________________________
  edgeN.setMap(grid2d->getEdges().getRowMap());
  edgeSurf.setMap(grid2d->getEdges().getRowMap());
  cellVolume.setMap(grid2d->getNodes().getMapRef());
  
  subCellVolume.setMap(subCells.getMapRef());
  subCellN.setMap(subCells.getMapRef());
  
  
  //Allocate_______________________________________________________________________________________
  UInt cell, el, ed;
  UInt id1, id2, id3;
  point3d P1, P2, P3;
  point3d N, L, edgeL, elementL;
  
  
  //3d support_____________________________________________________________________________________
  for(UInt i=1; i <= subCells.size(); ++i)
  {
    //Surface interface data
    id1 = subCells(i).getCid(1);
    id2 = subCells(i).getCid(2);
    id3 = subCells(i).getCid(3);
    
    P1 = dualPoints(id1);
    P2 = dualPoints(id2);
    P3 = dualPoints(id3);
    
    N = (P2-P1) ^ (P3-P1);
    
    //Element data
    el = subCellToElement(i);
    
    id1 = grid2d->getElementL(el).getCid(1);
    id2 = grid2d->getElementL(el).getCid(2);
    id3 = grid2d->getElementL(el).getCid(3);
    
    P1 = grid2d->getNodeL(id1);
    P2 = grid2d->getNodeL(id2);
    P3 = grid2d->getNodeL(id3);
    
    elementL = (P2-P1) ^ (P3-P1);
    elementL = elementL / elementL.norm2();
    
    //Creazione volume interfaccia
    subCellVolume(i) = 0.5 * N.norm2();
    
    //Creazione normale
    subCellN(i) = elementL;
  }

  
  //2d Support_____________________________________________________________________________________
  for(UInt i=1; i <= interfaces.size(); ++i)
  {
    //Interface data
    ed = interfaceToEdge(i);
    el = interfaceToElement(i);
    
    //Normale elemento
    id1 = grid2d->getElementL(el).getCid(1);
    id2 = grid2d->getElementL(el).getCid(2);
    id3 = grid2d->getElementL(el).getCid(3);
    
    P1 = grid2d->getNodeL(id1);
    P2 = grid2d->getNodeL(id2);
    P3 = grid2d->getNodeL(id3);
    
    N = (P2-P1) ^ (P3-P1);
    N /= N.norm2();
    
    //Normale all' interfaccia
    id1 = interfaces(i).getCid(1);
    id2 = interfaces(i).getCid(2);
    
    P1 = dualPoints(id1);
    P2 = dualPoints(id2);
    
    L = (P2-P1) ^ N;
    
    //Direzione dello spigolo
    id1 = grid2d->getEdgeL(ed).getCid(1);
    id2 = grid2d->getEdgeL(ed).getCid(2);
    
    P1 = grid2d->getNodeL(id1);
    P2 = grid2d->getNodeL(id2);
    
    edgeL = P2 - P1;
    edgeL = edgeL / edgeL.norm2();
    
    //Creazione normale
    if((L*edgeL) >= 0.0)
    { edgeN(ed) = edgeN(ed) + L; }
    else
    { edgeN(ed) = edgeN(ed) - L; }
  }
  
  
  //Normals resizing_______________________________________________________________________________
  assert(grid2d->getNumEdges() == edgeN.size());
  for(UInt i=1; i <= grid2d->getNumEdges(); ++i)
  {
    edgeSurf(i) = edgeN(i).norm2();
    edgeN(i)    = edgeN(i) / edgeN(i).norm2();
  }
  
  
  //Calcolo dei volumi di cella____________________________________________________________________
  for(UInt i=1; i <= subCells.size(); ++i)
  {
    cell = subCellToNode(i);
    
    id1 = subCells(i).getCid(1);
    id2 = subCells(i).getCid(2);
    id3 = subCells(i).getCid(3);
    
    P1 = dualPoints(id1);
    P2 = dualPoints(id2);
    P3 = dualPoints(id3);
    
    N = (P2-P1)^(P3-P1);
    
    cellVolume(cell) += 0.5 * N.norm2();
  }
}

template<typename ELMAP, typename NODEMAP> 
void
dual2d<linearTriangle,ELMAP,NODEMAP>::
clearDual()
{
  dualCreated = false;
  
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
}

template<typename ELMAP, typename NODEMAP> 
void
dual2d<linearTriangle,ELMAP,NODEMAP>::
clearConnect()
{
  connectCreated = false;
  
  subCellToElement.clear();
  subCellToNode.clear();  
  interfaceToEdge.clear();
  interfaceToElement.clear();
}

template<typename ELMAP, typename NODEMAP> 
void
dual2d<linearTriangle,ELMAP,NODEMAP>::
clearFiniteVolume()
{
  finiteVolumeCreated = false;
  
  cellVolume.clear();
  subCellVolume.clear();
  subCellN.clear();
  edgeN.clear();
  edgeSurf.clear();
}      

template<typename ELMAP, typename NODEMAP> 
void
dual2d<linearTriangle,ELMAP,NODEMAP>::
clear()
{
  dualCreated         = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
  
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
  
  subCellToElement.clear();
  subCellToNode.clear();  
  interfaceToEdge.clear();
  interfaceToElement.clear();
  
  cellVolume.clear();
  subCellVolume.clear();
  subCellN.clear();
  edgeN.clear();
  edgeSurf.clear();
}



//_______________________________________________________________________________________________________
// LINEAR QUAD SPECILIZATION
//-------------------------------------------------------------------------------------------------------
/*! Support class for the implementation of vertex-centered finite volume schemes. Every quad is decomposed in sub-cells and then aggregated to the mesh vertices.
Also the interfaces between the cells are created.

<b>Data</b>
<ol>
<li> dualPoints = the dual points: they are the union of the original grid vertices and of the barycentric nodes of the elements and of the edges;
<li> dualCells = the grid subcells. Every quad is divided into 4 subcells. Every subcell is a right-hand quad. Every subcell is linked to a vertex of the 
original triangle, therefore each of the original 4 vertices of the quad has one subcell. Every subcell is built joining a vertex of the original 
quad to as couple of edges midpoints and the barycenter of the quad. The order of the nodes is decided in order to create a right-hand quad however the first
node is always a node of the original grid;
<li> interfaces = the interfaces between the subcells. The finite volume method requires the definition of interfaces between cells. In our case the cell is the
union of the subcells sharing a given vertex of the original mesh. In the quad case the interfaces are made of lines. Denoting with p1, p2, the nodes
of the interface the node p1 is on ad edge and the node p2 is the barycenter of the element.
</ol>

<b>Connectivity data</b>
<ol>
<li> subCellToElement: the connectivity of the subCells to the father element;
<li> subCellToNode: the connetivity of the subCell to a node;
<li> interfaceToEdge: the connectivity of the interface to the edge;
<li> interfaceToElement: the connectivity of the interface to the element
</ol>
  
<b>Finite volume data</b>
<ol>
<li> cellVolume: the aggregated volume of the cell;
<li> subCellVolume: the volume of each sub cell
<li> subCellN: the orthogonal vector to the plane of the subcell
<li> edgeN: the mean normal vector to the interface. The normal is obtained as a weighted normal vector of the nornals of the quads that compose the interface.
Every interface is associated to an edge.
<li> edgeSurf: the surface of the interface
</ol>
*/
template<typename ELMAP, typename NODEMAP>
class dual2d<linearQuad,ELMAP,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef linearQuad              GEOSHAPE2D;
    typedef linearQuad::GEOBSHAPE   GEOSHAPE1D;
    typedef geoElement<GEOSHAPE2D>  GEOELEMENT2D;
    typedef geoElement<GEOSHAPE1D>  GEOELEMENT1D;
    
    typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>    MESH2D;
    typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP> CONNECT2D;
    
    typedef pVect<point3d,NODEMAP>              NODESVECT;
    typedef pGraph<GEOELEMENT2D,ELMAP,NODEMAP>  GRAPH2D;
    typedef pGraph<GEOELEMENT1D,ELMAP,NODEMAP>  GRAPH1D;
    
    typedef pVect<Real,ELMAP>       CELLVOLUMESVECT;
    typedef pVect<point3d,ELMAP>    EDGENSVECT;
    typedef pVect<point3d,ELMAP>    SUBCELLNSVECT;
    typedef pVect<Real,ELMAP>       EDGESSURFVECT;
    
    typedef pVect<UInt,ELMAP>       SUBCELL_TO_ELEMENT;
    typedef pVect<UInt,ELMAP>       SUBCELL_TO_NODE;
    typedef pVect<UInt,ELMAP>       INTERFACE_TO_EDGE;
    typedef pVect<UInt,ELMAP>       INTERFACE_TO_ELEMENT;
    //@}
      
    /*! @name Internal flags */ //@{ 
  public:
    bool commDevLoaded;
    bool grid2dLoaded;
    bool connectGrid2dLoaded;
    bool dualCreated;
    bool connectCreated;
    bool finiteVolumeCreated;
    //@}
      
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<communicator> commDev;
    Teuchos::RCP<MESH2D>       grid2d;
    Teuchos::RCP<CONNECT2D>    connectGrid2d;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    geoMapInterface<GEOSHAPE2D>  refShape2d;
    NODESVECT dualPoints;
    GRAPH2D   subCells;
    GRAPH1D   interfaces;
    //@}
    
     /*! @name Finite Volumes data */ //@{ 
  public:
    CELLVOLUMESVECT cellVolume;
    CELLVOLUMESVECT subCellVolume;
    SUBCELLNSVECT   subCellN;
    EDGENSVECT      edgeN;
    EDGESSURFVECT   edgeSurf;
    //@}
    
    /*! @name Internal connectivity data */ //@{ 
  public:
    SUBCELL_TO_ELEMENT   subCellToElement;
    SUBCELL_TO_NODE      subCellToNode;  
    INTERFACE_TO_EDGE    interfaceToEdge;
    INTERFACE_TO_ELEMENT interfaceToElement;
    //@}
    
    /*! @name Constructors */ //@{ 
  public:
    dual2d();
    dual2d(const Teuchos::RCP<communicator> & CommDev);
    dual2d(communicator & CommDev);
    void setMesh2d(const Teuchos::RCP<MESH2D> & Grid2d);
    void setMesh2d(const MESH2D & Grid2d);
    void setConnect2d(const Teuchos::RCP<CONNECT2D> & Connect2d);
    void setConnect2d(const CONNECT2D & Connect2d);
    //@}
    
    /*! @name Funzioni interne e buildup */ //@{ 
  public:    
    /*! The dual nodes is the union of the original nodes + the elements barycentric nodes + the edges ones.
    Given the id of a geoItem (such as an element, face, etc) this function produces the local id of the dual node associated.
    \param flag = 1 the original vertices, = 2 the element barycenter, = 3 the edge one
    \param lid the local numbering of the geoItem (such as an element, face, etc) */
    UInt getLocalNodeId(const UInt & flag, const UInt & lid) const;
      
    /*! Build dual mesh */
    void buildDualMesh();
      
    /*! Build finite volume data. Requires an explicit construction of the dual mesh */
    void buildFiniteVolumeData();
    
    /*! Build finite volume data without explicitly constructing the dual grid.
    The connections and the finite volume data are build. No global ownership is built*/
    //TODO : this function should be reimplemented and enhanced when the data needed by the 
    // finite volume solvers will become clear - this is due in 2013
    void buildFiniteVolumeDataOnly();
    
    /*! Clear the dual data only */
    void clearDual();
    
     /*! Clear cennection data only */
    void clearConnect();
      
    /*! Clear the finite volume data only */
    void clearFiniteVolume();
      
    /*! Clear all the data */
    void clear();  
    //@}
};

template<typename ELMAP, typename NODEMAP>
dual2d<linearQuad,ELMAP,NODEMAP>::
dual2d()
{
  commDevLoaded       = false;
  grid2dLoaded        = false;
  connectGrid2dLoaded = false;
  dualCreated         = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
}

template<typename ELMAP, typename NODEMAP>
dual2d<linearQuad,ELMAP,NODEMAP>::
dual2d(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded       = true;
  grid2dLoaded        = false;
  connectGrid2dLoaded = false;
  dualCreated         = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
  
  commDev = CommDev;
}

template<typename ELMAP, typename NODEMAP>
dual2d<linearQuad,ELMAP,NODEMAP>::
dual2d(communicator & CommDev)
{
  commDevLoaded       = true;
  grid2dLoaded        = false;
  connectGrid2dLoaded = false;
  dualCreated         = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
  
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename ELMAP, typename NODEMAP>
void
dual2d<linearQuad,ELMAP,NODEMAP>::
setMesh2d(const Teuchos::RCP<MESH2D> & Grid2d)
{
  grid2dLoaded = true;
  grid2d       = Grid2d;
}

template<typename ELMAP, typename NODEMAP>
void
dual2d<linearQuad,ELMAP,NODEMAP>::
setMesh2d(const MESH2D & Grid2d)
{
  grid2dLoaded = true;
  grid2d       = Teuchos::rcp(new MESH2D(Grid2d));
}

template<typename ELMAP, typename NODEMAP>
void
dual2d<linearQuad,ELMAP,NODEMAP>::
setConnect2d(const Teuchos::RCP<CONNECT2D> & Connect2d)
{
  connectGrid2dLoaded = true;
  connectGrid2d       = Connect2d;
}

template<typename ELMAP, typename NODEMAP>
void
dual2d<linearQuad,ELMAP,NODEMAP>::
setConnect2d(const CONNECT2D & Connect2d)
{
  connectGrid2dLoaded = true;
  connectGrid2d       = Teuchos::rcp(new CONNECT2D(Connect2d));
}

template<typename ELMAP, typename NODEMAP> 
UInt
dual2d<linearQuad,ELMAP,NODEMAP>::
getLocalNodeId(const UInt & flag, const UInt & lid) const
{ 
  return( UInt(flag > 1) * grid2d->getNumNodes()
        + UInt(flag > 2) * grid2d->getNumElements()
        + lid );
}

template<typename ELMAP, typename NODEMAP> 
void
dual2d<linearQuad,ELMAP,NODEMAP>::
buildDualMesh()
{
  //Assert_________________________________________________________________________________________
  assert(commDevLoaded);
  assert(grid2dLoaded);
  assert(connectGrid2dLoaded);
  
  dualCreated    = true;
  connectCreated = true;
  
  
  //Clearing_______________________________________________________________________________________
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
  
  interfaceToElement.clear();
  interfaceToEdge.clear();
  
  subCellToElement.clear();
  subCellToNode.clear();  
  
  
  //Num local and global elements__________________________________________________________________
  UInt numLocalElements = grid2d->getNumElements();
  UInt numLocalEdges    = grid2d->getNumEdges();
  
  pMapGlobalManip<NODEMAP> vrMapManip(commDev);
  pMapGlobalManip<ELMAP>   elMapManip(commDev);
  pMapGlobalManip<ELMAP>   edMapManip(commDev);
  
  pMap<NODEMAP> vrMap = grid2d->getNodes().getMapRef();
  pMap<ELMAP>   elMap = grid2d->getElements().getRowMap();
  pMap<ELMAP>   edMap = grid2d->getEdges().getRowMap();
  
  UInt numGlobalNodes    = vrMapManip.sizeG(vrMap);
  UInt numGlobalElements = elMapManip.sizeG(elMap);
  
  
  //Allocations____________________________________________________________________________________
  UInt pid = commDev->rank();
  
  UInt    id1, id2, id3, id4;
  point3d P1, P2, P3, P4, PM;
  
  GEOELEMENT2D subCell(true);
  GEOELEMENT1D interface(true);
  
  UInt numSubCells = 1;
  UInt numInterfaces = 1;
  
  UInt ba, locStray, elGid, geoId, edgeId, numPoints;
  sVect<UInt> p(4), e(4);
  
  ELMAP   elMapItem;
  NODEMAP nodeMapItem;
  
  
  //Nodes construction_____________________________________________________________________________
  dualPoints = grid2d->getNodes();   //Orginal nodes
  numPoints  = dualPoints.size() + 1;
  
  for(UInt i=1; i <= numLocalElements; ++i)    //Nodi baricentrici
  {
    id1 = grid2d->getElementL(i).getCid(1);
    id2 = grid2d->getElementL(i).getCid(2);
    id3 = grid2d->getElementL(i).getCid(3);
    id4 = grid2d->getElementL(i).getCid(4);
    
    P1  = grid2d->getNodeL(id1);
    P2  = grid2d->getNodeL(id2);
    P3  = grid2d->getNodeL(id3);
    P4  = grid2d->getNodeL(id4);
    
    PM  = (P1 + P2 + P3 + P4) / 4.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes +  grid2d->getElements().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  for(UInt i=1; i <= numLocalEdges; ++i)       //Nodi degli spigoli 
  {
    id1 = grid2d->getEdgeL(i).getCid(1);
    id2 = grid2d->getEdgeL(i).getCid(2);
    
    P1  = grid2d->getNodeL(id1);
    P2  = grid2d->getNodeL(id2);
    
    PM = (P1 + P2) / 2.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes + numGlobalElements + grid2d->getEdges().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  //Map ownership fixing
  pMap<NODEMAP> tempNodeMap = dualPoints.getMapRef();
  
  nodesMapFixer(tempNodeMap,*commDev);
  dualPoints.setMap(tempNodeMap);
  dualPoints.updateFinder();
  
  //Checking
  pVectGlobalManip<point3d,NODEMAP> dualNodesChecker(commDev);
  assert(dualNodesChecker.check(dualPoints));
  
  
  
  //Sub Cells construction_________________________________________________________________________
  for(UInt i=1; i <= numLocalElements; ++i)
  {
    elMapItem = grid2d->getElements().getRowMapL(i);
    elGid     = grid2d->getElements().getRowMapL(i).getGid();
     
    //Triangle nodes
    for(UInt j=1; j <= 4; ++j) 
    {
      p(j) = grid2d->getElementL(i).getCid(j);
    }
    
    //Hexa barycenter
    ba = getLocalNodeId(2,i);
    
    //Nodes on the Hexa edges
    for(UInt j=1; j <= 4; ++j)
    {
      e(j) = connectGrid2d->getElementToEdge(i,j);
      e(j) = getLocalNodeId(3,e(j));
    }
    
    //Subcells creation
    locStray = 1;
    elGid    = grid2d->getElements().getRowMapL(i).getGid();
    geoId    = grid2d->getElements().getItemL(i).getGeoId();
    
    elMapItem.setPid(pid);
    subCell.setGeoId(geoId);
    
    //Subcell 1--------------------------------------------
    subCell.setCid(1,p(1));
    subCell.setCid(2,e(1));
    subCell.setCid(3,ba);
    subCell.setCid(4,e(4));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(1));
    
    //Subcell 2--------------------------------------------
    subCell.setCid(1,p(2));
    subCell.setCid(2,e(2));
    subCell.setCid(3,ba);
    subCell.setCid(4,e(1));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(2));
    
    //Subcell 3--------------------------------------------
    subCell.setCid(1,p(3));
    subCell.setCid(2,e(3));
    subCell.setCid(3,ba);
    subCell.setCid(4,e(2));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(3));
    
    //Subcell 4--------------------------------------------
    subCell.setCid(1,p(4));
    subCell.setCid(2,e(4));
    subCell.setCid(3,ba);
    subCell.setCid(4,e(3));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(4));
    
    
    //Interfaces creation__________________________________________________________________________
    locStray = 1;
    
    //Interface 1------------------------------------------
    edgeId = 1;
    
    interface.setCid(1,e(1));
    interface.setCid(2,ba);
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid2d->getElementToEdge(i,edgeId));
    interfaceToElement.push_back(elMapItem, i);
    
    //Interface 2------------------------------------------
    edgeId = 2;
    
    interface.setCid(1,e(2));
    interface.setCid(2,ba);
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid2d->getElementToEdge(i,edgeId));
    interfaceToElement.push_back(elMapItem, i);
    
    //Interface 3------------------------------------------
    edgeId = 3;
    
    interface.setCid(1,e(3));
    interface.setCid(2,ba);
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid2d->getElementToEdge(i,edgeId));
    interfaceToElement.push_back(elMapItem, i);
    
    //Interface 4------------------------------------------
    edgeId = 4;
    
    interface.setCid(1,e(4));
    interface.setCid(2,ba);
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid2d->getElementToEdge(i,edgeId));
    interfaceToElement.push_back(elMapItem, i);
  }
  
  
  //Map fixing for subCells
  subCells.setColMap(dualPoints.getMapRef());
  subCells.updateRowFinder();
  subCells.updateColFinder();
  
  subCellToElement.setMap(subCells.getRowMap());
  subCellToElement.updateFinder();

  subCellToNode.setMap(subCells.getRowMap());
  subCellToNode.updateFinder();
  
  
  //MapFixing for interfaces
  interfaces.setColMap(dualPoints.getMapRef());
  interfaces.updateRowFinder();
  interfaces.updateColFinder();
  
  interfaceToEdge.setMap(interfaces.getRowMap());
  interfaceToEdge.updateFinder();
  
  interfaceToElement.setMap(interfaces.getRowMap());
  interfaceToElement.updateFinder();
  
  //Checking interfaces map
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> interfacesCheker(commDev);
  assert(interfacesCheker.check(interfaces));
}

template<typename ELMAP, typename NODEMAP> 
void
dual2d<linearQuad,ELMAP,NODEMAP>::
buildFiniteVolumeDataOnly()
{
  //Assert_________________________________________________________________________________________
  assert(commDevLoaded);
  assert(grid2dLoaded);
  assert(connectGrid2dLoaded);
  
  connectCreated      = true;
  finiteVolumeCreated = true;
  
  //Clearing_______________________________________________________________________________________
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
  
  interfaceToElement.clear();
  interfaceToEdge.clear();
  
  subCellToElement.clear();
  subCellToNode.clear();  
  
  
  //Num local and global elements__________________________________________________________________
  UInt numLocalElements = grid2d->getNumElements();
  UInt numLocalEdges    = grid2d->getNumEdges();
  
  pMapGlobalManip<NODEMAP> vrMapManip(commDev);
  pMapGlobalManip<ELMAP>   elMapManip(commDev);
  pMapGlobalManip<ELMAP>   edMapManip(commDev);
  
  pMap<NODEMAP> vrMap = grid2d->getNodes().getMapRef();
  pMap<ELMAP>   elMap = grid2d->getElements().getRowMap();
  pMap<ELMAP>   edMap = grid2d->getEdges().getRowMap();
  
  UInt numGlobalNodes    = vrMapManip.sizeG(vrMap);
  UInt numGlobalElements = elMapManip.sizeG(elMap);
  
  
  //Allocations____________________________________________________________________________________
  UInt pid = commDev->rank();
  
  UInt    id1, id2, id3, id4;
  point3d P1, P2, P3, P4, PM;
  
  GEOELEMENT2D subCell(true);
  GEOELEMENT1D interface(true);
  
  UInt numSubCells = 1;
  UInt numInterfaces = 1;
  
  UInt ba, locStray, elGid, geoId, edgeId, numPoints;
  sVect<UInt> p(4), e(4);
  
  ELMAP   elMapItem;
  NODEMAP nodeMapItem;
  
  
  //Nodes construction_____________________________________________________________________________
  dualPoints = grid2d->getNodes();   //Orginal nodes
  numPoints  = dualPoints.size() + 1;
  
  for(UInt i=1; i <= numLocalElements; ++i)    //Nodi baricentrici
  {
    id1 = grid2d->getElementL(i).getCid(1);
    id2 = grid2d->getElementL(i).getCid(2);
    id3 = grid2d->getElementL(i).getCid(3);
    id4 = grid2d->getElementL(i).getCid(4);
    
    P1  = grid2d->getNodeL(id1);
    P2  = grid2d->getNodeL(id2);
    P3  = grid2d->getNodeL(id3);
    P4  = grid2d->getNodeL(id4);
    
    PM  = (P1 + P2 + P3 + P4) / 4.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes +  grid2d->getElements().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  for(UInt i=1; i <= numLocalEdges; ++i)       //Nodi degli spigoli 
  {
    id1 = grid2d->getEdgeL(i).getCid(1);
    id2 = grid2d->getEdgeL(i).getCid(2);
    
    P1  = grid2d->getNodeL(id1);
    P2  = grid2d->getNodeL(id2);
    
    PM = (P1 + P2) / 2.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes + numGlobalElements + grid2d->getEdges().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  
  //Sub Cells construction_________________________________________________________________________
  for(UInt i=1; i <= numLocalElements; ++i)
  {
    elGid = grid2d->getElements().getRowMapL(i).getGid();
    
    //Triangle nodes
    for(UInt j=1; j <= 4; ++j) 
    {
      p(j) = grid2d->getElementL(i).getCid(j);
    }
    
    //Hexa barycenter
    ba = getLocalNodeId(2,i);
    
    //Nodes on the Hexa edges
    for(UInt j=1; j <= 4; ++j)
    {
      e(j) = connectGrid2d->getElementToEdge(i,j);
      e(j) = getLocalNodeId(3,e(j));
    }
    
    //Subcells creation
    locStray = 1;
    elGid    = grid2d->getElements().getRowMapL(i).getGid();
    geoId    = grid2d->getElements().getItemL(i).getGeoId();
    
    elMapItem.setPid(pid);
    subCell.setGeoId(geoId);
    
    //Subcell 1--------------------------------------------
    subCell.setCid(1,p(1));
    subCell.setCid(2,e(1));
    subCell.setCid(3,ba);
    subCell.setCid(4,e(4));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(1));
    
    //Subcell 2--------------------------------------------
    subCell.setCid(1,p(2));
    subCell.setCid(2,e(2));
    subCell.setCid(3,ba);
    subCell.setCid(4,e(1));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(2));
    
    //Subcell 3--------------------------------------------
    subCell.setCid(1,p(3));
    subCell.setCid(2,e(3));
    subCell.setCid(3,ba);
    subCell.setCid(4,e(2));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(3));
    
    //Subcell 4--------------------------------------------
    subCell.setCid(1,p(4));
    subCell.setCid(2,e(4));
    subCell.setCid(3,ba);
    subCell.setCid(4,e(3));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(4));
    
    
    //Interfaces creation__________________________________________________________________________
    locStray = 1;
    
    //Interface 1------------------------------------------
    edgeId = 1;
    
    interface.setCid(1,e(1));
    interface.setCid(2,ba);
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid2d->getElementToEdge(i,edgeId));
    interfaceToElement.push_back(elMapItem, i);
    
    //Interface 2------------------------------------------
    edgeId = 2;
    
    interface.setCid(1,e(2));
    interface.setCid(2,ba);
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid2d->getElementToEdge(i,edgeId));
    interfaceToElement.push_back(elMapItem, i);
    
    //Interface 3------------------------------------------
    edgeId = 3;
    
    interface.setCid(1,e(3));
    interface.setCid(2,ba);
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid2d->getElementToEdge(i,edgeId));
    interfaceToElement.push_back(elMapItem, i);
    
    //Interface 4------------------------------------------
    edgeId = 4;
    
    interface.setCid(1,e(4));
    interface.setCid(2,ba);
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(4 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid2d->getElementToEdge(i,edgeId));
    interfaceToElement.push_back(elMapItem, i);
  }
  
  
  //Map fixing for subCells
  pMap<ELMAP> tempSubCellsMap = subCells.getRowMap();

  subCells.setColMap(dualPoints.getMapRef());
  subCells.updateRowFinder();
  subCells.updateColFinder();
  
  subCellToElement.setMap(tempSubCellsMap);
  subCellToElement.updateFinder();

  subCellToNode.setMap(tempSubCellsMap);
  subCellToNode.updateFinder();
  
  
  //MapFixing for interfaces
  pMap<ELMAP> tempInterfacesMap = interfaces.getRowMap();

  interfaces.setColMap(dualPoints.getMapRef());
  interfaces.updateRowFinder();
  interfaces.updateColFinder();
  
  interfaceToEdge.setMap(tempInterfacesMap);
  interfaceToEdge.updateFinder();
  
  interfaceToElement.setMap(tempInterfacesMap);
  interfaceToElement.updateFinder();
  
  
  //Clearing_______________________________________________________________________________________
  buildFiniteVolumeData();
  
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
}

template<typename ELMAP, typename NODEMAP> 
void
dual2d<linearQuad,ELMAP,NODEMAP>::
buildFiniteVolumeData()
{
  //Assert_________________________________________________________________________________________
  assert(commDevLoaded);
  assert(grid2dLoaded);
  assert(connectGrid2dLoaded);
  
  finiteVolumeCreated = true;
  
  //Clearing_______________________________________________________________________________________
  cellVolume.clear();
  edgeN.clear();
  edgeSurf.clear();
  
  subCellVolume.clear();
  subCellN.clear();
  
  
  //Resizing_______________________________________________________________________________________
  edgeN.resize(grid2d->getNumEdges());
  edgeSurf.resize(grid2d->getNumEdges());
  cellVolume.resize(grid2d->getNumNodes());
  
  subCellVolume.resize(subCells.size());
  subCellN.resize(subCells.size());
  
  
  //Load maps______________________________________________________________________________________
  edgeN.setMap(grid2d->getEdges().getRowMap());
  edgeSurf.setMap(grid2d->getEdges().getRowMap());
  cellVolume.setMap(grid2d->getNodes().getMapRef());
  
  subCellVolume.setMap(subCells.getMapRef());
  subCellN.setMap(subCells.getMapRef());
  
  
  //Allocate_______________________________________________________________________________________
  UInt el, ed, cell;
  UInt id1, id2, id3, id4;
  Real surface;
  point3d N, L, P1, P2, edgeL;
  sVect<point3d> points2d(4);
  
  geoMapSupport2d<GEOSHAPE2D> geoSupport2d;
  
  
  //3d support_____________________________________________________________________________________
  for(UInt i=1; i <= subCells.size(); ++i)
  {
    //Surface interface data
    id1 = subCells(i).getCid(1);
    id2 = subCells(i).getCid(2);
    id3 = subCells(i).getCid(3);
    id4 = subCells(i).getCid(4);
    
    points2d(1) = dualPoints(id1);
    points2d(2) = dualPoints(id2);
    points2d(3) = dualPoints(id3);
    points2d(4) = dualPoints(id4);
    
    geoSupport2d.setPoints(points2d);
    geoSupport2d.surfaceNormal<STANDARD,3>(N,surface);
    
    //Creazione volume interfaccia
    subCellVolume(i) = surface;
    
    //Creazione normale
    subCellN(i) = N;
  }
    
  
  //2d Support_____________________________________________________________________________________
  for(UInt i=1; i <= interfaces.size(); ++i)
  {
    //Interface data
    ed = interfaceToEdge(i);
    el = interfaceToElement(i);
    
    //Normale elemento
    id1 = grid2d->getElementL(el).getCid(1);
    id2 = grid2d->getElementL(el).getCid(2);
    id3 = grid2d->getElementL(el).getCid(3);
    id4 = grid2d->getElementL(el).getCid(4);
    
    points2d(1) = grid2d->getNodeL(id1);
    points2d(2) = grid2d->getNodeL(id2);
    points2d(3) = grid2d->getNodeL(id3);
    points2d(4) = grid2d->getNodeL(id4);
    
    geoSupport2d.setPoints(points2d);
    geoSupport2d.surfaceNormal<STANDARD,3>(N,surface);
    
    //Normale all' interfaccia
    id1 = interfaces(i).getCid(1);
    id2 = interfaces(i).getCid(2);
    
    P1 = dualPoints(id1);
    P2 = dualPoints(id2);
    
    L = (P2-P1) ^ N;
    
    //Direzione dello spigolo
    id1 = grid2d->getEdgeL(ed).getCid(1);
    id2 = grid2d->getEdgeL(ed).getCid(2);
    
    P1 = grid2d->getNodeL(id1);
    P2 = grid2d->getNodeL(id2);
    
    edgeL = P2 - P1;
    edgeL = edgeL / edgeL.norm2();
    
    //Creazione normale
    if((L * edgeL) >= 0.0)
    { edgeN(ed) = edgeN(ed) + L; }
    else
    { edgeN(ed) = edgeN(ed) - L; }
  }
  
  
  //Normals resizing_______________________________________________________________________________
  assert(grid2d->getNumEdges() == edgeN.size());
  for(UInt i=1; i <= grid2d->getNumEdges(); ++i)
  {
    edgeSurf(i) = edgeN(i).norm2();
    edgeN(i)    = edgeN(i) / edgeN(i).norm2();
  }
  
  
  //Calcolo dei volumi di cella____________________________________________________________________
  for(UInt i=1; i <= subCells.size(); ++i)
  {
    cell = subCellToNode(i);
    
    id1 = subCells(i).getCid(1);
    id2 = subCells(i).getCid(2);
    id3 = subCells(i).getCid(3);
    id4 = subCells(i).getCid(4);
    
    points2d(1) = dualPoints(id1);
    points2d(2) = dualPoints(id2);
    points2d(3) = dualPoints(id3);
    points2d(4) = dualPoints(id4);
    
    geoSupport2d.setPoints(points2d);    
    cellVolume(cell) += geoSupport2d.volume<STANDARD,3>();
  }
}

template<typename ELMAP, typename NODEMAP> 
void
dual2d<linearQuad,ELMAP,NODEMAP>::
clearDual()
{
  dualCreated = false;
  
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
}

template<typename ELMAP, typename NODEMAP> 
void
dual2d<linearQuad,ELMAP,NODEMAP>::
clearConnect()
{
  connectCreated = false;
  
  subCellToElement.clear();
  subCellToNode.clear();  
  interfaceToEdge.clear();
  interfaceToElement.clear();
}

template<typename ELMAP, typename NODEMAP> 
void
dual2d<linearQuad,ELMAP,NODEMAP>::
clearFiniteVolume()
{
  finiteVolumeCreated = false;
  
  cellVolume.clear();
  subCellVolume.clear();
  subCellN.clear();
  edgeN.clear();
  edgeSurf.clear();
}      

template<typename ELMAP, typename NODEMAP> 
void
dual2d<linearQuad,ELMAP,NODEMAP>::
clear()
{
  dualCreated         = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
  
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
  
  subCellToElement.clear();
  subCellToNode.clear();  
  interfaceToEdge.clear();
  interfaceToElement.clear();
  
  cellVolume.clear();
  subCellVolume.clear();
  subCellN.clear();
  edgeN.clear();
  edgeSurf.clear();
}


#endif
