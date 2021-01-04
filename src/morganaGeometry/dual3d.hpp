/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef DUAL3D_HPP
#define DUAL3D_HPP

#include "algorithm"
#include "traitsGeometry.hpp"

#include "mesh3d.hpp"
#include "connect3d.hpp"
#include "geoMapSupport2d.hpp"
#include "geoMapSupport3d.hpp"

using namespace std;


//_______________________________________________________________________________________________________
// DUMMY IMPLEMENTATION
//-------------------------------------------------------------------------------------------------------

/*! Dual mesh 2d topology and geometry - dummy implementation */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class dual3d
{ };



//_______________________________________________________________________________________________________
// LINEAR TETRA SPECILIZATION
//-------------------------------------------------------------------------------------------------------

/*! Dual mesh 3d topology and geometry - \c linearTetra implementation. Support class for the implementation of vertex-centered finite volume schemes.
Every tetrahedra is decomposed in sub-cells and then aggregated to the mesh vertices.
Also the interfaces between the cells are created.

<b>Data</b>
<ol>
<li> dualPoints = the dual points: they are the union of the original grid vertices and of the barycentric nodes of the elements, of the faces and of the edges;
<li> dualCells = the grid subcells. Every tetrahedra is divided into 4 x 6 subcelss. Every subcell is a tetrahedron. Every subcell is linked to a vertex of the 
original tetrahedron, therefore each of the original 4 vertices of the tetrahedron has 6 subcells. Every subcell is built joining a vertex of the original 
tetrahedron to: the element barycenter, the barycenter of a face and the barycenter of an edge (the face should be adjacent with the edge). Given the set of
the nodes p1, p2, p3, p4 of the subcell the node p4 is a vertex of the original tetra and the remaining nodes are arranged so that [(p2-p1)^(p3-p1)] * (p1-p4) > 0; 
<li> interfaces = the interfaces between the subcells. The finite volume method require the definition of interaces between cells. In our case the cell is the
union of the subcells sharing a given vertex of the original mesh. In the tetra case the interfaces are made on triangles. Denoting with p1, p2, p3 the nodes
of the interface the node p1 is on ad edge, the node p2 rests on a face and the node p3 is the barycenter of the element.
</ol>

<b>Connectivity data</b>
<ol>
<li> subCellToElement: the connectivity of the subCells to the father element;
<li> subCellToNode: the connetivity of the subCell to a node;
<li> interfaceToEdge: the connectivity of the interface to the edge;
</ol>
  
<b>Finite volume data</b>
<ol>
<li> cellVolume: the volume of the cell;
<li> edgeN: the mean normal vector to the interface. The normal is obtained as a weighted normal vector of the nornals of the triangles that compose the interface.
Every interface is associated to an edge.
<li> edgeSurf: the surface of the interface
</ol>
*/
template<typename ELMAP, typename NODEMAP>
class dual3d<linearTetra,ELMAP,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef linearTetra             GEOSHAPE3D;
    typedef linearTetra::GEOBSHAPE  GEOSHAPE2D;
    typedef GEOSHAPE2D::GEOBSHAPE   GEOSHAPE1D;
    typedef geoElement<GEOSHAPE3D>  GEOELEMENT3D;
    typedef geoElement<GEOSHAPE2D>  GEOELEMENT2D;
    typedef geoElement<GEOSHAPE1D>  GEOELEMENT1D;
    
    typedef mesh3d<GEOSHAPE3D,ELMAP,NODEMAP>    MESH3D;
    typedef connect3d<GEOSHAPE3D,ELMAP,NODEMAP> CONNECT3D;
    
    typedef pVect<point3d,NODEMAP>              NODESVECT;
    typedef pGraph<GEOELEMENT3D,ELMAP,NODEMAP>  GRAPH3D;
    typedef pGraph<GEOELEMENT2D,ELMAP,NODEMAP>  GRAPH2D;
    
    typedef pVect<Real,ELMAP>       CELLVOLUMESVECT;
    typedef pVect<point3d,ELMAP>    EDGENSVECT;
    typedef pVect<Real,ELMAP>       EDGESSURFVECT;
    
    typedef pVect<UInt,ELMAP>       SUBCELL_TO_ELEMENT;
    typedef pVect<UInt,ELMAP>       SUBCELL_TO_NODE;
    typedef pVect<UInt,ELMAP>       INTERFACE_TO_EDGE;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 3;
    //@}
      
    /*! @name Internal flags */ //@{ 
  public:
    bool commDevLoaded;
    bool grid3dLoaded;
    bool connectGrid3dLoaded;
    bool dualCreated;
    bool connectCreated;
    bool finiteVolumeCreated;
    //@}
      
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<communicator> commDev;
    Teuchos::RCP<MESH3D>       grid3d;
    Teuchos::RCP<CONNECT3D>    connectGrid3d;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    geoMapInterface<GEOSHAPE3D>  refShape3d;
    NODESVECT dualPoints;
    GRAPH3D   subCells;
    GRAPH2D   interfaces;
    //@}
    
     /*! @name Finite Volumes data */ //@{ 
  public:
    CELLVOLUMESVECT cellVolume;
    EDGENSVECT      edgeN;
    EDGESSURFVECT   edgeSurf;
    //@}
    
    /*! @name Internal connectivity data */ //@{ 
  public:
    SUBCELL_TO_ELEMENT subCellToElement;
    SUBCELL_TO_NODE    subCellToNode;  
    INTERFACE_TO_EDGE  interfaceToEdge;
    //@}
    
    /*! @name Constructors */ //@{ 
  public:
    dual3d();
    dual3d(const Teuchos::RCP<communicator> & CommDev);
    dual3d(communicator & CommDev);
    void setMesh3d(const Teuchos::RCP<MESH3D> & Grid3d);
    void setMesh3d(const MESH3D & Grid3d);
    void setConnect3d(const Teuchos::RCP<CONNECT3D> & Connect3d);
    void setConnect3d(const CONNECT3D & Connect3d);
    //@}
    
    /*! @name Funzioni interne e buildup */ //@{ 
  public:    
    /*! The dual nodes is the union of the original nodes + the elements barycentric nodes + the faces barycentric nodes + the edges ones.
    Given the id of a geoItem (such as an element, face, etc) this function produces the local id of the dual node associated.
    \param flag = 1 the original vertices, = 2 the element barycenter, = 3 the face one, = 4 the edge one
    \param lid the local numbering of the geoItem (such as an element, face, etc) */
    UInt getLocalNodeId(const UInt & flag, const UInt & lid) const;
      
    /*! Build dual mesh. Explicitly builds the dual mesh */
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
dual3d<linearTetra,ELMAP,NODEMAP>::
dual3d()
{
  commDevLoaded       = false;
  grid3dLoaded        = false;
  connectGrid3dLoaded = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
  dualCreated         = false;
}

template<typename ELMAP, typename NODEMAP>
dual3d<linearTetra,ELMAP,NODEMAP>::
dual3d(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded       = true;
  grid3dLoaded        = false;
  connectGrid3dLoaded = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
  dualCreated         = false;
  
  commDev = CommDev;
}

template<typename ELMAP, typename NODEMAP>
dual3d<linearTetra,ELMAP,NODEMAP>::
dual3d(communicator & CommDev)
{
  commDevLoaded       = true;
  grid3dLoaded        = false;
  connectGrid3dLoaded = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
  dualCreated         = false;
  
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename ELMAP, typename NODEMAP>
void
dual3d<linearTetra,ELMAP,NODEMAP>::
setMesh3d(const Teuchos::RCP<MESH3D> & Grid3d)
{
  grid3dLoaded = true;
  grid3d       = Grid3d;
}

template<typename ELMAP, typename NODEMAP>
void
dual3d<linearTetra,ELMAP,NODEMAP>::
setMesh3d(const MESH3D & Grid3d)
{
  grid3dLoaded = true;
  grid3d       = Teuchos::rcp(new MESH3D(Grid3d));
}

template<typename ELMAP, typename NODEMAP>
void
dual3d<linearTetra,ELMAP,NODEMAP>::
setConnect3d(const Teuchos::RCP<CONNECT3D> & Connect3d)
{
  connectGrid3dLoaded = true;
  connectGrid3d       = Connect3d;
}

template<typename ELMAP, typename NODEMAP>
void
dual3d<linearTetra,ELMAP,NODEMAP>::
setConnect3d(const CONNECT3D & Connect3d)
{
  connectGrid3dLoaded = true;
  connectGrid3d       = Teuchos::rcp(new CONNECT3D(Connect3d));
}

template<typename ELMAP, typename NODEMAP> 
UInt
dual3d<linearTetra,ELMAP,NODEMAP>::
getLocalNodeId(const UInt & flag, const UInt & lid) const
{ 
  return( UInt(flag > 1) * grid3d->getNumNodes()
        + UInt(flag > 2) * grid3d->getNumElements()
        + UInt(flag > 3) * grid3d->getNumFaces()
        + lid );
}

template<typename ELMAP, typename NODEMAP> 
void
dual3d<linearTetra,ELMAP,NODEMAP>::
buildDualMesh()
{
  //Assert_________________________________________________________________________________________
  assert(commDevLoaded);
  assert(grid3dLoaded);
  assert(connectGrid3dLoaded);
  
  connectCreated = true;
  dualCreated    = true;
  
  
  //Clearing_______________________________________________________________________________________
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
  
  subCellToElement.clear();
  subCellToNode.clear();  
  interfaceToEdge.clear();
  
  
  //Num local and global elements__________________________________________________________________
  UInt numLocalNodes    = grid3d->getNumNodes();
  UInt numLocalElements = grid3d->getNumElements();
  UInt numLocalFaces    = grid3d->getNumFaces();
  UInt numLocalEdges    = grid3d->getNumEdges();
  
  pMapGlobalManip<NODEMAP> vrMapManip(commDev);
  pMapGlobalManip<ELMAP>   elMapManip(commDev);
  pMapGlobalManip<ELMAP>   faMapManip(commDev);
  pMapGlobalManip<ELMAP>   edMapManip(commDev);
  
  pMap<NODEMAP> vrMap = grid3d->getNodes().getMapRef();
  pMap<ELMAP>   elMap = grid3d->getElements().getRowMap();
  pMap<ELMAP>   faMap = grid3d->getFaces().getRowMap();
  pMap<ELMAP>   edMap = grid3d->getEdges().getRowMap();
  
  UInt numGlobalNodes    = vrMapManip.sizeG(vrMap);
  UInt numGlobalElements = elMapManip.sizeG(elMap);
  UInt numGlobalFaces    = faMapManip.sizeG(faMap);
  
  
  //Reserving space________________________________________________________________________________
  dualPoints.reserve(numLocalNodes + numLocalElements + numLocalFaces + numLocalEdges);
  subCells.reserve(24 * numLocalElements);
  subCellToElement.reserve(24 * numLocalElements);
  subCellToNode.reserve(24 * numLocalElements);
  interfaces.reserve(12 * numLocalElements);
  interfaceToEdge.reserve(12 * numLocalElements);
  
  
  //Allocations____________________________________________________________________________________
  UInt pid = commDev->rank();
  
  UInt    id1, id2, id3, id4;
  point3d P1, P2, P3, P4, PM;
  
  GEOELEMENT3D subCell(true);
  sVect<UInt>  f(4), e(6), p(4);
  UInt ba;
  sVect<point3d> F(4), E(6), P(4);
  point3d BA;
  
  Real test;
  UInt no1, no2, edgeId, elGid, locStray, geoId;
  
  GEOELEMENT2D interface(true);
  UInt numSubCells = 1;
  UInt numInterfaces = 1;
  UInt numPoints;
  
  ELMAP   elMapItem;
  NODEMAP nodeMapItem;
  
  
  //Nodes construction_____________________________________________________________________________
  dualPoints = grid3d->getNodes();   //Orginal nodes
  numPoints  = dualPoints.size() + 1;
  
  for(UInt i=1; i <= numLocalElements; ++i)    //Nodi baricentrici
  {
    id1 = grid3d->getElementL(i).getCid(1);
    id2 = grid3d->getElementL(i).getCid(2);
    id3 = grid3d->getElementL(i).getCid(3);
    id4 = grid3d->getElementL(i).getCid(4);
    
    P1  = grid3d->getNodeL(id1);
    P2  = grid3d->getNodeL(id2);
    P3  = grid3d->getNodeL(id3);
    P4  = grid3d->getNodeL(id4);
    
    PM  = (P1 + P2 + P3 + P4) / 4.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes +  grid3d->getElements().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  for(UInt i=1; i <= numLocalFaces; ++i) 	//Nodi delle facce
  {
    id1 = grid3d->getFaceL(i).getCid(1);
    id2 = grid3d->getFaceL(i).getCid(2);
    id3 = grid3d->getFaceL(i).getCid(3);
    
    P1  = grid3d->getNodeL(id1);
    P2  = grid3d->getNodeL(id2);
    P3  = grid3d->getNodeL(id3);
    
    PM  = (P1 + P2 + P3) / 3.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes + numGlobalElements +  grid3d->getFaces().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  for(UInt i=1; i <= numLocalEdges; ++i)       //Nodi degli spigoli 
  {
    id1 = grid3d->getEdgeL(i).getCid(1);
    id2 = grid3d->getEdgeL(i).getCid(2);
    
    P1  = grid3d->getNodeL(id1);
    P2  = grid3d->getNodeL(id2);
    
    PM = (P1 + P2) / 2.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes + numGlobalElements + numGlobalFaces +  grid3d->getEdges().getRowMapL(i).getGid());
    
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
  
  
  
  //Subcells and interfaces construction___________________________________________________________
  for(UInt i=1; i <= numLocalElements; ++i)
  {
    elMapItem = grid3d->getElements().getRowMapL(i);
    elGid     = grid3d->getElements().getRowMapL(i).getGid();
    
    //Tetrahedra nodes
    for(UInt j=1; j <= 4; ++j) 
    {
      p(j) = grid3d->getElementL(i).getCid(j);
      P(j) = dualPoints(p(j));
    }
    
    //Tetrahedra barycenter
    ba = getLocalNodeId(2,i);
    BA = dualPoints(ba);
    
    //Nodes on the tetra faces
    for(UInt j=1; j <= 4; ++j)
    {
      f(j) = connectGrid3d->getElementToFace(i,j);
      f(j) = getLocalNodeId(3,f(j));
      F(j) = dualPoints(f(j));
    }
    
    //Nodes on the tetra edges
    for(UInt j=1; j <= 6; ++j)
    {
      e(j) = connectGrid3d->getElementToEdge(i,j);
      e(j) = getLocalNodeId(4,e(j));
      E(j) = dualPoints(e(j));
    }
    
    //Subcells creation
    locStray = 1;
    elGid    = grid3d->getElements().getRowMapL(i).getGid();
    geoId    = grid3d->getElements().getItemL(i).getGeoId();
    
    for(UInt j=1; j <= 4; ++j) //Face cycle
    {
      for(UInt k=1; k <= 3; ++k) //Edges of the face
      {
	edgeId = refShape3d.faceToEdge(j,k);
	no1    = refShape3d.edgeToPoint(edgeId,1); //first node on the edge
	no2    = refShape3d.edgeToPoint(edgeId,2); //second node on the edge
	
	//First subcell
	test = point3d::dot((E(edgeId) - BA) ^ (F(j) - BA), BA - P(no1));
	
	if(test > 0.0)
	{
	  subCell.setCid(1,ba);
	  subCell.setCid(2,e(edgeId));
	  subCell.setCid(3,f(j));
	  subCell.setCid(4,p(no1));
	}
	else
	{
	  subCell.setCid(1,ba);
	  subCell.setCid(2,f(j));
	  subCell.setCid(3,e(edgeId));
	  subCell.setCid(4,p(no1));
	}
	  
	elMapItem.setPid(pid);
	elMapItem.setLid(numSubCells);
	elMapItem.setGid(24 * (elGid-1) + locStray);
	
	subCell.setGeoId(geoId);
	  
	locStray++;
	numSubCells++;
	  
	subCells.push_back(elMapItem,subCell);
	subCellToElement.push_back(elMapItem,i);
	subCellToNode.push_back(elMapItem,p(no1));
	
	
	//Seconda sottocella 
	test = point3d::dot((F(j) - BA) ^ (E(edgeId) - BA), BA - P(no2));
	
	if(test > 0.0)
	{
	  subCell.setCid(1,ba);
	  subCell.setCid(2,f(j));
	  subCell.setCid(3,e(edgeId));
	  subCell.setCid(4,p(no2));
	}
	else
	{
	  subCell.setCid(1,ba);
	  subCell.setCid(2,e(edgeId));
	  subCell.setCid(3,f(j));
	  subCell.setCid(4,p(no2));
	}
	
	elMapItem.setPid(pid);
	elMapItem.setLid(numSubCells);
	elMapItem.setGid(24 * (elGid-1) + locStray);
	
	subCell.setGeoId(geoId);
	  
	locStray++;
	numSubCells++;
	  
	subCells.push_back(elMapItem,subCell);
	subCellToElement.push_back(elMapItem,i);
	subCellToNode.push_back(elMapItem,p(no2));
      }
    }
    
    
    //Interfaces creation
    locStray = 1;
    
    for(UInt j=1; j <= 4; ++j) //faces
    {
      for(UInt k=1; k <= 3; ++k) //edges on the face
      {
	edgeId = refShape3d.faceToEdge(j,k);
	
	interface.setCid(1,e(edgeId));
	interface.setCid(2,f(j));
	interface.setCid(3,ba);
	interface.setGeoId(geoId);
	
	elMapItem.setPid(pid);
	elMapItem.setLid(numInterfaces);
	elMapItem.setGid(12 * (elGid-1) + locStray);
	
	numInterfaces++;
	locStray++;
	
	interfaces.push_back(elMapItem,interface);
	interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
      }
    }
    
    assert(locStray == 13);
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
  pGraphGlobalManip<GEOELEMENT3D,ELMAP,NODEMAP> subCellsCheker(commDev);
  assert(subCellsCheker.check(subCells));
  
  
  //MapFixing for interfaces
  interfaces.setColMap(dualPoints.getMapRef());
  interfaces.updateRowFinder();
  interfaces.updateColFinder();
  
  interfaceToEdge.setMap(interfaces.getRowMap());
  interfaceToEdge.updateFinder();
  
  //Checking interfaces map
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> interfacesCheker(commDev);
  assert(interfacesCheker.check(interfaces));
}

template<typename ELMAP, typename NODEMAP> 
void
dual3d<linearTetra,ELMAP,NODEMAP>::
buildFiniteVolumeData()
{
  //Assert_________________________________________________________________________________________
  assert(commDevLoaded);
  assert(grid3dLoaded);
  assert(connectGrid3dLoaded);
  
  finiteVolumeCreated = true;
  
  //Clearing_______________________________________________________________________________________
  cellVolume.clear();
  edgeN.clear();
  edgeSurf.clear();
  
  //Resizing_______________________________________________________________________________________
  edgeN.resize(grid3d->getNumEdges());
  edgeSurf.resize(grid3d->getNumEdges());
  cellVolume.resize(grid3d->getNumNodes());
  
  //Load maps______________________________________________________________________________________
  edgeN.setMap(grid3d->getEdges().getRowMap());
  edgeSurf.setMap(grid3d->getEdges().getRowMap());
  cellVolume.setMap(grid3d->getNodes().getMapRef());

  
  
  //Allocate_______________________________________________________________________________________
  UInt cell;
  UInt id1, id2, id3, id4;
  UInt edgeId;
  point3d P1, P2, P3, P4, PM, PN;
  point3d N, edgeL, elementL;
  
  
  //Interfaces data________________________________________________________________________________
  for(UInt i=1; i <= interfaces.size(); ++i)
  {
    //Interface data
    edgeId = interfaceToEdge(i);
    
    id1 = interfaces(i).getCid(1);
    id2 = interfaces(i).getCid(2);
    id3 = interfaces(i).getCid(3);
    
    P1 = dualPoints(id1);
    P2 = dualPoints(id2);
    P3 = dualPoints(id3);
    
    N = (P2-P1) ^ (P3-P1);
    
    //Edge data
    id1 = grid3d->getEdgeL(edgeId).getCid(1);
    id2 = grid3d->getEdgeL(edgeId).getCid(2);
    
    P1 = grid3d->getNodeL(id1);
    P2 = grid3d->getNodeL(id2);
    
    edgeL = P2 - P1;
    edgeL = edgeL / edgeL.norm2();
    
    //Creazione normale
    if((N * edgeL) >= 0.0)
    {
      edgeN(edgeId) = edgeN(edgeId) + N;
    }
    else
    {
      edgeN(edgeId) = edgeN(edgeId) - N;
    }
  }
  
  //Normal normalization___________________________________________________________________________
  assert(grid3d->getNumEdges() == edgeN.size());
  for(UInt i=1; i <= grid3d->getNumEdges(); ++i)
  {
    edgeSurf(i) = 0.5 * edgeN(i).norm2();
    edgeN(i)    = edgeN(i) / edgeN(i).norm2();
  }
  
  //Compute sub cell volumes_______________________________________________________________________
  for(UInt i=1; i <= subCells.size(); ++i)
  {
    cell = subCellToNode(i);
    
    id1 = subCells(i).getCid(1);
    id2 = subCells(i).getCid(2);
    id3 = subCells(i).getCid(3);
    id4 = subCells(i).getCid(4);
    
    P1 = dualPoints(id1);
    P2 = dualPoints(id2);
    P3 = dualPoints(id3);
    P4 = dualPoints(id4);
    
    PM = (P2-P1)^(P3-P1);
    PN = (P4-P1)* (1.0/6.0);
    cellVolume(cell) += fabs( point3d::dot(PM,PN) );
  }
}

template<typename ELMAP, typename NODEMAP> 
void
dual3d<linearTetra,ELMAP,NODEMAP>::
buildFiniteVolumeDataOnly()
{
  //Assert_________________________________________________________________________________________
  assert(commDevLoaded);
  assert(grid3dLoaded);
  assert(connectGrid3dLoaded);
  
  dualCreated         = false;
  connectCreated      = true;
  finiteVolumeCreated = false;
  
  
  //Clearing_______________________________________________________________________________________
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
  
  subCellToElement.clear();
  subCellToNode.clear();  
  interfaceToEdge.clear();
  
  
  //Num local and global elements__________________________________________________________________
  UInt numLocalNodes    = grid3d->getNumNodes();
  UInt numLocalElements = grid3d->getNumElements();
  UInt numLocalFaces    = grid3d->getNumFaces();
  UInt numLocalEdges    = grid3d->getNumEdges();
  
  pMapGlobalManip<NODEMAP> vrMapManip(commDev);
  pMapGlobalManip<ELMAP>   elMapManip(commDev);
  pMapGlobalManip<ELMAP>   faMapManip(commDev);
  pMapGlobalManip<ELMAP>   edMapManip(commDev);
  
  pMap<NODEMAP> vrMap = grid3d->getNodes().getMapRef();
  pMap<ELMAP>   elMap = grid3d->getElements().getRowMap();
  pMap<ELMAP>   faMap = grid3d->getFaces().getRowMap();
  pMap<ELMAP>   edMap = grid3d->getEdges().getRowMap();
  
  UInt numGlobalNodes    = vrMapManip.sizeG(vrMap);
  UInt numGlobalElements = elMapManip.sizeG(elMap);
  UInt numGlobalFaces    = faMapManip.sizeG(faMap);
  
  
  //Reserving space________________________________________________________________________________
  dualPoints.reserve(numLocalNodes + numLocalElements + numLocalFaces + numLocalEdges);
  subCells.reserve(24 * numLocalElements);
  subCellToElement.reserve(24 * numLocalElements);
  subCellToNode.reserve(24 * numLocalElements);
  interfaces.reserve(12 * numLocalElements);
  interfaceToEdge.reserve(12 * numLocalElements);
  
  
  //Allocations____________________________________________________________________________________
  UInt pid = commDev->rank();
  
  UInt    id1, id2, id3, id4;
  point3d P1, P2, P3, P4, PM;
  
  GEOELEMENT3D subCell(true);
  sVect<UInt>  f(4), e(6), p(4);
  UInt ba;
  sVect<point3d> F(4), E(6), P(4);
  point3d BA;
  
  Real test;
  UInt no1, no2, edgeId, elGid, locStray, geoId;
  
  GEOELEMENT2D interface(true);
  UInt numSubCells = 1;
  UInt numInterfaces = 1;
  UInt numPoints;
  
  ELMAP   elMapItem;
  NODEMAP nodeMapItem;
  
  
  //Nodes construction_____________________________________________________________________________
  dualPoints = grid3d->getNodes();   //Orginal nodes
  numPoints  = dualPoints.size() + 1;
  
  for(UInt i=1; i <= numLocalElements; ++i)    //Nodi baricentrici
  {
    id1 = grid3d->getElementL(i).getCid(1);
    id2 = grid3d->getElementL(i).getCid(2);
    id3 = grid3d->getElementL(i).getCid(3);
    id4 = grid3d->getElementL(i).getCid(4);
    
    P1  = grid3d->getNodeL(id1);
    P2  = grid3d->getNodeL(id2);
    P3  = grid3d->getNodeL(id3);
    P4  = grid3d->getNodeL(id4);
    
    PM  = (P1 + P2 + P3 + P4) / 4.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes +  grid3d->getElements().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  for(UInt i=1; i <= numLocalFaces; ++i) 	//Nodi delle facce
  {
    id1 = grid3d->getFaceL(i).getCid(1);
    id2 = grid3d->getFaceL(i).getCid(2);
    id3 = grid3d->getFaceL(i).getCid(3);
    
    P1  = grid3d->getNodeL(id1);
    P2  = grid3d->getNodeL(id2);
    P3  = grid3d->getNodeL(id3);
    
    PM  = (P1 + P2 + P3) / 3.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes + numGlobalElements +  grid3d->getFaces().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  for(UInt i=1; i <= numLocalEdges; ++i)       //Nodi degli spigoli 
  {
    id1 = grid3d->getEdgeL(i).getCid(1);
    id2 = grid3d->getEdgeL(i).getCid(2);
    
    P1  = grid3d->getNodeL(id1);
    P2  = grid3d->getNodeL(id2);
    
    PM = (P1 + P2) / 2.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes + numGlobalElements + numGlobalFaces +  grid3d->getEdges().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }  
  
  
  //Subcells and interfaces construction___________________________________________________________
  for(UInt i=1; i <= numLocalElements; ++i)
  {
    elGid = grid3d->getElements().getRowMapL(i).getGid();
    
    //Tetrahedra nodes
    for(UInt j=1; j <= 4; ++j) 
    {
      p(j) = grid3d->getElementL(i).getCid(j);
      P(j) = dualPoints(p(j));
    }
    
    //Tetrahedra barycenter
    ba = getLocalNodeId(2,i);
    BA = dualPoints(ba);
    
    //Nodes on the tetra faces
    for(UInt j=1; j <= 4; ++j)
    {
      f(j) = connectGrid3d->getElementToFace(i,j);
      f(j) = getLocalNodeId(3,f(j));
      F(j) = dualPoints(f(j));
    }
    
    //Nodes on the tetra edges
    for(UInt j=1; j <= 6; ++j)
    {
      e(j) = connectGrid3d->getElementToEdge(i,j);
      e(j) = getLocalNodeId(4,e(j));
      E(j) = dualPoints(e(j));
    }
    
    //Subcells creation
    locStray = 1;
    elGid    = grid3d->getElements().getRowMapL(i).getGid();
    geoId    = grid3d->getElements().getItemL(i).getGeoId();
    
    for(UInt j=1; j <= 4; ++j) //Face cycle
    {
      for(UInt k=1; k <= 3; ++k) //Edges of the face
      {
	edgeId = refShape3d.faceToEdge(j,k);
	no1    = refShape3d.edgeToPoint(edgeId,1); //first node on the edge
	no2    = refShape3d.edgeToPoint(edgeId,2); //second node on the edge
	
	//First subcell
	test = point3d::dot((E(edgeId) - BA) ^ (F(j) - BA), BA - P(no1));
	
	if(test > 0.0)
	{
	  subCell.setCid(1,ba);
	  subCell.setCid(2,e(edgeId));
	  subCell.setCid(3,f(j));
	  subCell.setCid(4,p(no1));
	}
	else
	{
	  subCell.setCid(1,ba);
	  subCell.setCid(2,f(j));
	  subCell.setCid(3,e(edgeId));
	  subCell.setCid(4,p(no1));
	}
	  
	elMapItem.setPid(pid);
	elMapItem.setLid(numSubCells);
	elMapItem.setGid(24 * (elGid-1) + locStray);
	
	subCell.setGeoId(geoId);
	  
	locStray++;
	numSubCells++;
	  
	subCells.push_back(elMapItem,subCell);
	subCellToElement.push_back(elMapItem,i);
	subCellToNode.push_back(elMapItem,p(no1));
	
	
	//Seconda sottocella 
	test = point3d::dot((F(j) - BA) ^ (E(edgeId) - BA), BA - P(no2));
	
	if(test > 0.0)
	{
	  subCell.setCid(1,ba);
	  subCell.setCid(2,f(j));
	  subCell.setCid(3,e(edgeId));
	  subCell.setCid(4,p(no2));
	}
	else
	{
	  subCell.setCid(1,ba);
	  subCell.setCid(2,e(edgeId));
	  subCell.setCid(3,f(j));
	  subCell.setCid(4,p(no2));
	}
	
	elMapItem.setPid(pid);
	elMapItem.setLid(numSubCells);
	elMapItem.setGid(24 * (elGid-1) + locStray);
	
	subCell.setGeoId(geoId);
	  
	locStray++;
	numSubCells++;
	  
	subCells.push_back(elMapItem,subCell);
	subCellToElement.push_back(elMapItem,i);
	subCellToNode.push_back(elMapItem,p(no2));
      }
    }
    
    
    //Interfaces creation
    locStray = 1;
    
    for(UInt j=1; j <= 4; ++j) //faces
    {
      for(UInt k=1; k <= 3; ++k) //edges on the face
      {
	edgeId = refShape3d.faceToEdge(j,k);
	
	interface.setCid(1,e(edgeId));
	interface.setCid(2,f(j));
	interface.setCid(3,ba);
	interface.setGeoId(geoId);
	
	elMapItem.setPid(pid);
	elMapItem.setLid(numInterfaces);
	elMapItem.setGid(12 * (elGid-1) + locStray);
	
	numInterfaces++;
	locStray++;
	
	interfaces.push_back(elMapItem,interface);
	interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
      }
    }
    
    assert(locStray == 13);
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
  
  
  //Building finite volume data and clear the remaining dual data____________________________________________________ 
  buildFiniteVolumeData();
  
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
}

template<typename ELMAP, typename NODEMAP> 
void
dual3d<linearTetra,ELMAP,NODEMAP>::
clearDual()
{
  dualCreated = false;
  
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
}

template<typename ELMAP, typename NODEMAP> 
void
dual3d<linearTetra,ELMAP,NODEMAP>::
clearConnect()
{
  connectCreated = false;
  
  subCellToElement.clear();
  subCellToNode.clear();
  interfaceToEdge.clear();
}

template<typename ELMAP, typename NODEMAP> 
void
dual3d<linearTetra,ELMAP,NODEMAP>::
clearFiniteVolume()
{
  finiteVolumeCreated = false;
  
  cellVolume.clear();
  edgeN.clear();
  edgeSurf.clear();
}

template<typename ELMAP, typename NODEMAP> 
void
dual3d<linearTetra,ELMAP,NODEMAP>::
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
    
  cellVolume.clear();
  edgeN.clear();
  edgeSurf.clear();
}



//_______________________________________________________________________________________________________
// LINEAR HEXA SPECILIZATION
//-------------------------------------------------------------------------------------------------------
/*! Support class for the implementation of vertex-centered finite volume schemes. Every hexahedra is decomposed in sub-cells and then aggregated to the mesh vertices.
Also the interfaces between the cells are created.

<b>Data</b>
<ol>
<li> dualPoints = the dual points: they are the union of the original grid vertices and of the barycentric nodes of the elements, of the faces and of the edges;
<li> dualCells = the grid subcells. Every hexa is divided into 8 subcelss. Every subcell is a regular hexa with a right-hand node numbering.
Every subcell is linked to a vertex of the original hexa, therefore each of the original 8 vertices of the hexa has one subcell.
Every subcell is built joining a vertex of the original hexa with three edges, three faces point and the hexa barycenter. The first of the nodes list is original point,
the rest of the nodes are listed so that the subcell-hexa fulfill the the right hand convenction. In all the cases the pattern of the nodes is:
originalNode - edgeNode - faceNode - edgeNode - edgeNode - faceNode - barycenter - faceNode;
<li> interfaces = the interfaces between the subcells. The finite volume method require the definition of interaces between cells. In our case the cell is the
union of the subcells sharing a given vertex of the original mesh. In the hexa case the interfaces are made of quads. Denoting with p1, p2, p3, p4 the nodes
of the interface the node p1 is on ad edge, the node p2 rests on a face and the node p3 is the barycenter of the element and the last is also an edge node.
</ol>

<b>Connectivity data</b>
<ol>
<li> subCellToElement: the connectivity of the subCells to the father element;
<li> subCellToNode: the connetivity of the subCell to a node;
<li> interfaceToEdge: the connectivity of the interface to the edge;
</ol>
  
<b>Finite volume data</b>
<ol>
<li> cellVolume: the volume of the cell;
<li> edgeN: the mean normal vector to the interface. The normal is obtained as a weighted normal vector of the nornals of the quads that compose the interface.
Every interface is associated to an edge.
<li> edgeSurf: the surface of the interface
</ol>
*/
template<typename ELMAP, typename NODEMAP>
class dual3d<linearHexa,ELMAP,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef linearHexa              GEOSHAPE3D;
    typedef linearHexa::GEOBSHAPE   GEOSHAPE2D;
    typedef GEOSHAPE2D::GEOBSHAPE   GEOSHAPE1D;
    typedef geoElement<GEOSHAPE3D>  GEOELEMENT3D;
    typedef geoElement<GEOSHAPE2D>  GEOELEMENT2D;
    typedef geoElement<GEOSHAPE1D>  GEOELEMENT1D;
    
    typedef mesh3d<GEOSHAPE3D,ELMAP,NODEMAP>    MESH3D;
    typedef connect3d<GEOSHAPE3D,ELMAP,NODEMAP> CONNECT3D;
    
    typedef pVect<point3d,NODEMAP>              NODESVECT;
    typedef pGraph<GEOELEMENT3D,ELMAP,NODEMAP>  GRAPH3D;
    typedef pGraph<GEOELEMENT2D,ELMAP,NODEMAP>  GRAPH2D;
    
    typedef pVect<Real,ELMAP>       CELLVOLUMESVECT;
    typedef pVect<point3d,ELMAP>    EDGENSVECT;
    typedef pVect<Real,ELMAP>       EDGESSURFVECT;
    
    typedef pVect<UInt,ELMAP>       SUBCELL_TO_ELEMENT;
    typedef pVect<UInt,ELMAP>       SUBCELL_TO_NODE;
    typedef pVect<UInt,ELMAP>       INTERFACE_TO_EDGE;
    //@}
      
    /*! @name Internal flags */ //@{ 
  public:
    bool commDevLoaded;
    bool grid3dLoaded;
    bool connectGrid3dLoaded;
    bool dualCreated;
    bool connectCreated;
    bool finiteVolumeCreated;
    //@}
      
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<communicator> commDev;
    Teuchos::RCP<MESH3D>       grid3d;
    Teuchos::RCP<CONNECT3D>    connectGrid3d;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    geoMapInterface<GEOSHAPE3D>  refShape3d;
    NODESVECT dualPoints;
    GRAPH3D   subCells;
    GRAPH2D   interfaces;
    //@}
    
     /*! @name Finite Volumes data */ //@{ 
  public:
    CELLVOLUMESVECT cellVolume;
    EDGENSVECT      edgeN;
    EDGESSURFVECT   edgeSurf;
    //@}
    
    /*! @name Internal connectivity data */ //@{ 
  public:
    SUBCELL_TO_ELEMENT subCellToElement;
    SUBCELL_TO_NODE    subCellToNode;  
    INTERFACE_TO_EDGE  interfaceToEdge;
    //@}
    
    /*! @name Constructors */ //@{ 
  public:
    dual3d();
    dual3d(const Teuchos::RCP<communicator> & CommDev);
    dual3d(communicator & CommDev);
    void setMesh3d(const Teuchos::RCP<MESH3D> & Grid3d);
    void setMesh3d(const MESH3D & Grid3d);
    void setConnect3d(const Teuchos::RCP<CONNECT3D> & Connect3d);
    void setConnect3d(const CONNECT3D & Connect3d);
    //@}
    
    /*! @name Funzioni interne e buildup */ //@{ 
  public:    
    /*! The dual nodes is the union of the original nodes + the elements barycentric nodes + the faces barycentric nodes + the edges ones.
    Given the id of a geoItem (such as an element, face, etc) this function produces the local id of the dual node associated.
    \param flag = 1 the original vertices, = 2 the element barycenter, = 3 the face one, = 4 the edge one
    \param lid the local numbering of the geoItem (such as an element, face, etc) */
    UInt getLocalNodeId(const UInt & flag, const UInt & lid) const;
      
    /*! Build dual mesh. Explicitly builds the dual mesh */
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
dual3d<linearHexa,ELMAP,NODEMAP>::
dual3d()
{
  commDevLoaded       = false;
  grid3dLoaded        = false;
  connectGrid3dLoaded = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
  dualCreated         = false;
}

template<typename ELMAP, typename NODEMAP>
dual3d<linearHexa,ELMAP,NODEMAP>::
dual3d(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded       = true;
  grid3dLoaded        = false;
  connectGrid3dLoaded = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
  dualCreated         = false;
  
  commDev = CommDev;
}

template<typename ELMAP, typename NODEMAP>
dual3d<linearHexa,ELMAP,NODEMAP>::
dual3d(communicator & CommDev)
{
  commDevLoaded       = true;
  grid3dLoaded        = false;
  connectGrid3dLoaded = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
  dualCreated         = false;
  
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename ELMAP, typename NODEMAP>
void
dual3d<linearHexa,ELMAP,NODEMAP>::
setMesh3d(const Teuchos::RCP<MESH3D> & Grid3d)
{
  grid3dLoaded = true;
  grid3d       = Grid3d;
}

template<typename ELMAP, typename NODEMAP>
void
dual3d<linearHexa,ELMAP,NODEMAP>::
setMesh3d(const MESH3D & Grid3d)
{
  grid3dLoaded = true;
  grid3d       = Teuchos::rcp(new MESH3D(Grid3d));
}

template<typename ELMAP, typename NODEMAP>
void
dual3d<linearHexa,ELMAP,NODEMAP>::
setConnect3d(const Teuchos::RCP<CONNECT3D> & Connect3d)
{
  connectGrid3dLoaded = true;
  connectGrid3d       = Connect3d;
}

template<typename ELMAP, typename NODEMAP>
void
dual3d<linearHexa,ELMAP,NODEMAP>::
setConnect3d(const CONNECT3D & Connect3d)
{
  connectGrid3dLoaded = true;
  connectGrid3d       = Teuchos::rcp(new CONNECT3D(Connect3d));
}

template<typename ELMAP, typename NODEMAP> 
UInt
dual3d<linearHexa,ELMAP,NODEMAP>::
getLocalNodeId(const UInt & flag, const UInt & lid) const
{ 
  return( UInt(flag > 1) * grid3d->getNumNodes()
        + UInt(flag > 2) * grid3d->getNumElements()
        + UInt(flag > 3) * grid3d->getNumFaces()
        + lid );
}

template<typename ELMAP, typename NODEMAP> 
void
dual3d<linearHexa,ELMAP,NODEMAP>::
buildDualMesh()
{
  //Assert_________________________________________________________________________________________
  assert(commDevLoaded);
  assert(grid3dLoaded);
  assert(connectGrid3dLoaded);
  
  connectCreated = true;
  
  
  //Clearing_______________________________________________________________________________________
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
  
  subCellToElement.clear();
  subCellToNode.clear();  
  interfaceToEdge.clear();
  
  
  //Num local and global elements__________________________________________________________________
  UInt numLocalNodes    = grid3d->getNumNodes();
  UInt numLocalElements = grid3d->getNumElements();
  UInt numLocalFaces    = grid3d->getNumFaces();
  UInt numLocalEdges    = grid3d->getNumEdges();
  
  pMapGlobalManip<NODEMAP> vrMapManip(commDev);
  pMapGlobalManip<ELMAP>   elMapManip(commDev);
  pMapGlobalManip<ELMAP>   faMapManip(commDev);
  pMapGlobalManip<ELMAP>   edMapManip(commDev);
  
  pMap<NODEMAP> vrMap = grid3d->getNodes().getMapRef();
  pMap<ELMAP>   elMap = grid3d->getElements().getRowMap();
  pMap<ELMAP>   faMap = grid3d->getFaces().getRowMap();
  pMap<ELMAP>   edMap = grid3d->getEdges().getRowMap();
  
  UInt numGlobalNodes    = vrMapManip.sizeG(vrMap);
  UInt numGlobalElements = elMapManip.sizeG(elMap);
  UInt numGlobalFaces    = faMapManip.sizeG(faMap);
  
  
  //Reserving space________________________________________________________________________________
  dualPoints.reserve(numLocalNodes + numLocalElements + numLocalFaces + numLocalEdges);
  subCells.reserve(8 * numLocalElements);
  subCellToElement.reserve(8 * numLocalElements);
  subCellToNode.reserve(8 * numLocalElements);
  interfaces.reserve(12 * numLocalElements);
  interfaceToEdge.reserve(12 * numLocalElements);
  
  
  //Allocations____________________________________________________________________________________
  UInt pid = commDev->rank();
  
  UInt    id1, id2, id3, id4, id5, id6, id7, id8;
  point3d P1, P2, P3, P4, P5, P6, P7, P8, PM;
  
  GEOELEMENT3D subCell(true);
  sVect<UInt>  f(6), e(12), p(8);
  UInt ba;
  point3d BA;
  
  UInt edgeId, elGid, locStray, geoId;
  
  GEOELEMENT2D interface(true);
  UInt numSubCells = 1;
  UInt numInterfaces = 1;
  UInt numPoints;
  
  ELMAP   elMapItem;
  NODEMAP nodeMapItem;
  
  
  //Nodes construction_____________________________________________________________________________
  dualPoints = grid3d->getNodes();   //Orginal nodes
  numPoints  = dualPoints.size() + 1;
  
  for(UInt i=1; i <= numLocalElements; ++i)    //Nodi baricentrici
  {
    id1 = grid3d->getElementL(i).getCid(1);
    id2 = grid3d->getElementL(i).getCid(2);
    id3 = grid3d->getElementL(i).getCid(3);
    id4 = grid3d->getElementL(i).getCid(4);
    
    id5 = grid3d->getElementL(i).getCid(5);
    id6 = grid3d->getElementL(i).getCid(6);
    id7 = grid3d->getElementL(i).getCid(7);
    id8 = grid3d->getElementL(i).getCid(8);
    
    P1  = grid3d->getNodeL(id1);
    P2  = grid3d->getNodeL(id2);
    P3  = grid3d->getNodeL(id3);
    P4  = grid3d->getNodeL(id4);
    
    P5  = grid3d->getNodeL(id5);
    P6  = grid3d->getNodeL(id6);
    P7  = grid3d->getNodeL(id7);
    P8  = grid3d->getNodeL(id8);
    
    PM  = (P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8) / 8.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes +  grid3d->getElements().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  for(UInt i=1; i <= numLocalFaces; ++i) 	//Nodi delle facce
  {
    id1 = grid3d->getFaceL(i).getCid(1);
    id2 = grid3d->getFaceL(i).getCid(2);
    id3 = grid3d->getFaceL(i).getCid(3);
    id4 = grid3d->getFaceL(i).getCid(4);
    
    P1  = grid3d->getNodeL(id1);
    P2  = grid3d->getNodeL(id2);
    P3  = grid3d->getNodeL(id3);
    P4  = grid3d->getNodeL(id4);
    
    PM  = (P1 + P2 + P3 + P4) / 4.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes + numGlobalElements +  grid3d->getFaces().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  for(UInt i=1; i <= numLocalEdges; ++i)       //Nodi degli spigoli 
  {
    id1 = grid3d->getEdgeL(i).getCid(1);
    id2 = grid3d->getEdgeL(i).getCid(2);
    
    P1  = grid3d->getNodeL(id1);
    P2  = grid3d->getNodeL(id2);
    
    PM = (P1 + P2) / 2.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes + numGlobalElements + numGlobalFaces +  grid3d->getEdges().getRowMapL(i).getGid());
    
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
  
  
  for(UInt i=1; i <= numLocalElements; ++i)
  {
    elMapItem = grid3d->getElements().getRowMapL(i);
    elGid     = grid3d->getElements().getRowMapL(i).getGid();
    
    //Hexa nodes
    for(UInt j=1; j <= 8; ++j) 
    {
      p(j) = grid3d->getElementL(i).getCid(j);
    }
    
    //Hexa barycenter
    ba = getLocalNodeId(2,i);
    BA = dualPoints(ba);
    
    //Nodes on the Hexa faces
    for(UInt j=1; j <= 6; ++j)
    {
      f(j) = connectGrid3d->getElementToFace(i,j);
      f(j) = getLocalNodeId(3,f(j));
    }
    
    //Nodes on the Hexa edges
    for(UInt j=1; j <= 12; ++j)
    {
      e(j) = connectGrid3d->getElementToEdge(i,j);
      e(j) = getLocalNodeId(4,e(j));
    }
    
    
    //Subcells creation
    locStray = 1;
    elGid    = grid3d->getElements().getRowMapL(i).getGid();
    geoId    = grid3d->getElements().getItemL(i).getGeoId();
    
    elMapItem.setPid(pid);
    subCell.setGeoId(geoId);
    
    //Subcell 1--------------------------------------------
    subCell.setCid(1,p(1));
    subCell.setCid(2,e(1));
    subCell.setCid(3,f(1));
    subCell.setCid(4,e(4));
    subCell.setCid(5,e(5));
    subCell.setCid(6,f(3));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(2));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(1));
    
    //Subcell 2--------------------------------------------
    subCell.setCid(1,p(2));
    subCell.setCid(2,e(2));
    subCell.setCid(3,f(1));
    subCell.setCid(4,e(1));
    subCell.setCid(5,e(6));
    subCell.setCid(6,f(4));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(3));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(2));
    
    //Subcell 3--------------------------------------------
    subCell.setCid(1,p(3));
    subCell.setCid(2,e(3));
    subCell.setCid(3,f(1));
    subCell.setCid(4,e(2));
    subCell.setCid(5,e(7));
    subCell.setCid(6,f(5));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(4));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(3));
    
    //Subcell 4--------------------------------------------
    subCell.setCid(1,p(4));
    subCell.setCid(2,e(4));
    subCell.setCid(3,f(1));
    subCell.setCid(4,e(3));
    subCell.setCid(5,e(8));
    subCell.setCid(6,f(2));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(5));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(4));
    
    //Subcell 5--------------------------------------------
    subCell.setCid(1,p(5));
    subCell.setCid(2,e(12));
    subCell.setCid(3,f(6));
    subCell.setCid(4,e(9));
    subCell.setCid(5,e(5));
    subCell.setCid(6,f(2));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(3));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(5));
    
    //Subcell 6--------------------------------------------
    subCell.setCid(1,p(6));
    subCell.setCid(2,e(9));
    subCell.setCid(3,f(6));
    subCell.setCid(4,e(10));
    subCell.setCid(5,e(6));
    subCell.setCid(6,f(3));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(4));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(6));
    
    //Subcell 7--------------------------------------------
    subCell.setCid(1,p(7));
    subCell.setCid(2,e(10));
    subCell.setCid(3,f(6));
    subCell.setCid(4,e(11));
    subCell.setCid(5,e(7));
    subCell.setCid(6,f(4));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(5));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(7));
    
    //Subcell 8--------------------------------------------
    subCell.setCid(1,p(8));
    subCell.setCid(2,e(11));
    subCell.setCid(3,f(6));
    subCell.setCid(4,e(12));
    subCell.setCid(5,e(8));
    subCell.setCid(6,f(5));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(2));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(8));
    
    
    //Interfaces creation__________________________________________________________________________
    locStray = 1;
    
    //Interface 1------------------------------------------
    edgeId = 1;
    
    interface.setCid(1,e(1));
    interface.setCid(2,f(1));
    interface.setCid(3,ba);
    interface.setCid(4,f(3));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 2------------------------------------------
    edgeId = 2;
    
    interface.setCid(1,e(2));
    interface.setCid(2,f(1));
    interface.setCid(3,ba);
    interface.setCid(4,f(4));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 3------------------------------------------
    edgeId = 3;
    
    interface.setCid(1,e(3));
    interface.setCid(2,f(1));
    interface.setCid(3,ba);
    interface.setCid(4,f(5));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 4------------------------------------------
    edgeId = 4;
    
    interface.setCid(1,e(4));
    interface.setCid(2,f(1));
    interface.setCid(3,ba);
    interface.setCid(4,f(2));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 5------------------------------------------
    edgeId = 5;
    
    interface.setCid(1,e(5));
    interface.setCid(2,f(3));
    interface.setCid(3,ba);
    interface.setCid(4,f(2));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 6------------------------------------------
    edgeId = 6;
    
    interface.setCid(1,e(6));
    interface.setCid(2,f(4));
    interface.setCid(3,ba);
    interface.setCid(4,f(3));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 7------------------------------------------
    edgeId = 7;
    
    interface.setCid(1,e(7));
    interface.setCid(2,f(5));
    interface.setCid(3,ba);
    interface.setCid(4,f(4));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 8------------------------------------------
    edgeId = 8;
    
    interface.setCid(1,e(8));
    interface.setCid(2,f(2));
    interface.setCid(3,ba);
    interface.setCid(4,f(5));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 9------------------------------------------
    edgeId = 9;
    
    interface.setCid(1,e(9));
    interface.setCid(2,f(3));
    interface.setCid(3,ba);
    interface.setCid(4,f(6));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 10-----------------------------------------
    edgeId = 10;
    
    interface.setCid(1,e(10));
    interface.setCid(2,f(4));
    interface.setCid(3,ba);
    interface.setCid(4,f(6));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 11-----------------------------------------
    edgeId = 11;
    
    interface.setCid(1,e(11));
    interface.setCid(2,f(5));
    interface.setCid(3,ba);
    interface.setCid(4,f(6));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 12-----------------------------------------
    edgeId = 12;
    
    interface.setCid(1,e(12));
    interface.setCid(2,f(2));
    interface.setCid(3,ba);
    interface.setCid(4,f(6));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
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
  pGraphGlobalManip<GEOELEMENT3D,ELMAP,NODEMAP> subCellsCheker(commDev);
  assert(subCellsCheker.check(subCells));
  
  
  //MapFixing for interfaces
  interfaces.setColMap(dualPoints.getMapRef());
  interfaces.updateRowFinder();
  interfaces.updateColFinder();
  
  interfaceToEdge.setMap(interfaces.getRowMap());
  interfaceToEdge.updateFinder();
  
  //Checking interfaces map
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> interfacesCheker(commDev);
  assert(interfacesCheker.check(interfaces));
}

template<typename ELMAP, typename NODEMAP> 
void
dual3d<linearHexa,ELMAP,NODEMAP>::
buildFiniteVolumeData()
{
  //Assert_________________________________________________________________________________________
  assert(commDevLoaded);
  assert(grid3dLoaded);
  assert(connectGrid3dLoaded);
  
  finiteVolumeCreated = true;
  
  //Clearing_______________________________________________________________________________________
  cellVolume.clear();
  edgeN.clear();
  edgeSurf.clear();
  
  //Resizing_______________________________________________________________________________________
  edgeN.resize(grid3d->getNumEdges());
  edgeSurf.resize(grid3d->getNumEdges());
  cellVolume.resize(grid3d->getNumNodes());
  
  //Load maps______________________________________________________________________________________
  edgeN.setMap(grid3d->getEdges().getRowMap());
  edgeSurf.setMap(grid3d->getEdges().getRowMap());
  cellVolume.setMap(grid3d->getNodes().getMapRef());
  
  //Allocate_______________________________________________________________________________________
  UInt id1, id2, id3, id4, id5, id6, id7, id8;
  UInt edgeId, cell;
  Real surface;
  point3d P1, P2, N, edgeL, elementL;
  sVect<point3d> points2d(4);
  sVect<point3d> points3d(8);
  geoMapSupport2d<GEOSHAPE2D> geoSupport2d;
  geoMapSupport3d<GEOSHAPE3D> geoSupport3d;
  
  
  //Interfaces data________________________________________________________________________________
  for(UInt i=1; i <= interfaces.size(); ++i)
  {
    //Interface data
    edgeId = interfaceToEdge(i);
    
    id1 = interfaces(i).getCid(1);
    id2 = interfaces(i).getCid(2);
    id3 = interfaces(i).getCid(3);
    id4 = interfaces(i).getCid(4);
    
    points2d(1) = dualPoints(id1);
    points2d(2) = dualPoints(id2);
    points2d(3) = dualPoints(id3);
    points2d(4) = dualPoints(id4);
     
    geoSupport2d.setPoints(points2d);
    geoSupport2d.surfaceNormal<STANDARD,3>(N,surface);
    
    
    id1 = grid3d->getEdgeL(edgeId).getCid(1);
    id2 = grid3d->getEdgeL(edgeId).getCid(2);
    
    P1 = grid3d->getNodeL(id1);
    P2 = grid3d->getNodeL(id2);
    
    edgeL = P2 - P1;
    
    
    //Creazione normale
    if((N * edgeL) >= 0.0)
    {
      edgeN(edgeId)    += N * surface;
      edgeSurf(edgeId) += surface;
    }
    else
    {
      edgeN(edgeId)    -= N * surface;
      edgeSurf(edgeId) += surface;
    }
  }
  
  
  //Normal normalization___________________________________________________________________________
  assert(grid3d->getNumEdges() == edgeN.size());
  for(UInt i=1; i <= grid3d->getNumEdges(); ++i)
  {
    edgeN(i) = edgeN(i) / edgeSurf(i);
  }
  
  
  //Calcolo dei volumi di cella____________________________________________________________________
  for(UInt i=1; i <= subCells.size(); ++i)
  {
    cell = subCellToNode(i);
    
    id1 = subCells(i).getCid(1);
    id2 = subCells(i).getCid(2);
    id3 = subCells(i).getCid(3);
    id4 = subCells(i).getCid(4);
    id5 = subCells(i).getCid(5);
    id6 = subCells(i).getCid(6);
    id7 = subCells(i).getCid(7);
    id8 = subCells(i).getCid(8);
    
    points3d(1) = dualPoints(id1);
    points3d(2) = dualPoints(id2);
    points3d(3) = dualPoints(id3);
    points3d(4) = dualPoints(id4);
    points3d(5) = dualPoints(id5);
    points3d(6) = dualPoints(id6);
    points3d(7) = dualPoints(id7);
    points3d(8) = dualPoints(id8);
    
    geoSupport3d.setPoints(points3d);    
    cellVolume(cell) += geoSupport3d.volume<STANDARD,3>();
  }
  
  
  //Clearing of all the remaining 
  cellVolume.clear();
  edgeN.clear();
  edgeSurf.clear();
}

template<typename ELMAP, typename NODEMAP> 
void
dual3d<linearHexa,ELMAP,NODEMAP>::
buildFiniteVolumeDataOnly()
{
  //Assert_________________________________________________________________________________________
  assert(commDevLoaded);
  assert(grid3dLoaded);
  assert(connectGrid3dLoaded);
  
  dualCreated         = false;
  connectCreated      = true;
  finiteVolumeCreated = true;
  
  
  //Clearing_______________________________________________________________________________________
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
  
  subCellToElement.clear();
  subCellToNode.clear();  
  interfaceToEdge.clear();
  
  
  //Num local and global elements__________________________________________________________________
  UInt numLocalNodes    = grid3d->getNumNodes();
  UInt numLocalElements = grid3d->getNumElements();
  UInt numLocalFaces    = grid3d->getNumFaces();
  UInt numLocalEdges    = grid3d->getNumEdges();
  
  pMapGlobalManip<NODEMAP> vrMapManip(commDev);
  pMapGlobalManip<ELMAP>   elMapManip(commDev);
  pMapGlobalManip<ELMAP>   faMapManip(commDev);
  pMapGlobalManip<ELMAP>   edMapManip(commDev);
  
  pMap<NODEMAP> vrMap = grid3d->getNodes().getMapRef();
  pMap<ELMAP>   elMap = grid3d->getElements().getRowMap();
  pMap<ELMAP>   faMap = grid3d->getFaces().getRowMap();
  pMap<ELMAP>   edMap = grid3d->getEdges().getRowMap();
  
  UInt numGlobalNodes    = vrMapManip.sizeG(vrMap);
  UInt numGlobalElements = elMapManip.sizeG(elMap);
  UInt numGlobalFaces    = faMapManip.sizeG(faMap);
  
  
  //Reserving space________________________________________________________________________________
  dualPoints.reserve(numLocalNodes + numLocalElements + numLocalFaces + numLocalEdges);
  subCells.reserve(8 * numLocalElements);
  subCellToElement.reserve(8 * numLocalElements);
  subCellToNode.reserve(8 * numLocalElements);
  interfaces.reserve(12 * numLocalElements);
  interfaceToEdge.reserve(12 * numLocalElements);
  
  
  //Allocations____________________________________________________________________________________
  UInt pid = commDev->rank();
  
  UInt    id1, id2, id3, id4, id5, id6, id7, id8;
  point3d P1, P2, P3, P4, P5, P6, P7, P8, PM;
  
  GEOELEMENT3D subCell(true);
  sVect<UInt>  f(6), e(12), p(8);
  UInt ba;
  point3d BA;
  
  UInt edgeId, elGid, locStray, geoId;
  
  GEOELEMENT2D interface(true);
  UInt numSubCells = 1;
  UInt numInterfaces = 1;
  UInt numPoints;
  
  ELMAP   elMapItem;
  NODEMAP nodeMapItem;
  
  
  //Nodes construction_____________________________________________________________________________
  dualPoints = grid3d->getNodes();   //Orginal nodes
  numPoints  = dualPoints.size() + 1;
  
  for(UInt i=1; i <= numLocalElements; ++i)    //Nodi baricentrici
  {
    id1 = grid3d->getElementL(i).getCid(1);
    id2 = grid3d->getElementL(i).getCid(2);
    id3 = grid3d->getElementL(i).getCid(3);
    id4 = grid3d->getElementL(i).getCid(4);
    
    id5 = grid3d->getElementL(i).getCid(5);
    id6 = grid3d->getElementL(i).getCid(6);
    id7 = grid3d->getElementL(i).getCid(7);
    id8 = grid3d->getElementL(i).getCid(8);
    
    P1  = grid3d->getNodeL(id1);
    P2  = grid3d->getNodeL(id2);
    P3  = grid3d->getNodeL(id3);
    P4  = grid3d->getNodeL(id4);
    
    P5  = grid3d->getNodeL(id5);
    P6  = grid3d->getNodeL(id6);
    P7  = grid3d->getNodeL(id7);
    P8  = grid3d->getNodeL(id8);
    
    PM  = (P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8) / 8.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes +  grid3d->getElements().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  for(UInt i=1; i <= numLocalFaces; ++i) 	//Nodi delle facce
  {
    id1 = grid3d->getFaceL(i).getCid(1);
    id2 = grid3d->getFaceL(i).getCid(2);
    id3 = grid3d->getFaceL(i).getCid(3);
    id4 = grid3d->getFaceL(i).getCid(4);
    
    P1  = grid3d->getNodeL(id1);
    P2  = grid3d->getNodeL(id2);
    P3  = grid3d->getNodeL(id3);
    P4  = grid3d->getNodeL(id4);
    
    PM  = (P1 + P2 + P3 + P4) / 4.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes + numGlobalElements +  grid3d->getFaces().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  for(UInt i=1; i <= numLocalEdges; ++i)       //Nodi degli spigoli 
  {
    id1 = grid3d->getEdgeL(i).getCid(1);
    id2 = grid3d->getEdgeL(i).getCid(2);
    
    P1  = grid3d->getNodeL(id1);
    P2  = grid3d->getNodeL(id2);
    
    PM = (P1 + P2) / 2.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes + numGlobalElements + numGlobalFaces +  grid3d->getEdges().getRowMapL(i).getGid());
    
    dualPoints.push_back(nodeMapItem,PM);
    numPoints++;
  }
  
  
  for(UInt i=1; i <= numLocalElements; ++i)
  {
    elGid = grid3d->getElements().getRowMapL(i).getGid();
    
    //Hexa nodes
    for(UInt j=1; j <= 8; ++j) 
    {
      p(j) = grid3d->getElementL(i).getCid(j);
    }
    
    //Hexa barycenter
    ba = getLocalNodeId(2,i);
    BA = dualPoints(ba);
    
    //Nodes on the Hexa faces
    for(UInt j=1; j <= 6; ++j)
    {
      f(j) = connectGrid3d->getElementToFace(i,j);
      f(j) = getLocalNodeId(3,f(j));
    }
    
    //Nodes on the Hexa edges
    for(UInt j=1; j <= 12; ++j)
    {
      e(j) = connectGrid3d->getElementToEdge(i,j);
      e(j) = getLocalNodeId(4,e(j));
    }
    
    
    //Subcells creation
    locStray = 1;
    elGid    = grid3d->getElements().getRowMapL(i).getGid();
    geoId    = grid3d->getElements().getItemL(i).getGeoId();
    
    elMapItem.setPid(pid);
    subCell.setGeoId(geoId);
    
    //Subcell 1--------------------------------------------
    subCell.setCid(1,p(1));
    subCell.setCid(2,e(1));
    subCell.setCid(3,f(1));
    subCell.setCid(4,e(4));
    subCell.setCid(5,e(5));
    subCell.setCid(6,f(3));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(2));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(1));
    
    //Subcell 2--------------------------------------------
    subCell.setCid(1,p(2));
    subCell.setCid(2,e(2));
    subCell.setCid(3,f(1));
    subCell.setCid(4,e(1));
    subCell.setCid(5,e(6));
    subCell.setCid(6,f(4));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(3));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(2));
    
    //Subcell 3--------------------------------------------
    subCell.setCid(1,p(3));
    subCell.setCid(2,e(3));
    subCell.setCid(3,f(1));
    subCell.setCid(4,e(2));
    subCell.setCid(5,e(7));
    subCell.setCid(6,f(5));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(4));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(3));
    
    //Subcell 4--------------------------------------------
    subCell.setCid(1,p(4));
    subCell.setCid(2,e(4));
    subCell.setCid(3,f(1));
    subCell.setCid(4,e(3));
    subCell.setCid(5,e(8));
    subCell.setCid(6,f(2));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(5));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(4));
    
    //Subcell 5--------------------------------------------
    subCell.setCid(1,p(5));
    subCell.setCid(2,e(12));
    subCell.setCid(3,f(6));
    subCell.setCid(4,e(9));
    subCell.setCid(5,e(5));
    subCell.setCid(6,f(2));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(3));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(5));
    
    //Subcell 6--------------------------------------------
    subCell.setCid(1,p(6));
    subCell.setCid(2,e(9));
    subCell.setCid(3,f(6));
    subCell.setCid(4,e(10));
    subCell.setCid(5,e(6));
    subCell.setCid(6,f(3));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(4));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(6));
    
    //Subcell 7--------------------------------------------
    subCell.setCid(1,p(7));
    subCell.setCid(2,e(10));
    subCell.setCid(3,f(6));
    subCell.setCid(4,e(11));
    subCell.setCid(5,e(7));
    subCell.setCid(6,f(4));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(5));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(7));
    
    //Subcell 8--------------------------------------------
    subCell.setCid(1,p(8));
    subCell.setCid(2,e(11));
    subCell.setCid(3,f(6));
    subCell.setCid(4,e(12));
    subCell.setCid(5,e(8));
    subCell.setCid(6,f(5));
    subCell.setCid(7,ba);
    subCell.setCid(8,f(2));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(8 * (elGid-1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(8));
    
    
    //Interfaces creation__________________________________________________________________________
    locStray = 1;
    
    //Interface 1------------------------------------------
    edgeId = 1;
    
    interface.setCid(1,e(1));
    interface.setCid(2,f(1));
    interface.setCid(3,ba);
    interface.setCid(4,f(3));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 2------------------------------------------
    edgeId = 2;
    
    interface.setCid(1,e(2));
    interface.setCid(2,f(1));
    interface.setCid(3,ba);
    interface.setCid(4,f(4));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 3------------------------------------------
    edgeId = 3;
    
    interface.setCid(1,e(3));
    interface.setCid(2,f(1));
    interface.setCid(3,ba);
    interface.setCid(4,f(5));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 4------------------------------------------
    edgeId = 4;
    
    interface.setCid(1,e(4));
    interface.setCid(2,f(1));
    interface.setCid(3,ba);
    interface.setCid(4,f(2));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 5------------------------------------------
    edgeId = 5;
    
    interface.setCid(1,e(5));
    interface.setCid(2,f(3));
    interface.setCid(3,ba);
    interface.setCid(4,f(2));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 6------------------------------------------
    edgeId = 6;
    
    interface.setCid(1,e(6));
    interface.setCid(2,f(4));
    interface.setCid(3,ba);
    interface.setCid(4,f(3));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 7------------------------------------------
    edgeId = 7;
    
    interface.setCid(1,e(7));
    interface.setCid(2,f(5));
    interface.setCid(3,ba);
    interface.setCid(4,f(4));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 8------------------------------------------
    edgeId = 8;
    
    interface.setCid(1,e(8));
    interface.setCid(2,f(2));
    interface.setCid(3,ba);
    interface.setCid(4,f(5));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 9------------------------------------------
    edgeId = 9;
    
    interface.setCid(1,e(9));
    interface.setCid(2,f(3));
    interface.setCid(3,ba);
    interface.setCid(4,f(6));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 10-----------------------------------------
    edgeId = 10;
    
    interface.setCid(1,e(10));
    interface.setCid(2,f(4));
    interface.setCid(3,ba);
    interface.setCid(4,f(6));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 11-----------------------------------------
    edgeId = 11;
    
    interface.setCid(1,e(11));
    interface.setCid(2,f(5));
    interface.setCid(3,ba);
    interface.setCid(4,f(6));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
    
    //Interface 12-----------------------------------------
    edgeId = 12;
    
    interface.setCid(1,e(12));
    interface.setCid(2,f(2));
    interface.setCid(3,ba);
    interface.setCid(4,f(6));
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(12 * (elGid-1) + locStray);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToEdge.push_back(elMapItem, connectGrid3d->getElementToEdge(i,edgeId));
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
  
  
  //Finite volume data creation and clearing_______________________________________________________
  buildFiniteVolumeData();
  
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
}

template<typename ELMAP, typename NODEMAP>
void
dual3d<linearHexa,ELMAP,NODEMAP>::
clearDual()
{
  dualCreated = false;
  
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
}
    
template<typename ELMAP, typename NODEMAP>
void
dual3d<linearHexa,ELMAP,NODEMAP>::
clearConnect()
{
  connectCreated = false;
  
  subCellToElement.clear();
  subCellToNode.clear();
  interfaceToEdge.clear();
}

template<typename ELMAP, typename NODEMAP> 
void
dual3d<linearHexa,ELMAP,NODEMAP>::
clearFiniteVolume()
{
  finiteVolumeCreated = false;
  
  cellVolume.clear();
  edgeN.clear();
  edgeSurf.clear();
}

template<typename ELMAP, typename NODEMAP> 
void
dual3d<linearHexa,ELMAP,NODEMAP>::
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
    
  cellVolume.clear();
  edgeN.clear();
  edgeSurf.clear();
}

#endif
