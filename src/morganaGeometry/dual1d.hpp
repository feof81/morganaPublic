/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef DUAL1D_HPP
#define DUAL1D_HPP

#include "algorithm"
#include "traitsGeometry.hpp"

#include "mesh1d.hpp"
#include "mesh2d.hpp"
#include "connect1d.hpp"
#include "connect2d.hpp"
#include "geoMapSupport1d.hpp"
#include "geoMapSupport2d.hpp"
#include <boost/concept_check.hpp>

using namespace std;


//_______________________________________________________________________________________________________
// DUMMY IMPLEMENTATION
//-------------------------------------------------------------------------------------------------------

/*! Dual mesh 1d topology and geometry - dummy implementation */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class dual1d
{ 
};



//_______________________________________________________________________________________________________
// LINEAR LINE SPECILIZATION
//-------------------------------------------------------------------------------------------------------

/*! Dual mesh 1d topology and geometry - \c linearLine implementation */
template<typename ELMAP, typename NODEMAP>
class dual1d<linearLine,ELMAP,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef linearLine              GEOSHAPE1D;
    typedef GEOSHAPE1D::GEOBSHAPE   GEOSHAPE0D;
    typedef geoElement<GEOSHAPE0D>  GEOELEMENT0D;
    typedef geoElement<GEOSHAPE1D>  GEOELEMENT1D;
    
    typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>    MESH1D;
    typedef connect1d<GEOSHAPE1D,ELMAP,NODEMAP> CONNECT1D;
    
    typedef pVect<point3d,NODEMAP>              NODESVECT;
    typedef pGraph<GEOELEMENT0D,ELMAP,NODEMAP>  GRAPH0D;
    typedef pGraph<GEOELEMENT1D,ELMAP,NODEMAP>  GRAPH1D;
    
    typedef pVect<Real,ELMAP>       CELLVOLUMESVECT;
    typedef pVect<point3d,ELMAP>    SUBCELLNSVECT;
    
    typedef pVect<UInt,ELMAP>       SUBCELL_TO_ELEMENT;
    typedef pVect<UInt,ELMAP>       SUBCELL_TO_NODE;
    typedef pVect<UInt,ELMAP>       INTERFACE_TO_ELEMENT;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 1;
    //@}
    
    /*! @name Internal flags */ //@{ 
  public:
    bool commDevLoaded;
    bool grid1dLoaded;
    bool connectGrid1dLoaded;
    bool connectCreated;
    bool finiteVolumeCreated;
    //@}
      
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<communicator> commDev;
    Teuchos::RCP<MESH1D>       grid1d;
    Teuchos::RCP<CONNECT1D>    connectGrid1d;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    geoMapInterface<GEOSHAPE1D>  refShape1d;
    NODESVECT dualPoints;
    GRAPH1D   subCells;
    GRAPH0D   interfaces;
    //@}
    
     /*! @name Finite Volumes data */ //@{ 
  public:
    CELLVOLUMESVECT subCellVolume;
    SUBCELLNSVECT   subCellN;
    //@}
    
    /*! @name Internal connectivity data */ //@{ 
  public:
    SUBCELL_TO_ELEMENT   subCellToElement;
    SUBCELL_TO_NODE      subCellToNode;  
    INTERFACE_TO_ELEMENT interfaceToElement;
    //@}
    
    /*! @name Constructors */ //@{ 
  public:
    dual1d();
    dual1d(const Teuchos::RCP<communicator> & CommDev);
    dual1d(communicator & CommDev);
    void setMesh1d(const Teuchos::RCP<MESH1D> & Grid2d);
    void setMesh1d(const MESH1D & Grid2d);
    void setConnect1d(const Teuchos::RCP<CONNECT1D> & Connect2d);
    void setConnect1d(const CONNECT1D & Connect2d);
    //@}
    
    /*! @name Funzioni interne e buildup */ //@{ 
  public:    
    /*! The dual nodes is the union of the original nodes + the elements barycentric nodes.
    Given the id of a geoItem (such as an element, face, etc) this function produces the local id of the dual node associated.
    \param flag = 1 the original vertices, = 2 the element barycenter
    \param lid the local numbering of the geoItem (such as an element, face, etc) */
    UInt getLocalNodeId(const UInt & flag, const UInt & lid) const;
      
    /*! Build dual mesh */
    void buildDualMesh();
    
    /*! Build finite volume data */
    template<typename GEOSHAPE2D>
    void buildFiniteVolumeData(const Teuchos::RCP<mesh2d<GEOSHAPE2D,ELMAP,NODEMAP> > & grid2d, const Teuchos::RCP<connect2d<GEOSHAPE2D,ELMAP,NODEMAP> > & connectGrid2d);
      
    /*! Build finite volume data */
    template<typename GEOSHAPE2D>
    void buildFiniteVolumeData(const mesh2d<GEOSHAPE2D,ELMAP,NODEMAP> & grid2d, const connect2d<GEOSHAPE2D,ELMAP,NODEMAP> & connectGrid2d);
      
    /*! Clear the finite volume data only */
    void clearFiniteVolume();
      
    /*! Clear all the data */
    void clear();
    //@}
};


template<typename ELMAP, typename NODEMAP>
dual1d<linearLine,ELMAP,NODEMAP>::
dual1d()
{
  commDevLoaded       = false;
  grid1dLoaded        = false;
  connectGrid1dLoaded = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
}

template<typename ELMAP, typename NODEMAP>
dual1d<linearLine,ELMAP,NODEMAP>::
dual1d(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded       = true;
  grid1dLoaded        = false;
  connectGrid1dLoaded = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
  
  commDev = CommDev;
}

template<typename ELMAP, typename NODEMAP>
dual1d<linearLine,ELMAP,NODEMAP>::
dual1d(communicator & CommDev)
{
  commDevLoaded       = true;
  grid1dLoaded        = false;
  connectGrid1dLoaded = false;
  connectCreated      = false;
  finiteVolumeCreated = false;
  
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename ELMAP, typename NODEMAP>
void
dual1d<linearLine,ELMAP,NODEMAP>::
setMesh1d(const Teuchos::RCP<MESH1D> & Grid1d)
{
  grid1dLoaded = true;
  grid1d       = Grid1d;
}

template<typename ELMAP, typename NODEMAP>
void
dual1d<linearLine,ELMAP,NODEMAP>::
setMesh1d(const MESH1D & Grid1d)
{
  grid1dLoaded = true;
  grid1d       = Teuchos::rcp(new MESH1D(Grid1d));
}

template<typename ELMAP, typename NODEMAP>
void
dual1d<linearLine,ELMAP,NODEMAP>::
setConnect1d(const Teuchos::RCP<CONNECT1D> & Connect1d)
{
  connectGrid1dLoaded = true;
  connectGrid1d       = Connect1d;
}

template<typename ELMAP, typename NODEMAP>
void
dual1d<linearLine,ELMAP,NODEMAP>::
setConnect1d(const CONNECT1D & Connect1d)
{
  connectGrid1dLoaded = true;
  connectGrid1d       = Teuchos::rcp(new CONNECT1D(Connect1d));
}

template<typename ELMAP, typename NODEMAP> 
UInt
dual1d<linearLine,ELMAP,NODEMAP>::
getLocalNodeId(const UInt & flag, const UInt & lid) const
{ 
  return( UInt(flag > 1) * grid1d->getNumNodes()
        + lid );
}

template<typename ELMAP, typename NODEMAP> 
void
dual1d<linearLine,ELMAP,NODEMAP>::
buildDualMesh()
{
  //Assert_________________________________________________________________________________________
  assert(commDevLoaded);
  assert(grid1dLoaded);
  assert(connectGrid1dLoaded);
  
  connectCreated = true;
  
  
  //Clearing_______________________________________________________________________________________
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
    
  subCellToElement.clear();
  subCellToNode.clear();
  interfaceToElement.clear();
  
  
  //Num local and global elements__________________________________________________________________
  UInt numLocalElements = grid1d->getNumElements();
  
  pMapGlobalManip<NODEMAP> vrMapManip(commDev);
  pMapGlobalManip<ELMAP>   elMapManip(commDev);
  
  pMap<NODEMAP> vrMap = grid1d->getNodes().getMapRef();
  pMap<ELMAP>   elMap = grid1d->getElements().getRowMap();
  
  UInt numGlobalNodes = vrMapManip.sizeG(vrMap);
  
  
  //Allocations____________________________________________________________________________________
  UInt pid = commDev->rank();
  
  UInt    id1, id2;
  point3d P1, P2, PM;
  
  GEOELEMENT1D subCell(true);
  GEOELEMENT0D interface(true);
  
  UInt numSubCells = 1;
  UInt numInterfaces = 1;
  
  UInt ba, locStray, elGid, geoId, numPoints;
  sVect<UInt> p(2);
  
  ELMAP   elMapItem;
  NODEMAP nodeMapItem;
  
  
  //Nodes construction_____________________________________________________________________________
  dualPoints = grid1d->getNodes();   //Orginal nodes
  numPoints  = dualPoints.size() + 1;
  
  for(UInt i=1; i <= numLocalElements; ++i)    //Nodi baricentrici
  {
    id1 = grid1d->getElementL(i).getCid(1);
    id2 = grid1d->getElementL(i).getCid(2);
    
    P1  = grid1d->getNodeL(id1);
    P2  = grid1d->getNodeL(id2);
    
    PM  = (P1 + P2) / 2.0;
    
    nodeMapItem.setPid(pid);
    nodeMapItem.setLid(numPoints);
    nodeMapItem.setGid(numGlobalNodes +  grid1d->getElements().getRowMapL(i).getGid());
    
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
    elGid = grid1d->getElements().getRowMapL(i).getGid();
    
    //Line nodes
    p(1) = grid1d->getElementL(i).getCid(1);
    p(2) = grid1d->getElementL(i).getCid(2);
    
    //Barycenter
    ba = getLocalNodeId(2,i);
    
    //Subcells creation
    locStray = 1;
    elGid    = grid1d->getElements().getRowMapL(i).getGid();
    geoId    = grid1d->getElements().getItemL(i).getGeoId();
    
    elMapItem.setPid(pid);
    subCell.setGeoId(geoId);
    
    //Subcell 1--------------------------------------------
    subCell.setCid(1,p(1));
    subCell.setCid(2,ba);
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(2 * (elGid - 1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(1));
    
    //Subcell 2--------------------------------------------
    subCell.setCid(1,ba);
    subCell.setCid(2,p(2));
    
    elMapItem.setLid(numSubCells);
    elMapItem.setGid(2 * (elGid - 1) + locStray);
    
    locStray++;
    numSubCells++;
    
    subCells.push_back(elMapItem,subCell);
    subCellToElement.push_back(elMapItem,i);
    subCellToNode.push_back(elMapItem,p(2));
    
    
    //Interfaces creation__________________________________________________________________________
    interface.setCid(1,ba);
    interface.setGeoId(geoId);
    
    elMapItem.setPid(pid);
    elMapItem.setLid(numInterfaces);
    elMapItem.setGid(elGid);
    
    numInterfaces++;
    locStray++;
    
    interfaces.push_back(elMapItem,interface);
    interfaceToElement.push_back(elMapItem, i);
  }
  

  //Map fixing for subCells
  pMap<ELMAP> tempSubCellsMap = subCells.getRowMap();
  nodesMapFixer(tempSubCellsMap,*commDev);
  
  subCells.setRowMap(tempSubCellsMap);
  subCells.setColMap(dualPoints.getMapRef());
  subCells.updateRowFinder();
  subCells.updateColFinder();
  
  subCellToElement.setMap(tempSubCellsMap);
  subCellToElement.updateFinder();

  subCellToNode.setMap(tempSubCellsMap);
  subCellToNode.updateFinder();
  
  //Checking for subCells
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> subCellsCheker(commDev);
  assert(subCellsCheker.check(subCells));
  
  
  //MapFixing for interfaces
  pMap<ELMAP> tempInterfacesMap = interfaces.getRowMap();
  nodesMapFixer(tempInterfacesMap,*commDev);
  
  interfaces.setRowMap(tempInterfacesMap);
  interfaces.setColMap(dualPoints.getMapRef());
  interfaces.updateRowFinder();
  interfaces.updateColFinder();
  
  interfaceToElement.setMap(tempInterfacesMap);
  interfaceToElement.updateFinder();
  
  //Checking interfaces map
  pGraphGlobalManip<GEOELEMENT0D,ELMAP,NODEMAP> interfacesCheker(commDev);
  assert(interfacesCheker.check(interfaces));
}

template<typename ELMAP, typename NODEMAP>
template<typename GEOSHAPE2D>
void
dual1d<linearLine,ELMAP,NODEMAP>::
buildFiniteVolumeData(const Teuchos::RCP<mesh2d<GEOSHAPE2D,ELMAP,NODEMAP> > & grid2d, const Teuchos::RCP<connect2d<GEOSHAPE2D,ELMAP,NODEMAP> > & connectGrid2d)
{
  //Allocations____________________________________________________________________________________
  geoMapInterface<GEOSHAPE1D> geo1d;
  geoMapSupport2d<GEOSHAPE2D> geo2d;
  
  UInt bVertId, belId, surfVertex, volVertex;
  UInt id1, id2, el1d, el2d, locEd;
  point3d P1, P2, N;
  point3d Ys = geo1d.getBarycenter();
  
  //2d data________________________________________________________________________________________
  sVect<UInt>  surfVertexToVertex(grid1d->getNumNodes());
  sVect<UInt>  surfElementToElement(grid1d->getNumElements());
  sVect<UInt>  surfElementToLocEdge(grid1d->getNumElements());
  
  for(UInt i=1; i <= grid2d->getNumNodes(); ++i)
  {
    if(connectGrid2d->getVertexIsBoundary(i))
    {
      bVertId = connectGrid2d->getVertexBVertex(i);
      assert(bVertId <= grid1d->getNumNodes());
      
      surfVertexToVertex(bVertId) = i;
    }
  }
  
  for(UInt i=1; i <= grid2d->getNumEdges(); ++i)
  {
    if(connectGrid2d->getEdgeIsBoundary(i))
    {
      belId = connectGrid2d->getEdgeBEdge(i);
      assert(belId <= grid1d->getNumElements());
      
      surfElementToElement(belId) = connectGrid2d->getEdgeToElement(i,1);
      surfElementToLocEdge(belId) = connectGrid2d->getEdgeToElement().getItemL(i).getSubLocId(1);
    }
  }
  
  //Clearing_______________________________________________________________________________________  
  subCellVolume.clear();
  subCellN.clear();
  
  //Resizing_______________________________________________________________________________________
  subCellVolume.resize(subCells.size());
  subCellN.resize(subCells.size());
  
  //Load maps______________________________________________________________________________________
  subCellVolume.setMap(subCells.getMapRef());
  subCellN.setMap(subCells.getMapRef());
  
  
  //Compute some quantities________________________________________________________________________
  for(UInt i=1; i <= subCells.size(); ++i)
  {
    //Surface interface data
    id1 = subCells(i).getCid(1);
    id2 = subCells(i).getCid(2);
    
    P1 = dualPoints(id1);
    P2 = dualPoints(id2);
    
    N = P1 - P2;
    subCellVolume(i) = N.norm2();
    
    //Line Normal
    surfVertex = subCellToNode(i);
    el1d       = subCellToElement(i);
    
    volVertex  = surfVertexToVertex(surfVertex);
    el2d       = surfElementToElement(el1d);
    locEd      = surfElementToLocEdge(el1d);
    
    geo2d.setPoints(grid2d->getElementNodesL(el2d));

    subCellN(i) = geo2d.normal(Ys,locEd);
  }
}

template<typename ELMAP, typename NODEMAP>
template<typename GEOSHAPE2D>
void
dual1d<linearLine,ELMAP,NODEMAP>::
buildFiniteVolumeData(const mesh2d<GEOSHAPE2D,ELMAP,NODEMAP> & grid2d, const connect2d<GEOSHAPE2D,ELMAP,NODEMAP> & connectGrid2d)
{
  //Allocations____________________________________________________________________________________
  geoMapInterface<GEOSHAPE1D> geo1d;
  geoMapSupport2d<GEOSHAPE2D> geo2d;
  
  UInt bVertId, belId, surfVertex, volVertex;
  UInt id1, id2, el1d, el2d, locEd;
  point3d P1, P2, N;
  point3d Ys = geo1d.getBarycenter();
  
  //2d data________________________________________________________________________________________
  sVect<UInt>  surfVertexToVertex(grid1d->getNumNodes());
  sVect<UInt>  surfElementToElement(grid1d->getNumElements());
  sVect<UInt>  surfElementToLocEdge(grid1d->getNumElements());
  
  for(UInt i=1; i <= grid2d.getNumNodes(); ++i)
  {
    if(connectGrid2d.getVertexIsBoundary(i))
    {
      bVertId = connectGrid2d.getVertexBVertex(i);
      assert(bVertId <= grid1d->getNumNodes());
      
      surfVertexToVertex(bVertId) = i;
    }
  }
  
  for(UInt i=1; i <= grid2d.getNumEdges(); ++i)
  {
    if(connectGrid2d.getEdgeIsBoundary(i))
    {
      belId = connectGrid2d.getEdgeBEdge(i);
      assert(belId <= grid1d->getNumElements());
      
      surfElementToElement(belId) = connectGrid2d.getEdgeToElement(i,1);
      surfElementToLocEdge(belId) = connectGrid2d.getEdgeToElement().getItemL(i).getSubLocId(1);
    }
  }
  
  //Clearing_______________________________________________________________________________________  
  subCellVolume.clear();
  subCellN.clear();
  
  //Resizing_______________________________________________________________________________________
  subCellVolume.resize(subCells.size());
  subCellN.resize(subCells.size());
  
  //Load maps______________________________________________________________________________________
  subCellVolume.setMap(subCells.getMapRef());
  subCellN.setMap(subCells.getMapRef());
  
  
  //Compute some quantities________________________________________________________________________
  for(UInt i=1; i <= subCells.size(); ++i)
  {
    //Surface interface data
    id1 = subCells(i).getCid(1);
    id2 = subCells(i).getCid(2);
    
    P1 = dualPoints(id1);
    P2 = dualPoints(id2);
    
    N = P1 - P2;
    subCellVolume(i) = N.norm2();
    
    //Line Normal
    surfVertex = subCellToNode(i);
    el1d       = subCellToElement(i);
    
    volVertex  = surfVertexToVertex(surfVertex);
    el2d       = surfElementToElement(el1d);
    locEd      = surfElementToLocEdge(el1d);
    
    geo2d.setPoints(grid2d.getElementNodesL(el2d));

    subCellN(i) = geo2d.normal(Ys,locEd);
  }
}

template<typename ELMAP, typename NODEMAP> 
void
dual1d<linearLine,ELMAP,NODEMAP>::
clearFiniteVolume()
{
  finiteVolumeCreated = false;
  
  subCellVolume.clear();
  subCellN.clear();
}      

template<typename ELMAP, typename NODEMAP> 
void
dual1d<linearLine,ELMAP,NODEMAP>::
clear()
{
  connectCreated      = false;
  finiteVolumeCreated = false;
  
  dualPoints.clear();
  subCells.clear();
  interfaces.clear();
  
  subCellToElement.clear();
  subCellToNode.clear();  
  interfaceToElement.clear();
  
  subCellVolume.clear();
  subCellN.clear();
}


#endif
