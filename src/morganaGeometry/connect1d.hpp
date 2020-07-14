/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef CONNECT1D_HPP
#define CONNECT1D_HPP

#include <map>
#include <set>

#include "pMapItem.h"
#include "pMapGlobalManip.h"

#include "pGraph.hpp"
#include "pGraphItem.h"
#include "pGraphItemOriented.h"
#include "pGraphItemSubLoc.h"
#include "pGraphGlobalManip.hpp"

#include "mesh1d.hpp"

using namespace std;


/*! Topological connection information of a 1d mesh */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class connect1d
{
    /*! @name Typedefs */ //@{
  public:
    typedef ELMAP    GRID_ELMAP;
    typedef NODEMAP  GRID_NODEMAP;
    typedef GEOSHAPE GRID_GEOSHAPE;
    
    typedef GEOSHAPE                GEOSHAPE1D;
    typedef geoElement<GEOSHAPE1D>  GEOELEMENT1D;
    
    typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>  MESH1D;
    
    typedef pGraph<pGraphItem,NODEMAP,NODEMAP>  VERTEX_TO_VERTEX;
    typedef pGraph<pGraphItem,NODEMAP,ELMAP>    VERTEX_TO_ELEMENT;
    typedef pGraph<pGraphItem,ELMAP,ELMAP>      ELEMENT_TO_ELEMENT;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 1;
    //@}
    
    /*! @name Internal control flags */ //@{
  public:
    bool connectCreated;
    bool commDevLoaded;
    bool grid1dLoaded;
    //@}
    
     /*! @name Links */ //@{
  public:
    Teuchos::RCP<const communicator> commDev;
    Teuchos::RCP<MESH1D> grid1d;
    //@}
    
    /*! @name Internal graphs */ //@{
  public:
    VERTEX_TO_VERTEX   vertexToVertex;
    VERTEX_TO_ELEMENT  vertexToElement;
    ELEMENT_TO_ELEMENT elementToElement;
    //@}
    
    /*! @name Reference shapes */ //@{
  public:
    geoMapInterface<GEOSHAPE1D> refShape1d;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    connect1d();
    connect1d(const Teuchos::RCP<const communicator> & CommDev);
    connect1d(const communicator & CommDev);
    connect1d(const connect1d & C);
    connect1d & operator=(const connect1d & C);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    void setMesh1d(const Teuchos::RCP<MESH1D> & Grid1d);
    //@}
    
    /*! @name Computing funcions */ //@{
  public:
    void buildConnectivity();
    void clearConnectivity();
    void clear();
    //@}
    
    /*! @name Get num functions */ //@{
  public:
    UInt getNumVertexToVertex(const UInt & i) const;
    UInt getNumVertexToElement(const UInt & i) const;
    UInt getNumElementToElement(const UInt & i) const;
    //@}
    
    /*! @name Get volume functions */ //@{
  public:
    const UInt & getVertexToVertex(const UInt & i, const UInt & j) const;
    const UInt & getVertexToElement(const UInt & i, const UInt & j) const;
    const UInt & getElementToElement(const UInt & i, const UInt & j) const;
    //@}
    
    /*! @name Dump graphs functions */ //@{
  public:
    const VERTEX_TO_VERTEX   & getVertexToVertex() const;
    const VERTEX_TO_ELEMENT  & getVertexToElement() const;
    const ELEMENT_TO_ELEMENT & getElementToElement() const;
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
connect1d()
{
  assert(GEOSHAPE::nDim == 1);
  
  connectCreated = false;
  commDevLoaded  = false;
  grid1dLoaded   = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
connect1d(const Teuchos::RCP<const communicator> & CommDev)
{
  assert(GEOSHAPE::nDim == 1);
  
  connectCreated = false;
  commDevLoaded  = true;
  grid1dLoaded  = false;
  
  commDev = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
connect1d(const communicator & CommDev)
{
  assert(GEOSHAPE::nDim == 1);
  
  connectCreated = false;
  commDevLoaded  = true;
  grid1dLoaded   = false;
  
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
connect1d(const connect1d & C)
{
  assert(GEOSHAPE::nDim == 1);
  
  connectCreated = C.connectCreated;
  commDevLoaded  = C.commDevLoaded;
  grid1dLoaded   = C.grid1dLoaded;

  commDev = C.commDev;
  grid1d  = C.grid1d;
  
  vertexToVertex   = C.vertexToVertex;
  vertexToElement  = C.vertexToElement;
  elementToElement = C.elementToElement;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
connect1d<GEOSHAPE,ELMAP,NODEMAP> &
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
operator=(const connect1d & C)
{
  assert(GEOSHAPE::nDim == 1);
  
  connectCreated = C.connectCreated;
  commDevLoaded  = C.commDevLoaded;
  grid1dLoaded   = C.grid1dLoaded;

  commDev = C.commDev;
  grid1d  = C.grid1d;
  
  vertexToVertex   = C.vertexToVertex;
  vertexToElement  = C.vertexToElement;
  elementToElement = C.elementToElement;
  
  return *this;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
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
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
buildConnectivity()
{
  //Control logic----------------------------------------------------
  assert(grid1dLoaded);
  connectCreated = true;
  
  //Clearing---------------------------------------------------------
  vertexToVertex.clear();
  vertexToElement.clear();
  elementToElement.clear();
  
  //Other allocations------------------------------------------------ 
  UInt id1, id2;
  
  //Other connecting allocations-------------------------------------
  pGraphItem  vertexToVertexItem;
  pGraphItem  elementToElementItem;
  pGraphItem  vertexToElementItem;
  
  
  //Vertex to Vertex connectivity------------------------------------
  vertexToVertex.resize(grid1d->getNumNodes());
  
  for(UInt i=1; i <= grid1d->getNumElements(); ++i)
  {    
    id1  = grid1d->getElementL(i).getCid(1);
    id2  = grid1d->getElementL(i).getCid(2);
    
    vertexToVertexItem = vertexToVertex.getItemL(id1);
    vertexToVertexItem.push_back(id2);
    vertexToVertex.getItemL(id1) = vertexToVertexItem;
    
    vertexToVertexItem = vertexToVertex.getItemL(id2);
    vertexToVertexItem.push_back(id1);
    vertexToVertex.getItemL(id2) = vertexToVertexItem;
  }
  
  vertexToVertex.setRowMap(grid1d->getNodes().getMapRef());
  vertexToVertex.setColMap(grid1d->getNodes().getMapRef());
  
  vertexToVertex.updateRowFinder();
  vertexToVertex.updateColFinder();
  
  
  //Vertex to element connectivity-----------------------------------
  vertexToElement.resize(grid1d->getNumNodes());
  
  for(UInt i=1; i <= grid1d->getNumElements(); ++i)
  {
    for(UInt j=1; j <= refShape1d.getNumVertices(); j++)
    {
      id1 = grid1d->getElementL(i).getCid(j);
      
      vertexToElementItem = vertexToElement.getItemL(id1);
      vertexToElementItem.push_back(i);
      vertexToElement.getItemL(id1) = vertexToElementItem;
    }
  }
  
  vertexToElement.setRowMap(grid1d->getNodes().getMapRef());
  vertexToElement.setColMap(grid1d->getElements().getMapRef());
  
  vertexToElement.updateRowFinder();
  vertexToElement.updateColFinder();
  
  
  //Element to element connectivity----------------------------------
  elementToElement.resize(grid1d->getNumElements());
  
  for(UInt i=1; i <= grid1d->getNumNodes(); ++i)
  {
    vertexToElementItem = vertexToElement.getItemL(i);
    
    if(vertexToElementItem.size() == 2)
    {
      id1 = vertexToElementItem.getCid(1);
      id2 = vertexToElementItem.getCid(2);
      
      
      elementToElementItem = elementToElement.getItemL(id1);
      elementToElementItem.push_back(id2);
      elementToElement.getItemL(id1) = elementToElementItem;
      
      elementToElementItem = elementToElement.getItemL(id2);
      elementToElementItem.push_back(id1);
      elementToElement.getItemL(id2) = elementToElementItem;
    }
  }
  
  elementToElement.setRowMap(grid1d->getElements().getMapRef());
  elementToElement.setColMap(grid1d->getElements().getMapRef());
  
  elementToElement.updateRowFinder();
  elementToElement.updateColFinder();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
clearConnectivity()
{
  connectCreated = false;
  
  vertexToVertex.clear();
  vertexToElement.clear();
  elementToElement.clear();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
clear()
{
  connectCreated = false;
  
  vertexToVertex.clear();
  vertexToElement.clear();
  elementToElement.clear();
}



//_________________________________________________________________________________________________
// GET NUM FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumVertexToVertex(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= vertexToVertex.rowSizeL());
  return(vertexToVertex.getItemL(i).size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumVertexToElement(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= vertexToElement.size());
  return(vertexToElement.getItemL(i).size());
}


template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumElementToElement(const UInt & i) const
{
  assert(connectCreated);
  assert(i >= 1);
  assert(i <= elementToElement.rowSizeL());
  return(elementToElement.getItemL(i).size());
}



//_________________________________________________________________________________________________
// GET VOLUME FUINCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
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
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexToElement(const UInt & i, const UInt & j) const
{
  assert(connectCreated);
  assert(i >=1); assert(j >= 1);
  assert(i <= vertexToElement.size());
  assert(j <= vertexToElement.getItemL(i).size());
  return(vertexToElement.getCid_LL(i,j));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
getElementToElement(const UInt & i, const UInt & j) const
{
  assert(connectCreated);
  assert(i >=1); assert(j >= 1);
  assert(i <= elementToElement.rowSizeL());
  assert(j <= elementToElement.getItemL(i).size());
  return(elementToElement.getCid_LL(i,j));
}


//_________________________________________________________________________________________________
// DUMP GRAPH FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect1d<GEOSHAPE,ELMAP,NODEMAP>::VERTEX_TO_VERTEX &
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexToVertex() const
{
  assert(connectCreated);
  return(vertexToVertex);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect1d<GEOSHAPE,ELMAP,NODEMAP>::VERTEX_TO_ELEMENT &
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
getVertexToElement() const
{
  assert(connectCreated);
  return(vertexToElement);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename connect1d<GEOSHAPE,ELMAP,NODEMAP>::ELEMENT_TO_ELEMENT &
connect1d<GEOSHAPE,ELMAP,NODEMAP>::
getElementToElement() const
{
  assert(connectCreated);
  return(elementToElement);
}


#endif
