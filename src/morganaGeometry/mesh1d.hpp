/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MESH1D_HPP
#define MESH1D_HPP

#include "morganaTypes.hpp"

#include "pVect.hpp"
#include "pGraph.hpp"

#include "geoElement.hpp"
#include "geoMapInterface.hpp"


/*! The mesh 1d */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class mesh1d : public geoMapInterface<GEOSHAPE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef ELMAP    GRID_ELMAP;
    typedef NODEMAP  GRID_NODEMAP;
    typedef GEOSHAPE GRID_GEOSHAPE;
    
    typedef GEOSHAPE                            GEOSHAPE1D;
    typedef geoElement<GEOSHAPE1D>              GEOELEMENT1D;
    
    typedef pVect<point3d,NODEMAP>              NODESVECT;
    typedef pVect<UInt,NODEMAP>                 BOOLVECT;
    
    typedef pGraph<GEOELEMENT1D,ELMAP,NODEMAP>  GRAPH1D;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 1;
    //@}
    
    /*! @name Internal containers */ //@{
  public:
    NODESVECT nodes;
    BOOLVECT  isVertex;
    GRAPH1D   elements;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    mesh1d();
    mesh1d(const NODESVECT & Nodes, const BOOLVECT & IsVertex, const GRAPH1D & Elements);
    mesh1d(const NODESVECT & Nodes, const GRAPH1D & Elements);
    mesh1d(const mesh1d & Mesh);
    mesh1d & operator=(const mesh1d & Mesh);
    //@}
    
    /*! @name Internal flags and methods */ //@{
  public:
    UInt numVertices;
    bool verticesComputed;
    bool mapTransferred;
    MeshStandards standard;
    
    void computeNumVertices();
    void transferMap();
    const bool & getVerticesComputed() const;
    const bool & getMapTransferred() const;
    void setMeshStandard(const MeshStandards & Standard);
    const MeshStandards & getMeshStandard();
    //@}
    
     /*! @name Clearing */ //@{
  public:
    void clear();
    //@}
    
    /*! @name Counting functions */ //@{
  public:
    UInt getNumNodes() const;
    UInt getNumVertices() const;
    UInt getNumElements() const;
    //@}
    
    /*! @name Get lists funstions */ //@{
  public:
    const BOOLVECT  & getIsVertex();
    const NODESVECT & getNodes() const;
    const GRAPH1D   & getElements() const;
    //@}
    
    /*! @name Get item functions */ //@{
  public:
    const point3d & getNodeL(const UInt & lid) const;
    const point3d & getNodeG(const UInt & gid) const;
    const GEOELEMENT1D & getElementL(const UInt & lid) const;
    const GEOELEMENT1D & getElementG(const UInt & gid) const;
    //@}
    
    /*! @name Get points for geoitem functions */ //@{
  public:
    sVect<point3d> getElementNodesL(const UInt & lid) const;
    sVect<point3d> getElementNodesG(const UInt & gid) const;
    //@}
    
    /*! @name Set functions */ //@{
  public:
    void setNodes(const NODESVECT & Points, const BOOLVECT & IsVertex);
    void setNodes(const NODESVECT & Points);
    void setElements(const GRAPH1D & Elements);
    //@}
    
    /*! @name Set items functions */ //@{
  public:    
    void setNodeL(const UInt & lid, const point3d & p, const bool & IsVertex);
    void setNodeG(const UInt & gid, const point3d & p, const bool & IsVertex);
    void setElementL(const UInt & lid, const GEOELEMENT1D & element);
    void setElementG(const UInt & gid, const GEOELEMENT1D & element);
    //@}
    
    /*! @name Geometry ops */ //@{
  public:
    void scale(const Real & coeff);
    void offset(const point3d & P0);
    //@}
    
    /*! @name Printout */ //@{
  public:
    template<typename G, typename R, typename C>
    friend ostream & operator<<(ostream & f, const mesh1d<G,R,C> & V);
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
mesh1d()
{
  assert(GEOSHAPE::nDim == 1);
  
  numVertices = 0;
  verticesComputed = false;
  mapTransferred   = false;
  standard         = STDN;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
mesh1d(const NODESVECT & Nodes, const BOOLVECT & IsVertex, const GRAPH1D & Elements)
{
  assert(Nodes.size() == IsVertex.size());
  assert(GEOSHAPE::nDim == 1);
  
  nodes    = Nodes;
  isVertex = IsVertex;
  elements = Elements;
  
  computeNumVertices();
  transferMap();
  verticesComputed = true;
  mapTransferred   = true;
  standard         = STDN;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
mesh1d(const NODESVECT & Nodes, const GRAPH1D & Elements)
{
  assert(GEOSHAPE::nDim == 1);
  
  nodes    = Nodes;
  elements = Elements;
  
  for(UInt i=1; i <= Nodes.size(); ++i)
  {
    isVertex.push_back(true);
  }
  
  computeNumVertices();
  transferMap();
  verticesComputed = true;
  mapTransferred   = true;
  standard         = STDN;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
mesh1d(const mesh1d & Mesh)
{
  assert(GEOSHAPE::nDim == 1);
  
  nodes    = Mesh.nodes;
  isVertex = Mesh.isVertex;
  elements = Mesh.elements;
  standard = Mesh.standard;
  
  computeNumVertices();
  transferMap();
  verticesComputed = true;
  mapTransferred   = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
mesh1d<GEOSHAPE,ELMAP,NODEMAP> &
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
operator=(const mesh1d &Mesh)
{
  nodes    = Mesh.nodes;
  isVertex = Mesh.isVertex;
  elements = Mesh.elements;
  standard = Mesh.standard;
  
  computeNumVertices();
  transferMap();
  verticesComputed = true;
  mapTransferred   = true;
  
  return *this;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
setMeshStandard(const MeshStandards & Standard)
{
  standard = Standard;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const MeshStandards &
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getMeshStandard()
{
  return(standard);
}



//_________________________________________________________________________________________________
// COMPUTING AND CLEARING
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
computeNumVertices()
{
  numVertices = 0;
  
  for(UInt i=1; i<=isVertex.size(); ++i )
  {
    numVertices += UInt(isVertex(i));
  }
  
  verticesComputed = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
transferMap()
{
  elements.setColMap(nodes.getMapRef());
  
  mapTransferred = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
clear()
{
  nodes.clear();
  isVertex.clear();
  elements.clear();
  
  numVertices = 0;  
  verticesComputed = true;
  mapTransferred   = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const bool &
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getMapTransferred() const
{
  return(mapTransferred);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const bool &
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getVerticesComputed() const
{
  return(verticesComputed);
}



//_________________________________________________________________________________________________
// COUNTING
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumNodes() const
{
  return(nodes.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumVertices() const
{
  assert(verticesComputed);
  return(numVertices);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumElements() const
{
  return(elements.size());
}



//_________________________________________________________________________________________________
// GETLISTS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename mesh1d<GEOSHAPE,ELMAP,NODEMAP>::BOOLVECT &
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getIsVertex()
{
  return(isVertex);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename mesh1d<GEOSHAPE,ELMAP,NODEMAP>::NODESVECT &
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getNodes() const
{
  return(nodes);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename mesh1d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH1D &
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getElements() const
{
  assert(mapTransferred);
  return(elements);
}



//_________________________________________________________________________________________________
// GET ITEM
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getNodeL(const UInt & lid) const
{
  assert(lid <= nodes.size());
  return(nodes.getL(lid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getNodeG(const UInt & gid) const
{
  assert(nodes.isG(gid));
  return(nodes.getG(gid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const
typename mesh1d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT1D &
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getElementL(const UInt & lid) const
{
  assert(lid <= elements.size());
  return(elements.getItemL(lid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const
typename mesh1d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT1D &
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getElementG(const UInt & gid) const
{
  assert(elements.isRowG(gid));
  return(elements.getItemG(gid));
}



//_________________________________________________________________________________________________
// GET POINTS FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
sVect<point3d>
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getElementNodesL(const UInt & lid) const
{
  assert(mapTransferred);
  assert(lid <= elements.size());
  
  sVect<point3d> out(geoMapInterface<GEOSHAPE1D>::getNumPoints());
  
  for(UInt j=1; j <= geoMapInterface<GEOSHAPE1D>::getNumPoints(); ++j)
  {
    out(j) = nodes( elements.getCid_LL(lid,j) );
  }
  
  return(out);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
sVect<point3d>
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getElementNodesG(const UInt & gid) const
{
  assert(mapTransferred);
  assert(elements.isRowG(gid));
  
  sVect<point3d> out(geoMapInterface<GEOSHAPE1D>::getNumPoints());
  
  for(UInt j=1; j <= geoMapInterface<GEOSHAPE1D>::getNumPoints(); ++j)
  {
    out(j) = nodes( elements.getCid_GL(gid,j) );
  }
  
  return(out);  
}



//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
setNodes(const NODESVECT & Points, const BOOLVECT & IsVertex)
{
  nodes    = Points;
  isVertex = IsVertex;
  
  computeNumVertices();
  transferMap();
  verticesComputed = true;
  mapTransferred   = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
setNodes(const NODESVECT & Points)
{
  nodes = Points;
  
  isVertex.resize(nodes.size());
  
  for(UInt i=1; i <= nodes.size(); ++i)
  {
    isVertex(i) = true;
  }
  
  computeNumVertices();
  transferMap();
  verticesComputed = true;
  mapTransferred   = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
setElements(const GRAPH1D & Elements)
{
  elements = Elements;
  transferMap();
}



//_________________________________________________________________________________________________
// SET ITEMS
//-------------------------------------------------------------------------------------------------    
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>   
void
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
setNodeL(const UInt & lid, const point3d & p, const bool & IsVertex)
{
  assert(lid <= nodes.size());
  assert(nodes.size() == isVertex.size());
  
  nodes(lid)    = p;
  isVertex(lid) = IsVertex;
  
  verticesComputed = false;
  mapTransferred   = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>   
void
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
setNodeG(const UInt & gid, const point3d & p, const bool & IsVertex)
{
  assert(nodes.isRowG(gid));
  
  nodes.getG(gid)    = p;
  isVertex.getG(gid) = IsVertex;
  
  verticesComputed = false;
  mapTransferred   = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
setElementL(const UInt & lid, const GEOELEMENT1D & element)
{
  assert(lid <= elements.size());
  elements(lid) = element;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
setElementG(const UInt & gid, const GEOELEMENT1D & element)
{
  assert(elements.isRowG(gid));
  elements.getG(gid) = element;
}



//_________________________________________________________________________________________________
// GEOMETRY OPS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
scale(const Real & coeff)
{
  for(UInt i=1; i <= nodes.size(); ++i)
  {
    nodes(i) = nodes(i) * coeff;
  }
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh1d<GEOSHAPE,ELMAP,NODEMAP>::
offset(const point3d & P0)
{
  for(UInt i=1; i <= nodes.size(); ++i)
  {
    nodes(i) = nodes(i) + P0;
  }
}


//_________________________________________________________________________________________________
// PRINTOUT 
//-------------------------------------------------------------------------------------------------
template<typename G, typename R, typename C>
ostream & operator<<(ostream & f, const mesh1d<G,R,C> & V)
{
  f << "Mesh 1d" << endl;
  f << "Nodes    : " << V.nodes.size() << endl;
  f << "Elements : " << V.elements.size() << endl;
  
  return(f);
}


#endif
