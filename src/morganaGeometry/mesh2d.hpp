/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MESH2D_HPP
#define MESH2D_HPP

#include "morganaTypes.hpp"

#include "pVect.hpp"
#include "pGraph.hpp"

#include "geoElement.hpp"
#include "geoMapInterface.hpp"


/*! The mesh 2d */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class mesh2d : public geoMapInterface<GEOSHAPE>
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
    
    typedef pVect<point3d,NODEMAP>              NODESVECT;
    typedef pVect<UInt,NODEMAP>                 BOOLVECT;
    
    typedef pGraph<GEOELEMENT2D,ELMAP,NODEMAP>  GRAPH2D;
    typedef pGraph<GEOELEMENT1D,ELMAP,NODEMAP>  GRAPH1D;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 2;
    //@}
    
    /*! @name Internal containers */ //@{
  public:
    NODESVECT nodes;
    BOOLVECT  isVertex;
    GRAPH2D   elements;
    GRAPH1D   edges;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    mesh2d();
    mesh2d(const NODESVECT & Nodes, const BOOLVECT & IsVertex, const GRAPH2D & Elements);
    mesh2d(const NODESVECT & Nodes, const GRAPH2D & Elements);
    mesh2d(const mesh2d & Mesh);
    mesh2d & operator=(const mesh2d &Mesh);
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
    void clearEdges();
    //@}
    
    /*! @name Counting functions */ //@{
  public:
    UInt getNumNodes() const;
    UInt getNumVertices() const;
    UInt getNumElements() const;
    UInt getNumEdges() const;
    //@}
    
    /*! @name Get lists funstions */ //@{
  public:
    const BOOLVECT  & getIsVertex();
    const NODESVECT & getNodes() const;
    const GRAPH2D   & getElements() const;
    const GRAPH1D   & getEdges() const;
    //@}
    
    /*! @name Get item functions */ //@{
  public:
    const point3d & getNodeL(const UInt & lid) const;
    const point3d & getNodeG(const UInt & gid) const;
    const GEOELEMENT2D & getElementL(const UInt & lid) const;
    const GEOELEMENT2D & getElementG(const UInt & gid) const;
    const GEOELEMENT1D & getEdgeL(const UInt & lid) const;
    const GEOELEMENT1D & getEdgeG(const UInt & gid) const;
    //@}
    
    /*! @name Get points for geoitem functions */ //@{
  public:
    sVect<point3d> getElementNodesL(const UInt & lid) const;
    sVect<point3d> getElementNodesG(const UInt & gid) const;
    sVect<point3d> getEdgeNodesL(const UInt & lid) const;
    sVect<point3d> getEdgeNodesG(const UInt & gid) const;
    //@}
    
    /*! @name Set functions */ //@{
  public:
    void setNodes(const NODESVECT & Points, const BOOLVECT & IsVertex);
    void setNodes(const NODESVECT & Points);
    void setElements(const GRAPH2D & Elements);
    void setEdges(const GRAPH1D & Edges);
    //@}
    
    /*! @name Set items functions */ //@{
  public:    
    void setNodeL(const UInt & lid, const point3d & p, const bool & IsVertex);
    void setNodeG(const UInt & gid, const point3d & p, const bool & IsVertex);
    void setElementL(const UInt & lid, const GEOELEMENT2D & element);
    void setElementG(const UInt & gid, const GEOELEMENT2D & element);
    void setEdgeL(const UInt & lid, const GEOELEMENT1D & edge);
    void setEdgeG(const UInt & gid, const GEOELEMENT1D & edge);
    //@}
    
    /*! @name Geometry ops */ //@{
  public:
    void scale(const Real & coeff);
    //@}
    
    /*! @name Printout */ //@{
  public:
    template<typename G, typename R, typename C>
    friend ostream & operator<<(ostream & f, const mesh2d<G,R,C> & V);
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
mesh2d()
{
  assert(GEOSHAPE::nDim == 2);
  
  numVertices = 0;
  verticesComputed = false;
  mapTransferred   = false;
  standard         = STDN;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
mesh2d(const NODESVECT & Nodes, const BOOLVECT & IsVertex, const GRAPH2D & Elements)
{
  assert(Nodes.size() == IsVertex.size());
  assert(GEOSHAPE::nDim == 2);
  
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
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
mesh2d(const NODESVECT & Nodes, const GRAPH2D & Elements)
{
  assert(GEOSHAPE::nDim == 2);
  
  nodes    = Nodes;
  elements = Elements;
  
  for(UInt i=1; i <= Nodes.size(); ++i)
  { isVertex.push_back(nodes.getMapL(i),UInt(true)); }
  
  computeNumVertices();
  transferMap();
  verticesComputed = true;
  mapTransferred   = true;
  standard         = STDN;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
mesh2d(const mesh2d & Mesh)
{
  assert(GEOSHAPE::nDim == 2);
  
  nodes    = Mesh.nodes;
  isVertex = Mesh.isVertex;
  elements = Mesh.elements;
  edges    = Mesh.edges;
  standard = Mesh.standard;
  
  computeNumVertices();
  transferMap();
  verticesComputed = true;
  mapTransferred   = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
mesh2d<GEOSHAPE,ELMAP,NODEMAP> &
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
operator=(const mesh2d &Mesh)
{
  nodes    = Mesh.nodes;
  isVertex = Mesh.isVertex;
  elements = Mesh.elements;
  edges    = Mesh.edges;
  standard = Mesh.standard;
  
  computeNumVertices();
  transferMap();
  verticesComputed = true;
  mapTransferred   = true;
  
  return *this;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
setMeshStandard(const MeshStandards & Standard)
{
  standard = Standard;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const MeshStandards &
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getMeshStandard()
{
  return(standard);
}



//_________________________________________________________________________________________________
// COMPUTING AND CLEARING
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
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
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
transferMap()
{
  elements.setColMap(nodes.getMapRef());
  edges.setColMap(nodes.getMapRef());
  
  mapTransferred = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
clear()
{
  nodes.clear();
  isVertex.clear();
  elements.clear();
  edges.clear();
  
  numVertices = 0;  
  verticesComputed = true;
  mapTransferred   = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
clearEdges()
{
  edges.clear();
  mapTransferred = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const bool &
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getMapTransferred() const
{
  return(mapTransferred);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const bool &
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getVerticesComputed() const
{
  return(verticesComputed);
}



//_________________________________________________________________________________________________
// COUNTING
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumNodes() const
{
  return(nodes.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumVertices() const
{
  assert(verticesComputed);
  return(numVertices);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumElements() const
{
  return(elements.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumEdges() const
{
  return(edges.size());
}



//_________________________________________________________________________________________________
// GETLISTS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename mesh2d<GEOSHAPE,ELMAP,NODEMAP>::BOOLVECT &
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getIsVertex()
{
  return(isVertex);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename mesh2d<GEOSHAPE,ELMAP,NODEMAP>::NODESVECT &
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getNodes() const
{
  return(nodes);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename mesh2d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH2D &
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getElements() const
{
  assert(mapTransferred);
  return(elements);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename mesh2d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH1D &
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getEdges() const
{
  assert(mapTransferred);
  return(edges);
}



//_________________________________________________________________________________________________
// GET ITEM
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getNodeL(const UInt & lid) const
{
  assert(lid <= nodes.size());
  return(nodes.getL(lid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getNodeG(const UInt & gid) const
{
  assert(nodes.isG(gid));
  return(nodes.getG(gid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const
typename mesh2d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT2D &
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getElementL(const UInt & lid) const
{
  assert(lid <= elements.size());
  return(elements.getItemL(lid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const
typename mesh2d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT2D &
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getElementG(const UInt & gid) const
{
  assert(elements.isRowG(gid));
  return(elements.getItemG(gid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename mesh2d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT1D &
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeL(const UInt & lid) const
{
  assert(lid <= edges.size());
  return(edges.getItemL(lid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename mesh2d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT1D &
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeG(const UInt & gid) const
{
  assert(edges.isRowG(gid));
  return(edges.getItemG(gid));
}



//_________________________________________________________________________________________________
// GET POINTS FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
sVect<point3d>
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getElementNodesL(const UInt & lid) const
{
  assert(mapTransferred);
  assert(lid <= elements.size());
  
  sVect<point3d> out(geoMapInterface<GEOSHAPE2D>::getNumPoints());
  
  for(UInt j=1; j <= geoMapInterface<GEOSHAPE2D>::getNumPoints(); ++j)
  {
    out(j) = nodes( elements.getCid_LL(lid,j) );
  }
  
  return(out);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
sVect<point3d>
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getElementNodesG(const UInt & gid) const
{
  assert(mapTransferred);
  assert(elements.isRowG(gid));
  
  sVect<point3d> out(geoMapInterface<GEOSHAPE2D>::getNumPoints());
  
  for(UInt j=1; j <= geoMapInterface<GEOSHAPE2D>::getNumPoints(); ++j)
  {
    out(j) = nodes( elements.getCid_GL(gid,j) );
  }
  
  return(out);  
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
sVect<point3d>
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeNodesL(const UInt & lid) const
{
  assert(mapTransferred);
  assert(lid <= edges.size());
  
  sVect<point3d> out(geoMapInterface<GEOSHAPE1D>::getNumPoints());
  
  for(UInt j=1; j <= geoMapInterface<GEOSHAPE1D>::getNumPoints(); ++j)
  {
    out(j) = nodes( edges.getCid_LL(lid,j) );
  }
  
  return(out);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
sVect<point3d>
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeNodesG(const UInt & gid) const
{
  assert(mapTransferred);
  assert(edges.isRowG(gid));
  
  sVect<point3d> out(geoMapInterface<GEOSHAPE1D>::getNumPoints());
  
  for(UInt j=1; j <= geoMapInterface<GEOSHAPE1D>::getNumPoints(); ++j)
  {
    out(j) = nodes( edges.getCid_GL(gid,j) );
  }
  
  return(out);
}



//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
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
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
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
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
setElements(const GRAPH2D & Elements)
{
  elements = Elements;
  transferMap();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
setEdges(const GRAPH1D & Edges)
{
  edges = Edges;
  transferMap();
}



//_________________________________________________________________________________________________
// SET ITEMS
//-------------------------------------------------------------------------------------------------    
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>   
void
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
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
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
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
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
setElementL(const UInt & lid, const GEOELEMENT2D & element)
{
  assert(lid <= elements.size());
  elements(lid) = element;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
setElementG(const UInt & gid, const GEOELEMENT2D & element)
{
  assert(elements.isRowG(gid));
  elements.getG(gid) = element;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
setEdgeL(const UInt & lid, const GEOELEMENT1D & edge)
{
  assert(lid <= edges.size());
  edges(lid) = edge;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
setEdgeG(const UInt & gid, const GEOELEMENT1D & edge)
{
  assert(edges.isRowG(gid));
  edges.getG(gid) = edge;
}



//_________________________________________________________________________________________________
// GEOMETRY OPS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
mesh2d<GEOSHAPE,ELMAP,NODEMAP>::
scale(const Real & coeff)
{
  for(UInt i=1; i <= nodes.size(); ++i)
  {
    nodes(i) = nodes(i) * coeff;
  }
}



//_________________________________________________________________________________________________
// PRINTOUT 
//-------------------------------------------------------------------------------------------------
template<typename G, typename R, typename C>
ostream & operator<<(ostream & f, const mesh2d<G,R,C> & V)
{
  f << "Mesh 2d" << endl;
  f << "Nodes    : " << V.nodes.size() << endl;
  f << "Elements : " << V.elements.size() << endl;
  f << "Edges    : " << V.edges.size() << endl;
  
  return(f);
}

#endif
