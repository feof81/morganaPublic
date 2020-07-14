/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef SREFINE3D_HPP
#define SREFINE3D_HPP

#include "mesh3d.hpp"
#include "connect3d.hpp"
#include "pointElement.hpp"


/*! Edge ordering sub-class */
struct edgeComp
{
  typedef linearTetra                    GEOSHAPE3D;
  typedef typename GEOSHAPE3D::GEOBSHAPE GEOSHAPE2D;
  typedef typename GEOSHAPE2D::GEOBSHAPE GEOSHAPE1D;
  typedef pointElement<GEOSHAPE1D> POINT_EDGE;
  
  bool operator() (const POINT_EDGE & lhs, const POINT_EDGE & rhs) const
  {
    static Real lengthLhs, lengthRhs;
    
    lengthLhs = point3d::norm2(lhs.getPoint(1) - lhs.getPoint(2));
    lengthRhs = point3d::norm2(rhs.getPoint(1) - rhs.getPoint(2));
    
    return( (lengthLhs  < lengthRhs) ||
            (lengthLhs == lengthRhs) * (lhs < rhs) );
  }
};



/*! General empty class */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class sRefine3d
{
};

/*! Linear Tetra specialization */
template<typename ELMAP, typename NODEMAP>
class sRefine3d<linearTetra, ELMAP, NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef linearTetra                    GEOSHAPE3D;
    typedef typename GEOSHAPE3D::GEOBSHAPE GEOSHAPE2D;
    typedef typename GEOSHAPE2D::GEOBSHAPE GEOSHAPE1D;
 
    typedef geoElement<GEOSHAPE3D>  GEOELEMENT3D;
    typedef geoElement<GEOSHAPE2D>  GEOELEMENT2D;
    
    typedef sVect<point3d>     NODESVECT;
    typedef sVect<GEOELEMENT3D>  GRAPH3D;
    typedef sVect<GEOELEMENT2D>  GRAPH2D;
    
    typedef sVect<NODEMAP> NODEMAPVECT;
    typedef sVect<ELMAP>   GRAPH3DMAP;
    typedef sVect<ELMAP>   GRAPH2DMAP;
    
    typedef pointElement<GEOSHAPE1D>                POINT_EDGE;
    typedef std::multimap<POINT_EDGE,UInt,edgeComp> EDGES_MULTIMAP;
    
    typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>       MESH2D;
    typedef mesh3d<GEOSHAPE3D,ELMAP,NODEMAP>       MESH3D;
    typedef connect3d<GEOSHAPE3D,ELMAP,NODEMAP> CONNECT3D;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    Real tollH;
    UInt rFactor;
    std::set<UInt> elementsTbr;
    
    NODEMAPVECT nodesMap3d, nodesMap2d;
    NODESVECT   nodes3d, nodes2d;
    
    GRAPH3DMAP elementsMap3d;
    GRAPH3D    elements3d;
    
    GRAPH2DMAP elementsMap2d;
    GRAPH2D    elements2d;
    
    EDGES_MULTIMAP edges3d, edges2d;
    //@}
    
    /*! @name Control flags */ //@{
  public:
    bool meshLoaded;
    //@}
    
    /*! @name Constructors and io functions */ //@{
  public:
    sRefine3d();
    void upload(const MESH3D & grid3d,
                const MESH2D & grid2d);
    
    void upload(const Teuchos::RCP<MESH3D> & grid3d,
                const Teuchos::RCP<MESH2D> & grid2d);
    
    void download(MESH3D & grid3d,
                  MESH2D & grid2d);
    
    void download(Teuchos::RCP<MESH3D> & grid3d,
                  Teuchos::RCP<MESH2D> & grid2d);
    //@}
    
    /*! @name Internal functions */ //@{
  public:
    void splitElement2d(const POINT_EDGE & E,
                        const UInt & newNodeId,
                        const UInt & elId);
    
    void splitElement3d(const POINT_EDGE & E,
                        const UInt & newNodeId,
                        const UInt & elId);
    
    void splitEdge2d(const POINT_EDGE & E);
    
    void splitEdge3d(const POINT_EDGE & E);
    
    POINT_EDGE getLongestEdge3d(const UInt & r,
                                const POINT_EDGE & E);
    //@}
    
     /*! @name External functions - Serial functions */ //@{
  public:
    void setRefinementParams(const Real & TollH,
                             const UInt & RFactor,
                             const sVect<UInt> & ElList);
    
    bool              getEdgeTbr(POINT_EDGE & edge); //Tolto un CONST
    sVect<POINT_EDGE> getEdgesTbr(const UInt & maxLength);
    bool              refineEdge(const POINT_EDGE & edge) const;
    UInt              refineEdges(const sVect<POINT_EDGE> & edges) const;
    bool              refineLeb(const UInt & numSteps);
    //@}
    
    /*! @name Get functions */ //@{
  public:
    const NODESVECT & getNodes3d() const;
    const NODESVECT & getNodes2d() const;
    const GRAPH3D   & getElements3d() const;
    const GRAPH2D   & getElements2d() const;
    //@}
    
    /*! @name Print functions */ //@{
  public:
    void print();
    void printElements();
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND IO FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename ELMAP, typename NODEMAP>
sRefine3d<linearTetra, ELMAP, NODEMAP>::
sRefine3d()
{
}

template<typename ELMAP, typename NODEMAP>
void
sRefine3d<linearTetra, ELMAP, NODEMAP>::
upload(const MESH3D & grid3d,
       const MESH2D & grid2d)
{
  //Control flags----------------------------------------------------
  meshLoaded  = true;
  
  //Resize-----------------------------------------------------------
  nodes2d.resize(grid2d.getNumNodes());
  nodesMap2d.resize(grid2d.getNumNodes());
  
  nodes3d.resize(grid3d.getNumNodes());
  nodesMap3d.resize(grid3d.getNumNodes());
  
  elements2d.resize(grid2d.getNumElements());
  elementsMap2d.resize(grid2d.getNumElements());
  
  elements3d.resize(grid3d.getNumElements());
  elementsMap3d.resize(grid3d.getNumElements());
  
  edges2d.clear();
  edges3d.clear();
  
  //Fill the elements and nodes--------------------------------------
  for(UInt i=1; i <= grid2d.getNumNodes(); ++i)
  {
    nodes2d(i)    = grid2d.getNodeL(i);
    nodesMap2d(i) = grid2d.getNodes().getMapL(i);
  }
  
  for(UInt i=1; i <= grid3d.getNumNodes(); ++i)
  {
    nodes3d(i)    = grid3d.getNodeL(i);
    nodesMap3d(i) = grid3d.getNodes().getMapL(i);
  }
  
  for(UInt i=1; i <= grid2d.getNumElements(); ++i)
  {
    elements2d(i)    = grid2d.getElementL(i);
    elementsMap2d(i) = grid2d.getElements().getRowMapL(i);
  }
  
  for(UInt i=1; i <= grid3d.getNumElements(); ++i)
  {
    elements3d(i)    = grid3d.getElementL(i);
    elementsMap3d(i) = grid3d.getElements().getRowMapL(i);
  }
  
  //Fill the edges of the 3d mesh------------------------------------
  UInt locP;
  POINT_EDGE edge;
  sVect<point3d> elPoints;
  geoMapInterface<GEOSHAPE3D> interface3d;
  pair<POINT_EDGE,UInt> paio;
  
  for(UInt i=1; i <= grid3d.getNumElements(); ++i)
  {
    //Extract the nodes
    elPoints = grid3d.getElementNodesL(i);
    
    //Build the edges
    for(UInt e=1; e <= interface3d.getNumEdges(); ++e)
    {
      //Edge building
      for(UInt p=1; p <= interface3d.getPointsOnEdge(e); p++)
      {
	locP = interface3d.edgeToPoint(e,p);
	edge.setPoint(p,elPoints(locP));
      }
      
      edge.setGeoId(0);
      edge.reorder();
      
      //Insert the edge
      paio.first  = edge;
      paio.second = i;
      
      edges3d.insert(paio);
    }
  }
  
  //Fill the edges of the 2d mesh------------------------------------
  geoMapInterface<GEOSHAPE2D> interface2d;
  
  for(UInt i=1; i <= grid2d.getNumElements(); ++i)
  {
    //Extract the nodes
    elPoints = grid2d.getElementNodesL(i);
    
    //Build the edges
    for(UInt e=1; e <= interface2d.getNumEdges(); ++e)
    {
      //Edge building
      for(UInt p=1; p <= interface2d.getPointsOnEdge(e); p++)
      {
	locP = interface2d.edgeToPoint(e,p);
	edge.setPoint(p,elPoints(locP));
      }
      
      edge.setGeoId(0);
      edge.reorder();
      
      //Insert the edge
      paio.first  = edge;
      paio.second = i;
      
      edges2d.insert(paio);
    }
  }
}

template<typename ELMAP, typename NODEMAP>
void
sRefine3d<linearTetra, ELMAP, NODEMAP>::
upload(const Teuchos::RCP<MESH3D> & grid3d,
       const Teuchos::RCP<MESH2D> & grid2d)
{
  upload(*grid3d,
         *grid2d);
}

template<typename ELMAP, typename NODEMAP>
void
sRefine3d<linearTetra, ELMAP, NODEMAP>::
download(MESH3D & grid3d,
         MESH2D & grid2d) {
  //Grid 3d--------------------------------------------------------------------
  grid3d.clear();
  
  //Points
  pVect<point3d,NODEMAP> nodeList3d;
  nodeList3d.reserve(nodes3d.size());
  
  for(UInt i=1; i <= nodes3d.size(); ++i)
  { nodeList3d.push_back(nodes3d(i), nodesMap3d(i)); }
  
  //Elements
  pGraph<GEOELEMENT3D,ELMAP,NODEMAP> elList3d;
  elList3d.reserve(elements3d.size());
  
  for(UInt i=1; i <= elements3d.size(); ++i)
  { elList3d.push_back(elements3d(i), elementsMap3d(i)); }
  
  grid3d.setNodes(nodeList3d);
  grid3d.setElements(elList3d);
  grid3d.transferMap();
  
  //Grid 2d--------------------------------------------------------------------
  grid2d.clear();
  
  //Points
  pVect<point3d,NODEMAP> nodeList2d;
  nodeList2d.reserve(nodes2d.size());
  
  for(UInt i=1; i <= nodes2d.size(); ++i)
  { nodeList2d.push_back(nodes2d(i),nodesMap2d(i)); }
  
  //Elements
  pGraph<GEOELEMENT2D,ELMAP,NODEMAP> elList2d;
  elList2d.reserve(elements2d.size());
  
  for(UInt i=1; i <= elements2d.size(); ++i)
  { elList2d.push_back(elements2d(i),elementsMap2d(i)); }
  
  grid2d.setNodes(nodeList2d);
  grid2d.setElements(elList2d);
  grid2d.transferMap();
}

template<typename ELMAP, typename NODEMAP>
void
sRefine3d<linearTetra, ELMAP, NODEMAP>::
download(Teuchos::RCP<MESH3D> & grid3d,
         Teuchos::RCP<MESH2D> & grid2d)
{
  download(*grid3d, *grid2d);
}


//_________________________________________________________________________________________________
// INTERNAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename ELMAP, typename NODEMAP>
void
sRefine3d<linearTetra, ELMAP, NODEMAP>::
splitElement2d(const POINT_EDGE & E,
               const UInt & newNodeId,
               const UInt & elId)
{ 
  //Typedef------------------------------------------------
  typedef typename std::multimap<POINT_EDGE,UInt>::iterator ITER;
  typedef std::pair<ITER,ITER> DOUBLE_ITER;
  
  //Assert-------------------------------------------------
  assert(elId > 0);
  assert(elId <= elements2d.size());
  assert(newNodeId > 0);
  assert(newNodeId <= nodes2d.size());
  assert(meshLoaded);
  
  //Order the local nodes----------------------------------
  UInt id;
  point3d P;
  point3d P1 = E.getPoint(1);
  point3d P2 = E.getPoint(2);
  sVect<UInt> edNodes, elNodes(GEOSHAPE2D::numPoints);
  UInt indexIn = 1, indexOut = 3;
  
  for(UInt i = 1; i <= GEOSHAPE2D::numPoints; ++i)
  {
    id = elements2d(elId).getCid(i);
    P  = nodes2d(id);
    
    if( (!(P != P1)) || (!(P != P2)) )
    {
      elNodes(indexIn) = id;
      indexIn++;
      
      edNodes.push_back(i);
    }
    else
    {
      elNodes(indexOut) = id;
      indexOut++;
    }
  }
  
  //Delete the edges--------------------------------------- 
  DOUBLE_ITER ret;
  
  POINT_EDGE ed;
  ed.setGeoId(0);
  
  //Ed 1
  ed.setPoint(1, nodes2d(elNodes(1)));
  ed.setPoint(2, nodes2d(elNodes(2)));
  ed.reorder();
  ret = edges2d.equal_range(ed);
  
  for(ITER it = ret.first; it != ret.second; ++it)
  { if(it->second == elId) { edges2d.erase(it); break; } }
  
  //Ed 2 
  ed.setPoint(1, nodes2d(elNodes(2)));
  ed.setPoint(2, nodes2d(elNodes(3)));
  ed.reorder();
  ret = edges2d.equal_range(ed);
  
  for(ITER it = ret.first; it != ret.second; ++it)
  { if(it->second == elId) { edges2d.erase(it); break; } }
  
  //Ed 3
  ed.setPoint(1, nodes2d(elNodes(1)));
  ed.setPoint(2, nodes2d(elNodes(3)));
  ed.reorder();
  ret = edges2d.equal_range(ed);

  for(ITER it = ret.first; it != ret.second; ++it)
  { if(it->second == elId) { edges2d.erase(it); break; } }
  
  //New elements cids
  GEOELEMENT2D elA(elements2d(elId)), elB(elements2d(elId));
  elA.setCid(edNodes(2), newNodeId);
  elB.setCid(edNodes(1), newNodeId);
  
  //Insert the elements
  elements2d(elId) = elA;
  elements2d.push_back(elB);
  
  ELMAP elMap = elementsMap2d(elId);
  elMap.setLid(elements2d.size());
  elMap.setGid(elements2d.size());
  elementsMap2d.push_back(elMap);
  
  //Insert edges
  pair<POINT_EDGE,UInt> paio;
  
  //First element
  paio.second = elId;
  
  ed.setPoint(1, nodes2d(elNodes(1)));
  ed.setPoint(2, nodes2d(elNodes(3)));
  ed.reorder();
  paio.first = ed;
  edges2d.insert(paio);
  
  ed.setPoint(1, nodes2d(elNodes(1)));
  ed.setPoint(2, nodes2d(newNodeId));
  ed.reorder();
  paio.first = ed;
  edges2d.insert(paio);
  
  ed.setPoint(1, nodes2d(elNodes(3)));
  ed.setPoint(2, nodes2d(newNodeId));
  ed.reorder();
  paio.first = ed;
  edges2d.insert(paio);
  
  //Second element
  paio.second = elements2d.size();
  
  ed.setPoint(1, nodes2d(elNodes(2)));
  ed.setPoint(2, nodes2d(elNodes(3)));
  ed.reorder();
  paio.first = ed;
  edges2d.insert(paio);
  
  ed.setPoint(1, nodes2d(elNodes(2)));
  ed.setPoint(2, nodes2d(newNodeId));
  ed.reorder();
  paio.first = ed;
  edges2d.insert(paio);
  
  ed.setPoint(1, nodes2d(elNodes(3)));
  ed.setPoint(2, nodes2d(newNodeId));
  ed.reorder();
  paio.first = ed;
  edges2d.insert(paio);
}

template<typename ELMAP, typename NODEMAP>
void
sRefine3d<linearTetra, ELMAP, NODEMAP>::
splitElement3d(const POINT_EDGE & E,
               const UInt & newNodeId,
               const UInt & elId)
{
  //Typedef------------------------------------------------
  typedef typename std::multimap<POINT_EDGE,UInt>::iterator ITER;
  typedef std::pair<ITER,ITER> DOUBLE_ITER;
  typedef std::pair<Real,POINT_EDGE> PAIR_LEB;
  
  //Assert-------------------------------------------------
  assert(elId > 0);
  assert(elId <= elements3d.size());
  assert(newNodeId > 0);
  assert(newNodeId <= nodes3d.size());
  assert(meshLoaded);
  
  //Order the local nodes----------------------------------
  UInt id;
  point3d P;
  point3d P1 = E.getPoint(1);
  point3d P2 = E.getPoint(2);
  sVect<UInt> edNodes, elNodes(GEOSHAPE3D::numPoints);
  UInt indexIn = 1, indexOut = 3;
  
  for(UInt i = 1; i <= GEOSHAPE3D::numPoints; ++i)
  {
    id = elements3d(elId).getCid(i);
    P  = nodes3d(id);
    
    if( (!(P != P1)) || (!(P != P2)) )
    {
      elNodes(indexIn) = id;
      indexIn++;
      
      edNodes.push_back(i);
    }
    else
    {
      elNodes(indexOut) = id;
      indexOut++;
    }
  }
  
  //Delete the edges---------------------------------------
  DOUBLE_ITER ret;
  
  POINT_EDGE ed;
  ed.setGeoId(0);
  
  //Ed 1
  ed.setPoint(1, nodes3d(elNodes(1)));
  ed.setPoint(2, nodes3d(elNodes(2)));
  ed.reorder();
  ret = edges3d.equal_range(ed);
  
  for(ITER it = ret.first; it != ret.second; ++it)
  { if(it->second == elId) { edges3d.erase(it); break; } }
  
  //Ed 2
  ed.setPoint(1, nodes3d(elNodes(1)));
  ed.setPoint(2, nodes3d(elNodes(3)));
  ed.reorder();
  ret = edges3d.equal_range(ed);
  
  for(ITER it = ret.first; it != ret.second; ++it)
  { if(it->second == elId) { edges3d.erase(it); break; } }
  
  //Ed 3
  ed.setPoint(1, nodes3d(elNodes(1)));
  ed.setPoint(2, nodes3d(elNodes(4)));
  ed.reorder();
  ret = edges3d.equal_range(ed);
  
  for(ITER it = ret.first; it != ret.second; ++it)
  { if(it->second == elId) { edges3d.erase(it); break; } }
  
  //Ed 4
  ed.setPoint(1, nodes3d(elNodes(2)));
  ed.setPoint(2, nodes3d(elNodes(3)));
  ed.reorder();
  ret = edges3d.equal_range(ed);
  
  for(ITER it = ret.first; it != ret.second; ++it)
  { if(it->second == elId) { edges3d.erase(it); break; } }
  
  //Ed 5
  ed.setPoint(1, nodes3d(elNodes(2)));
  ed.setPoint(2, nodes3d(elNodes(4)));
  ed.reorder();
  ret = edges3d.equal_range(ed);
  
  for(ITER it = ret.first; it != ret.second; ++it)
  { if(it->second == elId) { edges3d.erase(it); break; } }
  
  //Ed 6
  ed.setPoint(1, nodes3d(elNodes(3)));
  ed.setPoint(2, nodes3d(elNodes(4)));
  ed.reorder();
  ret = edges3d.equal_range(ed);
  
  for(ITER it = ret.first; it != ret.second; ++it)
  { if(it->second == elId) { edges3d.erase(it); break; } }
  
  //New elements cids
  GEOELEMENT3D elA(elements3d(elId)), elB(elements3d(elId));
  elA.setCid(edNodes(2), newNodeId);
  elB.setCid(edNodes(1), newNodeId);
  
  //Insert the elements 
  elements3d(elId) = elA;
  elements3d.push_back(elB);
  
  ELMAP elMap = elementsMap3d(elId);
  elMap.setLid(elements3d.size());
  elMap.setGid(elements3d.size());
  elementsMap3d.push_back(elMap);
  
  //Insert the elements in the tbr list
  if(elementsTbr.count(elId) != 0)
  { elementsTbr.insert(elements3d.size()); }
  
  //Insert edges
  PAIR_LEB        pairLeb;
  pair<POINT_EDGE,UInt> paio;
  
  //First element------------------------------------------
  paio.second = elId;
  
  //Edge 1
  ed.setPoint(1, nodes3d(elNodes(1)));
  ed.setPoint(2, nodes3d(newNodeId));
  ed.reorder();
  paio.first = ed;
  edges3d.insert(paio);
  
  pairLeb.first  = point3d::norm2(nodes3d(elNodes(1)) - nodes3d(newNodeId));
  pairLeb.second = ed;
  
  //Edge 2
  ed.setPoint(1, nodes3d(newNodeId));
  ed.setPoint(2, nodes3d(elNodes(3)));
  ed.reorder();
  paio.first = ed;
  edges3d.insert(paio);
  
  pairLeb.first  = point3d::norm2(nodes3d(elNodes(3)) - nodes3d(newNodeId));
  pairLeb.second = ed;
  
  //Edge 3
  ed.setPoint(1, nodes3d(elNodes(1)));
  ed.setPoint(2, nodes3d(elNodes(3)));
  ed.reorder();
  paio.first = ed;
  edges3d.insert(paio);
  
  pairLeb.first  = point3d::norm2(nodes3d(elNodes(1)) - nodes3d(elNodes(3)));
  pairLeb.second = ed;
  
  //Edge 4
  ed.setPoint(1, nodes3d(elNodes(1)));
  ed.setPoint(2, nodes3d(elNodes(4)));
  ed.reorder();
  paio.first = ed;
  edges3d.insert(paio);
  
  pairLeb.first  = point3d::norm2(nodes3d(elNodes(1)) - nodes3d(elNodes(4)));
  pairLeb.second = ed;
  
  //Edge 5
  ed.setPoint(1, nodes3d(newNodeId));
  ed.setPoint(2, nodes3d(elNodes(4)));
  ed.reorder();
  paio.first = ed;
  edges3d.insert(paio);
  
  pairLeb.first  = point3d::norm2(nodes3d(newNodeId) - nodes3d(elNodes(4)));
  pairLeb.second = ed;
  
  //Edge 6
  ed.setPoint(1, nodes3d(elNodes(3)));
  ed.setPoint(2, nodes3d(elNodes(4)));
  ed.reorder();
  paio.first = ed;
  edges3d.insert(paio);
  
  pairLeb.first  = point3d::norm2(nodes3d(elNodes(3)) - nodes3d(elNodes(4)));
  pairLeb.second = ed;
  
  //Second element-----------------------------------------
  paio.second = elements3d.size();
  
  //Edge 1
  ed.setPoint(1, nodes3d(newNodeId));
  ed.setPoint(2, nodes3d(elNodes(2)));
  ed.reorder();
  paio.first = ed;
  edges3d.insert(paio);
  
  pairLeb.first  = point3d::norm2(nodes3d(newNodeId) - nodes3d(elNodes(2)));
  pairLeb.second = ed;
  
  //Edge 2
  ed.setPoint(1, nodes3d(elNodes(2)));
  ed.setPoint(2, nodes3d(elNodes(3)));
  ed.reorder();
  paio.first = ed;
  edges3d.insert(paio);
  
  pairLeb.first  = point3d::norm2(nodes3d(elNodes(2)) - nodes3d(elNodes(3)));
  pairLeb.second = ed;
  
  //Edge 3
  ed.setPoint(1, nodes3d(newNodeId));
  ed.setPoint(2, nodes3d(elNodes(3)));
  ed.reorder();
  paio.first = ed;
  edges3d.insert(paio);
  
  pairLeb.first  = point3d::norm2(nodes3d(newNodeId) - nodes3d(elNodes(3)));
  pairLeb.second = ed;
  
  //Edge 4
  ed.setPoint(1, nodes3d(newNodeId));
  ed.setPoint(2, nodes3d(elNodes(4)));
  ed.reorder();
  paio.first = ed;
  edges3d.insert(paio);
  
  pairLeb.first  = point3d::norm2(nodes3d(newNodeId) - nodes3d(elNodes(4)));
  pairLeb.second = ed;
  
  //Edge 5
  ed.setPoint(1, nodes3d(elNodes(2)));
  ed.setPoint(2, nodes3d(elNodes(4)));
  ed.reorder();
  paio.first = ed;
  edges3d.insert(paio);
  
  pairLeb.first  = point3d::norm2(nodes3d(elNodes(2)) - nodes3d(elNodes(4)));
  pairLeb.second = ed;
  
  //Edge 6
  ed.setPoint(1, nodes3d(elNodes(3)));
  ed.setPoint(2, nodes3d(elNodes(4)));
  ed.reorder();
  paio.first = ed;
  edges3d.insert(paio);
  
  pairLeb.first  = point3d::norm2(nodes3d(elNodes(3)) - nodes3d(elNodes(4)));
  pairLeb.second = ed;
}

template<typename ELMAP, typename NODEMAP>
void
sRefine3d<linearTetra, ELMAP, NODEMAP>::
splitEdge2d(const POINT_EDGE & E)
{
  //Typedef------------------------------------------------
  typedef typename std::multimap<POINT_EDGE,UInt>::iterator ITER;
  typedef std::pair<ITER,ITER> DOUBLE_ITER;
  
  assert(meshLoaded);
  
  //Get elements-------------------------------------------
  DOUBLE_ITER ret;
  sVect<UInt> elStar;
  ret = edges2d.equal_range(E);
  
  for(ITER it = ret.first; it != ret.second; ++it)
  { elStar.push_back(it->second); }
  
  //Add node-----------------------------------------------
  point3d Pm = (E.getPoint(1) + E.getPoint(2)) / 2.0;  
  nodes2d.push_back(Pm);
  
  UInt nodeId = nodes2d.size();
  
  NODEMAP nodeMap = nodesMap2d(nodeId-1);
  nodeMap.setLid(nodeId);
  nodeMap.setGid(nodeId);
  nodesMap2d.push_back(nodeMap);
  
  //Split the elements-------------------------------------
  for(UInt i=1; i <= elStar.size(); ++i)
  { splitElement2d(E, nodeId, elStar(i)); }
}

template<typename ELMAP, typename NODEMAP>
void
sRefine3d<linearTetra, ELMAP, NODEMAP>::
splitEdge3d(const POINT_EDGE & E)
{
  //Typedef------------------------------------------------
  typedef typename std::multimap<POINT_EDGE,UInt>::iterator ITER;
  typedef std::pair<ITER,ITER> DOUBLE_ITER;
  
  assert(meshLoaded);
  
  //Get elements-------------------------------------------
  DOUBLE_ITER ret;
  sVect<UInt> elStar;
  ret = edges3d.equal_range(E);
  
  for(ITER it = ret.first; it != ret.second; ++it)
  { elStar.push_back(it->second); }
  
  //Add node-----------------------------------------------
  point3d Pm = (E.getPoint(1) + E.getPoint(2)) / 2.0;  
  nodes3d.push_back(Pm);
  
  UInt nodeId = nodes3d.size();
  
  NODEMAP nodeMap = nodesMap3d(nodeId-1);
  nodeMap.setLid(nodeId);
  nodeMap.setGid(nodeId);
  nodesMap3d.push_back(nodeMap);
  
  //Split the elements-------------------------------------
  for(UInt i=1; i <= elStar.size(); ++i)
  { splitElement3d(E, nodeId, elStar(i)); }
}

template<typename ELMAP, typename NODEMAP>
typename sRefine3d<linearTetra, ELMAP, NODEMAP>::POINT_EDGE
sRefine3d<linearTetra, ELMAP, NODEMAP>::
getLongestEdge3d(const UInt & r,
                 const POINT_EDGE & E)
{
  //Typedef------------------------------------------------
  typedef typename std::multimap<POINT_EDGE,UInt>::iterator ITER;
  typedef std::pair<ITER,ITER> DOUBLE_ITER;
  
  //Assert-------------------------------------------------
  assert(meshLoaded);
  assert(edges3d.count(E) != 0);
  
  //Alloc--------------------------------------------------
  Real length;
  UInt elId, id1, id2, id3, id4;
  point3d P1, P2, P3, P4;
  POINT_EDGE outEdge(E), newEdge(E);
  DOUBLE_ITER ret;
  
  //Main loop----------------------------------------------
  for(UInt i=1; i <= r; ++i)
  {
    //data of the edge
    length  = point3d::norm2(outEdge.getPoint(1) - outEdge.getPoint(2));
    newEdge = outEdge;
    
    //List of the elements
    ret = edges3d.equal_range(outEdge);
    assert(edges3d.count(outEdge) != 0);
  
    //Loop on the elements
    for(ITER it = ret.first; it != ret.second; ++it)
    {
      elId = it->second;
    
      id1 = elements3d(elId).getCid(1);
      id2 = elements3d(elId).getCid(2);
      id3 = elements3d(elId).getCid(3);
      id4 = elements3d(elId).getCid(4);
    
      P1 = nodes3d(id1);
      P2 = nodes3d(id2);
      P3 = nodes3d(id3);
      P4 = nodes3d(id4);
    
      //Ed1
      if(length < point3d::norm2(P1 - P2))
      {
        newEdge.setPoint(1,P1);
        newEdge.setPoint(2,P2);
        newEdge.reorder();
	length = point3d::norm2(P1 - P2);
      }
    
      //Ed2
      if(length < point3d::norm2(P1 - P3))
      {
        newEdge.setPoint(1,P1);
        newEdge.setPoint(2,P3);
        newEdge.reorder();
	length = point3d::norm2(P1 - P3);
      }
    
      //Ed3
      if(length < point3d::norm2(P1 - P4))
      {
        newEdge.setPoint(1,P1);
        newEdge.setPoint(2,P4);
        newEdge.reorder();
	length = point3d::norm2(P1 - P4);
      }
    
      //Ed4
      if(length < point3d::norm2(P2 - P3))
      {
        newEdge.setPoint(1,P2);
        newEdge.setPoint(2,P3);
        newEdge.reorder();
	length = point3d::norm2(P2 - P3);
      }
    
      //Ed5
      if(length < point3d::norm2(P2 - P4))
      {
        newEdge.setPoint(1,P2);
        newEdge.setPoint(2,P4);
        newEdge.reorder();
	length = point3d::norm2(P2 - P4);
      }
    
      //Ed6
      if(length < point3d::norm2(P3 - P4))
      {
        newEdge.setPoint(1,P3);
        newEdge.setPoint(2,P4);
        newEdge.reorder();
	length = point3d::norm2(P3 - P4);
      }
    }
  
    if(! (newEdge != outEdge))
    { return(newEdge);}
    
    outEdge = newEdge;
  }
  
  return(outEdge);
}


//_________________________________________________________________________________________________
// EXTERNAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename ELMAP, typename NODEMAP>
void
sRefine3d<linearTetra, ELMAP, NODEMAP>::
setRefinementParams(const Real & TollH,
                    const UInt & RFactor,
                    const sVect<UInt> & ElList)
{
  assert(meshLoaded);
  
  tollH   = TollH;
  rFactor = RFactor;
  
  elementsTbr.clear();
  
  for(UInt i=1; i <= ElList.size(); ++i)
  { elementsTbr.insert(ElList(i)); }
}

template<typename ELMAP, typename NODEMAP>
bool
sRefine3d<linearTetra, ELMAP, NODEMAP>::
getEdgeTbr(POINT_EDGE & edge)
{
  typedef typename EDGES_MULTIMAP::reverse_iterator ITER;
  
  Real length;
  UInt elId;
  
  for(ITER iter = edges3d.rbegin(); iter != edges3d.rend(); ++iter)
  {
    edge = iter->first;
    elId = iter->second;
    
    length = point3d::norm2(edge.getPoint(1) - edge.getPoint(2));
    
    if( (elementsTbr.count(elId) != 0) && (length >= tollH))
    {
      edge = getLongestEdge3d(rFactor,edge);
      return(true);
    }
  }
  
  return(false);
}

template<typename ELMAP, typename NODEMAP>
sVect<typename sRefine3d<linearTetra, ELMAP, NODEMAP>::POINT_EDGE>
sRefine3d<linearTetra, ELMAP, NODEMAP>::
getEdgesTbr(const UInt & maxLength)
{
  typedef typename EDGES_MULTIMAP::reverse_iterator ITER;
  
  Real length;
  UInt elId;
  POINT_EDGE edge;
  sVect<POINT_EDGE> out;
  
  for(ITER iter = edges3d.rbegin(); iter != edges3d.rend(); ++iter)
  {
    edge = iter->first;
    elId = iter->second;
    
    length = point3d::norm2(edge.getPoint(1) - edge.getPoint(2));
    
    if(length >= tollH)
    {
      if( elementsTbr.count(elId) != 0)
      {
        edge = getLongestEdge3d(rFactor,edge);
        out.push_back(edge);
      }
    }
    else
    { break; }
    
    if(out.size() >= maxLength)
    { break; }
  }
  
  return(out);
}

template<typename ELMAP, typename NODEMAP>
bool
sRefine3d<linearTetra, ELMAP, NODEMAP>::
refineEdge(const POINT_EDGE & edge) const
{
  if(edges3d.count(edge))
  {
    splitEdge3d(edge);
    
    if(edges2d.count(edge))
    { splitEdge2d(edge); }
    
    return(true);
  }
  
  return(false);    
}

template<typename ELMAP, typename NODEMAP>
UInt
sRefine3d<linearTetra, ELMAP, NODEMAP>::
refineEdges(const sVect<POINT_EDGE> & edges) const
{
  UInt numSplit = 0;
  POINT_EDGE edge;
  
  for(UInt i=1; i <= edges.size(); ++i)
  {
    edge = edges(i);
    
    if(edges3d.count(edge))
    {
      splitEdge3d(edge);
      
      if(edges2d.count(edge))
      { splitEdge2d(edge); }
      
      numSplit++;
    }
  }
  
  return(numSplit);
}

template<typename ELMAP, typename NODEMAP>
bool
sRefine3d<linearTetra, ELMAP, NODEMAP>::
refineLeb(const UInt & numSteps)
{
  bool flag = false;
  POINT_EDGE edge;
  
  for(UInt i=1; i <= numSteps; ++i)
  {
    flag = getEdgeTbr(edge);
    
    if(flag)
    { refineEdge(edge); }
    else
    { return(flag); }
  }
  
  return(flag);
}



//_________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename ELMAP, typename NODEMAP>
const typename sRefine3d<linearTetra, ELMAP, NODEMAP>::NODESVECT &
sRefine3d<linearTetra, ELMAP, NODEMAP>::
getNodes3d() const
{
  return(nodes3d);
}

template<typename ELMAP, typename NODEMAP>
const typename sRefine3d<linearTetra, ELMAP, NODEMAP>::NODESVECT &
sRefine3d<linearTetra, ELMAP, NODEMAP>::
getNodes2d() const
{
  return(nodes2d);
}

template<typename ELMAP, typename NODEMAP> 
const typename sRefine3d<linearTetra, ELMAP, NODEMAP>::GRAPH3D &
sRefine3d<linearTetra, ELMAP, NODEMAP>::
getElements3d() const
{
  return(elements3d);
}

template<typename ELMAP, typename NODEMAP>
const typename sRefine3d<linearTetra, ELMAP, NODEMAP>::GRAPH2D &
sRefine3d<linearTetra, ELMAP, NODEMAP>::
getElements2d() const
{
  return(elements2d);
}


//_________________________________________________________________________________________________
// PRINT FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename ELMAP, typename NODEMAP>
void
sRefine3d<linearTetra, ELMAP, NODEMAP>::
print()
{
  typedef typename std::multimap<POINT_EDGE,UInt>::iterator ITER;

  cout << "POINT_EDGES 2D --------------------------------------------" << endl;
  for(ITER iter = edges2d.begin(); iter != edges2d.end(); ++iter)
  { cout << iter->second << endl << iter->first << endl; }
  
  cout << "POINT_EDGES 3D --------------------------------------------" << endl;
  for(ITER iter = edges3d.begin(); iter != edges3d.end(); ++iter)
  { cout << iter->second << endl << iter->first << endl; }
}

template<typename ELMAP, typename NODEMAP>
void
sRefine3d<linearTetra, ELMAP, NODEMAP>::
printElements()
{
  cout << "NODES 2D --------------------------------------------" << endl;
  for(UInt i=1; i <= nodes2d.size(); ++i)
  { cout << nodes2d(i); }
  cout << endl;
  
  cout << "ELEMENTS 2D -----------------------------------------" << endl;
  for(UInt i=1; i <= elements2d.size(); ++i)
  { cout << elements2d(i) << endl; }
  cout << endl;
  
  cout << "NODES 3D --------------------------------------------" << endl;
  for(UInt i=1; i <= nodes3d.size(); ++i)
  { cout << nodes3d(i); }
  cout << endl;
  
  cout << "ELEMENTS 3D -----------------------------------------" << endl;
  for(UInt i=1; i <= elements3d.size(); ++i)
  { cout << elements3d(i) << endl; }
  cout << endl;
}

#endif
