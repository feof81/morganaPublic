/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MESHREFINEUNIFORM2D_HPP
#define MESHREFINEUNIFORM2D_HPP

#include "pMapGlobalManip.h"
#include "traitsGeometry.hpp"
#include "connect2d.hpp"


//________________________________________________________________________________________________
// GENERAL EMPTY CLASS
//------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class meshRefineUniform2d
{
    /*! @name Typedefs */ //@{
  public:
    typedef mesh2d<GEOSHAPE,ELMAP,NODEMAP>       MESH2D;
    typedef connect2d<GEOSHAPE,ELMAP,NODEMAP> CONNECT2D;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    meshRefineUniform2d();
    meshRefineUniform2d(const Teuchos::RCP<const communicator> & CommDev);
    meshRefineUniform2d(const communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    //@}
};


template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshRefineUniform2d<GEOSHAPE,ELMAP,NODEMAP>::
meshRefineUniform2d()
{
  commDevLoaded = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshRefineUniform2d<GEOSHAPE,ELMAP,NODEMAP>::
meshRefineUniform2d(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshRefineUniform2d<GEOSHAPE,ELMAP,NODEMAP>::
meshRefineUniform2d(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshRefineUniform2d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshRefineUniform2d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}


//________________________________________________________________________________________________
// LINEAR TETRA SPECIALIZATION
//------------------------------------------------------------------------------------------------
template<typename ELMAP, typename NODEMAP>
class meshRefineUniform2d<linearTriangle,ELMAP,NODEMAP>
{   
    /*! @name Typedefs */ //@{
  public:
    typedef linearTriangle                     GEOSHAPE;
    typedef mesh2d<GEOSHAPE,ELMAP,NODEMAP>       MESH2D;
    typedef connect2d<GEOSHAPE,ELMAP,NODEMAP> CONNECT2D;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    meshRefineUniform2d();
    meshRefineUniform2d(const Teuchos::RCP<const communicator> & CommDev);
    meshRefineUniform2d(const communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    //@}
    
    /*! @name Refine functions */ //@{
  public:
    void refineUniform(      Teuchos::RCP<MESH2D>    & newMesh2d,
                       const Teuchos::RCP<MESH2D>    & oldMesh2d,  
                       const Teuchos::RCP<CONNECT2D> & oldConnect2d);
    
    void refineUniform(      MESH2D    & newMesh2d,
                       const MESH2D    & oldMesh2d,  
                       const CONNECT2D & oldConnect2d);
    //@}
};


template<typename ELMAP, typename NODEMAP>
meshRefineUniform2d<linearTriangle,ELMAP,NODEMAP>::
meshRefineUniform2d()
{
  commDevLoaded = false;
}

template<typename ELMAP, typename NODEMAP>
meshRefineUniform2d<linearTriangle,ELMAP,NODEMAP>::
meshRefineUniform2d(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename ELMAP, typename NODEMAP>
meshRefineUniform2d<linearTriangle,ELMAP,NODEMAP>::
meshRefineUniform2d(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename ELMAP, typename NODEMAP>
void
meshRefineUniform2d<linearTriangle,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename ELMAP, typename NODEMAP>
void
meshRefineUniform2d<linearTriangle,ELMAP,NODEMAP>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename ELMAP, typename NODEMAP>
void
meshRefineUniform2d<linearTriangle,ELMAP,NODEMAP>::
refineUniform(      Teuchos::RCP<MESH2D>    & newMesh2d,
              const Teuchos::RCP<MESH2D>    & oldMesh2d,  
              const Teuchos::RCP<CONNECT2D> & oldConnect2d)
{
  //Assert-----------------------------------------------------------
  assert(commDevLoaded);
  
  //Refine-----------------------------------------------------------
  refineUniform(*newMesh2d,
                *oldMesh2d,  
                *oldConnect2d);
}

template<typename ELMAP, typename NODEMAP>
void
meshRefineUniform2d<linearTriangle,ELMAP,NODEMAP>::
refineUniform(      MESH2D    & newMesh2d,
              const MESH2D    & oldMesh2d,  
              const CONNECT2D & oldConnect2d)
{
  //Assert-----------------------------------------------------------
  assert(commDevLoaded);
  
  //Typedefs---------------------------------------------------------
  typedef typename MESH2D::NODESVECT NODESVECT;
  typedef typename MESH2D::GRAPH2D     GRAPH2D;
  typedef geoElement<GEOSHAPE>      GEOELEMENT;
  
  //New nodes--------------------------------------------------------
  NODESVECT nodes = oldMesh2d.getNodes();
  
  UInt base = nodes.size();
  sVect<point3d> edgeNodes;
  
  for(UInt ed=1; ed <= oldMesh2d.getNumEdges(); ed++)
  {
    edgeNodes = oldMesh2d.getEdgeNodesL(ed);
    nodes.push_back( (edgeNodes(1) + edgeNodes(2))/2.0, NODEMAP(base + ed,0));
  }
  
  //Global numbering nodes-------------------------------------------
  pVectGlobalManip<point3d,NODEMAP> nodeManip(commDev);
  nodeManip.buildGlobalNumbering(nodes);
  
  //New elements-----------------------------------------------------
  UInt ed1, ed2, ed3;
  UInt node1, node2, node3, node4, node5, node6;
  UInt elLid, geoId;
  GEOELEMENT tri(true);
  GRAPH2D newElements;
  
  for(UInt el2d=1; el2d <= oldMesh2d.getNumElements(); el2d++)
  {
    //Extract edges
    ed1 = oldConnect2d.getElementToEdge(el2d,1);
    ed2 = oldConnect2d.getElementToEdge(el2d,2);
    ed3 = oldConnect2d.getElementToEdge(el2d,3);
    
    //Compute node lids
    node1  = oldMesh2d.getElementL(el2d).getCid(1);
    node2  = oldMesh2d.getElementL(el2d).getCid(2);
    node3  = oldMesh2d.getElementL(el2d).getCid(3);
    node4  = base + ed1;
    node5  = base + ed2;
    node6  = base + ed3;
    
    //Extract geoId
    geoId = oldMesh2d.getElementL(el2d).getGeoId();
    tri.setGeoId(geoId);
    
    //Build new elements
    elLid = (el2d-1) * 4 + 1;
    tri(1) = node1; tri(2) = node4; tri(3) = node6;
    newElements.push_back(tri, ELMAP(elLid,0));
    
    elLid = (el2d-1) * 4 + 2;
    tri(1) = node6; tri(2) = node4;  tri(3) = node5;
    newElements.push_back(tri, ELMAP(elLid,0));
    
    elLid = (el2d-1) * 4 + 3;
    tri(1) = node4; tri(2) = node2;  tri(3) = node5;
    newElements.push_back(tri, ELMAP(elLid,0));
    
    elLid = (el2d-1) * 4 + 4;
    tri(1) = node6; tri(2) = node5; tri(3) = node3;
    newElements.push_back(tri, ELMAP(elLid,0));
  }
  
  newElements.setColMap(nodes.getMapRef());
  
  //Global numbering elements----------------------------------------
  pGraphGlobalManip<GEOELEMENT,ELMAP,NODEMAP> elManip(commDev);
  elManip.buildGlobalNumbering(newElements);
  
  //Load data--------------------------------------------------------
  newMesh2d.clear();
  newMesh2d.setNodes(nodes);
  newMesh2d.setElements(newElements);
  newMesh2d.transferMap();
  newMesh2d.setMeshStandard(oldMesh2d.getMeshStandard());
}


#endif
