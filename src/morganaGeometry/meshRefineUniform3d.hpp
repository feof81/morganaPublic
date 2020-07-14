/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MESHREFINEUNIFORM3D_HPP
#define MESHREFINEUNIFORM3D_HPP

#include "pMapGlobalManip.h"
#include "traitsGeometry.hpp"
#include "connect3d.hpp"


//________________________________________________________________________________________________
// GENERAL EMPTY CLASS
//------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class meshRefineUniform3d
{
    /*! @name Typedefs */ //@{
  public:
    typedef mesh3d<GEOSHAPE,ELMAP,NODEMAP>       MESH3D;
    typedef connect3d<GEOSHAPE,ELMAP,NODEMAP> CONNECT3D;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    meshRefineUniform3d();
    meshRefineUniform3d(const Teuchos::RCP<const communicator> & CommDev);
    meshRefineUniform3d(const communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    //@}
};


template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshRefineUniform3d<GEOSHAPE,ELMAP,NODEMAP>::
meshRefineUniform3d()
{
  commDevLoaded = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshRefineUniform3d<GEOSHAPE,ELMAP,NODEMAP>::
meshRefineUniform3d(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshRefineUniform3d<GEOSHAPE,ELMAP,NODEMAP>::
meshRefineUniform3d(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshRefineUniform3d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshRefineUniform3d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}



//________________________________________________________________________________________________
// LINEAR TETRA SPECIALIZATION
//------------------------------------------------------------------------------------------------
template<typename ELMAP, typename NODEMAP>
class meshRefineUniform3d<linearTetra,ELMAP,NODEMAP>
{   
    /*! @name Typedefs */ //@{
  public:
    typedef linearTetra                        GEOSHAPE;
    typedef mesh3d<GEOSHAPE,ELMAP,NODEMAP>       MESH3D;
    typedef connect3d<GEOSHAPE,ELMAP,NODEMAP> CONNECT3D;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    meshRefineUniform3d();
    meshRefineUniform3d(const Teuchos::RCP<const communicator> & CommDev);
    meshRefineUniform3d(const communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    //@}
    
    /*! @name Refine functions */ //@{
  public:
    void refineUniform(      Teuchos::RCP<MESH3D>    & newMesh3d,
                       const Teuchos::RCP<MESH3D>    & oldMesh3d,  
                       const Teuchos::RCP<CONNECT3D> & oldConnect3d);
    
    void refineUniform(      MESH3D    & newMesh3d,
                       const MESH3D    & oldMesh3d,  
                       const CONNECT3D & oldConnect3d);
    //@}
};


template<typename ELMAP, typename NODEMAP>
meshRefineUniform3d<linearTetra,ELMAP,NODEMAP>::
meshRefineUniform3d()
{
  commDevLoaded = false;
}

template<typename ELMAP, typename NODEMAP>
meshRefineUniform3d<linearTetra,ELMAP,NODEMAP>::
meshRefineUniform3d(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename ELMAP, typename NODEMAP>
meshRefineUniform3d<linearTetra,ELMAP,NODEMAP>::
meshRefineUniform3d(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename ELMAP, typename NODEMAP>
void
meshRefineUniform3d<linearTetra,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename ELMAP, typename NODEMAP>
void
meshRefineUniform3d<linearTetra,ELMAP,NODEMAP>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename ELMAP, typename NODEMAP>
void
meshRefineUniform3d<linearTetra,ELMAP,NODEMAP>::
refineUniform(      Teuchos::RCP<MESH3D>    & newMesh3d,
              const Teuchos::RCP<MESH3D>    & oldMesh3d,  
              const Teuchos::RCP<CONNECT3D> & oldConnect3d)
{
  //Assert-----------------------------------------------------------
  assert(commDevLoaded);
  
  //Refine-----------------------------------------------------------
  refineUniform(*newMesh3d,
                *oldMesh3d,  
                *oldConnect3d);
}

template<typename ELMAP, typename NODEMAP>
void
meshRefineUniform3d<linearTetra,ELMAP,NODEMAP>::
refineUniform(      MESH3D    & newMesh3d,
              const MESH3D    & oldMesh3d,  
              const CONNECT3D & oldConnect3d)
{
  //Assert-----------------------------------------------------------
  assert(commDevLoaded);
  
  //Typedefs---------------------------------------------------------
  typedef typename MESH3D::NODESVECT NODESVECT;
  typedef typename MESH3D::GRAPH3D     GRAPH3D;
  typedef geoElement<GEOSHAPE>      GEOELEMENT;
  
  //New nodes--------------------------------------------------------
  NODESVECT nodes = oldMesh3d.getNodes();
  
  UInt base = nodes.size();
  sVect<point3d> edgeNodes;
  
  for(UInt ed=1; ed <= oldMesh3d.getNumEdges(); ed++)
  {
    edgeNodes = oldMesh3d.getEdgeNodesL(ed);
    nodes.push_back( (edgeNodes(1) + edgeNodes(2))/2.0, NODEMAP(base + ed,0));
  }
  
  //Global numbering nodes-------------------------------------------
  pVectGlobalManip<point3d,NODEMAP> nodeManip(commDev);
  nodeManip.buildGlobalNumbering(nodes);
  
  //New elements-----------------------------------------------------
  UInt ed1, ed2, ed3, ed4, ed5, ed6;
  UInt node1, node2, node3, node4, node5, node6, node7, node8, node9, node10;
  UInt elLid, geoId;
  GEOELEMENT tet(true);
  GRAPH3D newElements;
  
  for(UInt el3d=1; el3d <= oldMesh3d.getNumElements(); el3d++)
  {
    //Extract edges
    ed1 = oldConnect3d.getElementToEdge(el3d,1);
    ed2 = oldConnect3d.getElementToEdge(el3d,2);
    ed3 = oldConnect3d.getElementToEdge(el3d,3);
    ed4 = oldConnect3d.getElementToEdge(el3d,4);
    ed5 = oldConnect3d.getElementToEdge(el3d,5);
    ed6 = oldConnect3d.getElementToEdge(el3d,6);
    
    //Compute node lids
    node1  = oldMesh3d.getElementL(el3d).getCid(1);
    node2  = oldMesh3d.getElementL(el3d).getCid(2);
    node3  = oldMesh3d.getElementL(el3d).getCid(3);
    node4  = oldMesh3d.getElementL(el3d).getCid(4);
    node5  = base + ed1;
    node6  = base + ed2;
    node7  = base + ed3;
    node8  = base + ed4;
    node9  = base + ed5;
    node10 = base + ed6;
    
    //Extract geoId
    geoId = oldMesh3d.getElementL(el3d).getGeoId();
    tet.setGeoId(geoId);
    
    //Build new elements
    elLid = (el3d-1) * 8 + 1;
    tet(1) = node1; tet(2) = node8; tet(3) = node5; tet(4) = node7;
    newElements.push_back(tet, ELMAP(elLid,0));
    
    elLid = (el3d-1) * 8 + 2;
    tet(1) = node8; tet(2) = node4;  tet(3) = node9; tet(4) = node10;
    newElements.push_back(tet, ELMAP(elLid,0));
    
    elLid = (el3d-1) * 8 + 3;
    tet(1) = node5; tet(2) = node9;  tet(3) = node2; tet(4) = node6;
    newElements.push_back(tet, ELMAP(elLid,0));
    
    elLid = (el3d-1) * 8 + 4;
    tet(1) = node7; tet(2) = node10; tet(3) = node6; tet(4) = node3;
    newElements.push_back(tet, ELMAP(elLid,0));
    
    elLid = (el3d-1) * 8 + 5;
    tet(1) = node5; tet(2) = node7;  tet(3) = node8; tet(4) = node10;
    newElements.push_back(tet, ELMAP(elLid,0));
    
    elLid = (el3d-1) * 8 + 6;
    tet(1) = node5; tet(2) = node8;  tet(3) = node9; tet(4) = node10;
    newElements.push_back(tet, ELMAP(elLid,0));
    
    elLid = (el3d-1) * 8 + 7;
    tet(1) = node5; tet(2) = node9;  tet(3) = node6; tet(4) = node10;
    newElements.push_back(tet, ELMAP(elLid,0));
    
    elLid = (el3d-1) * 8 + 8;
    tet(1) = node5; tet(2) = node6;  tet(3) = node7; tet(4) = node10;
    newElements.push_back(tet, ELMAP(elLid,0));
  }
  
  newElements.setColMap(nodes.getMapRef());
  
  //Global numbering elements----------------------------------------
  pGraphGlobalManip<GEOELEMENT,ELMAP,NODEMAP> elManip(commDev);
  elManip.buildGlobalNumbering(newElements);
  
  //Load data--------------------------------------------------------
  newMesh3d.clear();
  newMesh3d.setNodes(nodes);
  newMesh3d.setElements(newElements);
  newMesh3d.transferMap();
  newMesh3d.setMeshStandard(oldMesh3d.getMeshStandard());
}

#endif
