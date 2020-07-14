/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef DUALMESH2D_HPP
#define DUALMESH2D_HPP

#include "mesh1d.hpp"
#include "dual2d.hpp"

using namespace std;


/*! Interface for the dual mesh 2d, see also \c dual2d */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class dualMesh2d : public dual2d<GEOSHAPE,ELMAP,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef GEOSHAPE                      GEOSHAPE2D;
    typedef typename GEOSHAPE::GEOBSHAPE  GEOSHAPE1D;
    typedef geoElement<GEOSHAPE2D>        GEOELEMENT2D;
    typedef geoElement<GEOSHAPE1D>        GEOELEMENT1D;
    
    typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>    MESH1D;
    typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>    MESH2D;
    
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
    
    typedef dual2d<GEOSHAPE,ELMAP,NODEMAP> DUAL2D;
    //@}
    
    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 2;
    //@}

    /*! @name Constructors */ //@{ 
  public:
    dualMesh2d();
    dualMesh2d(const Teuchos::RCP<communicator> & CommDev);
    dualMesh2d(communicator & CommDev);
    //@}
    
    /*! @name Get numbers */ //@{ 
  public:
    const UInt & getNumDualPoints() const;
    const UInt & getNumSubCells() const;
    const UInt & getNumInterfaces() const;
    //@}
    
    /*! @name Get numbers connect */ //@{ 
  public:
    const UInt & getNumSubCellToElement() const;
    const UInt & getNumSubCellToNode() const;
    const UInt & getNumInterfaceToEdge() const;
    const UInt & getNumInterfaceToElement() const;
    //@}
    
    /*! @name Get dual functions */ //@{
  public:
    const point3d      & getDualPointL(const UInt & lid) const;
    const point3d      & getDualPointG(const UInt & gid) const;
    const GEOELEMENT2D & getSubCellL(const UInt & lid) const;
    const GEOELEMENT2D & getSubCellG(const UInt & gid) const;
    const GEOELEMENT1D & getInterfaceL(const UInt & lid) const;
    const GEOELEMENT1D & getInterfaceG(const UInt & gid) const;
    //@}
    
    /*! @name Get dual connecting functions */ //@{
  public:
    const UInt & getSubCellToElement(const UInt & i) const;
    const UInt & getSubcellToNode(const UInt & i) const;  
    const UInt & getInterfaceToEdge(const UInt & i) const;
    const UInt & getInterfaceToElement(const UInt & i) const;
    //@}
    
    /*! @name Get finite volume functions */ //@{
  public:
    const Real    & getCellVolume(const UInt & i) const;
    const Real    & getSubCellVolume(const UInt & i) const;
    const point3d & getCellN(const UInt & i) const;
    const point3d & getEdgeN(const UInt & i) const;
    const Real    & getEdgeSurf(const UInt & i) const;
    //@}
    
    /*! @name Dump internal arrays */ //@{
  public:
    const NODESVECT & getDualPoints() const;
    const GRAPH2D   & getSubCells() const;
    const GRAPH1D   & getInterfaces() const;
    //@}
    
    /*! @name Dump internal connect arrays */ //@{
  public:
    const SUBCELL_TO_ELEMENT   & getSubCellToElement() const;
    const SUBCELL_TO_NODE      & getSubCellToNode() const;
    const INTERFACE_TO_EDGE    & getInterfaceToEdge() const;
    const INTERFACE_TO_ELEMENT & getInterfaceToElement() const;
    //@}
    
    /*! @name Dual meshes printout */ //@{
  public:
    MESH2D getSubCellMesh();
    MESH1D getInterfaceMesh();
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
dualMesh2d() : dual2d<GEOSHAPE,ELMAP,NODEMAP>()
{ }

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
dualMesh2d(const Teuchos::RCP<communicator> & CommDev) : dual2d<GEOSHAPE,ELMAP,NODEMAP>(CommDev)
{ }

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
dualMesh2d(communicator & CommDev) : dual2d<GEOSHAPE,ELMAP,NODEMAP>(CommDev)
{ }



//_________________________________________________________________________________________________
// GET NUMBERS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumDualPoints() const
{
  assert(DUAL2D::dualCreated);
  return(DUAL2D::dualPoints.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumSubCells() const
{
  assert(DUAL2D::dualCreated);
  return(DUAL2D::subCells.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumInterfaces() const
{
  assert(DUAL2D::dualCreated);
  return(DUAL2D::interfaces.size());
}



//_________________________________________________________________________________________________
// GET NUMBERS CONNECT
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumSubCellToElement() const
{
  assert(DUAL2D::connectCreated);
  return(DUAL2D::subCellToElement.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumSubCellToNode() const
{
  assert(DUAL2D::connectCreated);
  return(DUAL2D::subCellToNode.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumInterfaceToEdge() const
{
  assert(DUAL2D::connectCreated);
  return(DUAL2D::interfaceToEdge.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getNumInterfaceToElement() const
{
  assert(DUAL2D::connectCreated);
  return(DUAL2D::interfaceToElement.size());
}



//_________________________________________________________________________________________________
// GET DUAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getDualPointL(const UInt & lid) const
{
  assert(DUAL2D::dualCreated);
  assert(lid >= 1);
  assert(lid <= DUAL2D::dualPoints.size());
  return(DUAL2D::dualPoints.getL(lid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getDualPointG(const UInt & gid) const
{
  assert(DUAL2D::dualCreated);
  assert(DUAL2D::dualPoints.isG(gid));
  return(DUAL2D::dualPoints.getG(gid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT2D &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellL(const UInt & lid) const
{
  assert(DUAL2D::dualCreated);
  assert(lid >= 1);
  assert(lid <= DUAL2D::subCells.size());
  return(DUAL2D::subCells.getItemL(lid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT2D &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellG(const UInt & gid) const
{
  assert(DUAL2D::dualCreated);
  assert(DUAL2D::subCells.isColG(gid));
  return(DUAL2D::subCells.getItemG(gid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT1D &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceL(const UInt & lid) const
{
  assert(DUAL2D::dualCreated);
  assert(lid >= 1);
  assert(lid <= DUAL2D::interfaces.size());
  return(DUAL2D::interfaces.getItemL(lid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT1D &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceG(const UInt & gid) const
{
  assert(DUAL2D::dualCreated);
  assert(DUAL2D::interfaces.isColG(gid));
  return(DUAL2D::interfaces.getItemG(gid));
}



//_________________________________________________________________________________________________
// GET DUAL CONNECTING FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellToElement(const UInt & i) const
{
  assert(DUAL2D::connectCreated);
  assert(i >= 1);
  assert(i <= DUAL2D::subCellToElement.size());
  return(DUAL2D::subCellToElement(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getSubcellToNode(const UInt & i) const
{
  assert(DUAL2D::connectCreated);
  assert(i >= 1);
  assert(i <= DUAL2D::subcellToNode.size());
  return(DUAL2D::subcellToNode(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceToEdge(const UInt & i) const
{
  assert(DUAL2D::connectCreated);
  assert(i >= 1);
  assert(i <= DUAL2D::interfaceToEdge.size());
  return(DUAL2D::interfaceToEdge(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceToElement(const UInt & i) const
{
  assert(DUAL2D::connectCreated);
  assert(i >= 1);
  assert(i <= DUAL2D::interfaceToElement.size());
  return(DUAL2D::interfaceToElement(i));
}



//_________________________________________________________________________________________________
// GET FINITE VOLUME FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const Real &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getCellVolume(const UInt & i) const
{
  assert(DUAL2D::finiteVolumeCreated);
  assert(i >= 1);
  assert(i <= DUAL2D::cellVolume.size());
  return(DUAL2D::cellVolume(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const Real &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellVolume(const UInt & i) const
{
  assert(DUAL2D::finiteVolumeCreated);
  assert(i >= 1);
  assert(i <= DUAL2D::subCellVolume.size());
  return(DUAL2D::subCellVolume(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getCellN(const UInt & i) const
{
  assert(DUAL2D::finiteVolumeCreated);
  assert(i >= 1);
  assert(i <= DUAL2D::subCellN.size());
  return(DUAL2D::subCellN(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeN(const UInt & i) const
{
  assert(DUAL2D::finiteVolumeCreated);
  assert(i >= 1);
  assert(i <= DUAL2D::edgeN.size());
  return(DUAL2D::edgeN(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const Real &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeSurf(const UInt & i) const
{
  assert(DUAL2D::finiteVolumeCreated);
  assert(i >= 1);
  assert(i <= DUAL2D::edgeSurf.size());
  return(DUAL2D::edgeSurf(i));
}



//_________________________________________________________________________________________________
// DUMP INTERNAL ARRAYS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::NODESVECT &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getDualPoints() const
{
  assert(DUAL2D::dualCreated);
  return(DUAL2D::dualPoints);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH2D &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCells() const
{
  assert(DUAL2D::dualCreated);
  return(DUAL2D::subCells);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH1D &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaces() const
{
  assert(DUAL2D::dualCreated);
  return(DUAL2D::interfaces);
}

   
   
//_________________________________________________________________________________________________
// DUMP INTERNAL CONNECT ARRAYS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::SUBCELL_TO_ELEMENT &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellToElement() const
{
  assert(DUAL2D::connectCreated);
  return(DUAL2D::subCellToElement);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::SUBCELL_TO_NODE &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellToNode() const
{
  assert(DUAL2D::connectCreated);
  return(DUAL2D::subCellToNode);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::INTERFACE_TO_EDGE &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceToEdge() const
{
  assert(DUAL2D::connectCreated);
  return(DUAL2D::interfaceToEdge);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::INTERFACE_TO_ELEMENT &
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceToElement() const
{
  assert(DUAL2D::connectCreated);
  return(DUAL2D::interfaceToElement);
}



//_________________________________________________________________________________________________
// DUAL MESHES PRINTOUT
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
typename dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::MESH2D
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellMesh()
{
  assert(DUAL2D::dualCreated);
  
  MESH2D outMesh;
  outMesh.setNodes(DUAL2D::dualPoints);
  outMesh.setElements(DUAL2D::subCells);
  outMesh.transferMap();
  return(outMesh);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
typename dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::MESH1D
dualMesh2d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceMesh()
{
  assert(DUAL2D::dualCreated);
  
  MESH1D outMesh;
  outMesh.setNodes(DUAL2D::dualPoints);
  outMesh.setElements(DUAL2D::interfaces);
  outMesh.transferMap();
  return(outMesh);
}


#endif
