/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef DUALMESH3D_HPP
#define DUALMESH3D_HPP

#include "mesh2d.hpp"
#include "dual3d.hpp"

using namespace std;


/*! Interface for the dual mesh 3d, see also \c dual3d */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class dualMesh3d : public dual3d<GEOSHAPE,ELMAP,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef GEOSHAPE                        GEOSHAPE3D;
    typedef typename GEOSHAPE::GEOBSHAPE    GEOSHAPE2D;
    typedef typename GEOSHAPE2D::GEOBSHAPE  GEOSHAPE1D;
    
    typedef geoElement<GEOSHAPE3D>  GEOELEMENT3D;
    typedef geoElement<GEOSHAPE2D>  GEOELEMENT2D;
    typedef geoElement<GEOSHAPE1D>  GEOELEMENT1D;
    
    typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>    MESH2D;
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
    
    typedef dual3d<GEOSHAPE,ELMAP,NODEMAP> DUAL3D;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 3;
    //@}
  
     /*! @name Constructors */ //@{ 
  public:
    dualMesh3d();
    dualMesh3d(const Teuchos::RCP<communicator> & CommDev);
    dualMesh3d(communicator & CommDev);
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
    //@}
    
    /*! @name Get dual functions */ //@{
  public:
    const point3d      & getDualPointL(const UInt & lid) const;
    const point3d      & getDualPointG(const UInt & gid) const;
    const GEOELEMENT3D & getSubCellL(const UInt & lid) const;
    const GEOELEMENT3D & getSubCellG(const UInt & gid) const;
    const GEOELEMENT2D & getInterfaceL(const UInt & lid) const;
    const GEOELEMENT2D & getInterfaceG(const UInt & gid) const;
    //@}
    
    /*! @name Get dual connecting functions */ //@{
  public:
    const UInt & getSubCellToElement(const UInt & i) const;
    const UInt & getSubcellToNode(const UInt & i) const;  
    const UInt & getInterfaceToEdge(const UInt & i) const;
    //@}
    
    /*! @name Get finite volume functions */ //@{
  public:
    const Real    & getCellVolume(const UInt & i) const;
    const point3d & getEdgeN(const UInt & i) const;
    const Real    & getEdgeSurf(const UInt & i) const;
    //@}
    
    /*! @name Dump internal arrays */ //@{
  public:
    const NODESVECT & getDualPoints() const;
    const GRAPH3D   & getSubCells() const;
    const GRAPH2D   & getInterfaces() const;
    //@}
      
    /*! @name Dump internal connect arrays */ //@{
  public:
    const SUBCELL_TO_ELEMENT & getSubCellToElement() const;
    const SUBCELL_TO_NODE    & getSubCellToNode() const;
    const INTERFACE_TO_EDGE  & getInterfaceToEdge() const;
    //@}
    
    /*! @name Dual meshes printout */ //@{
  public:
    MESH3D getSubCellMesh();
    MESH2D getInterfaceMesh();
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
dualMesh3d() : dual3d<GEOSHAPE,ELMAP,NODEMAP>()
{ }

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
dualMesh3d(const Teuchos::RCP<communicator> & CommDev) : dual3d<GEOSHAPE,ELMAP,NODEMAP>(CommDev)
{ }

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
dualMesh3d(communicator & CommDev) : dual3d<GEOSHAPE,ELMAP,NODEMAP>(CommDev)
{ }


//_________________________________________________________________________________________________
// GET NUMBERS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getNumDualPoints() const
{
  assert(DUAL3D::dualCreated);
  return(DUAL3D::dualPoints.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getNumSubCells() const
{
  assert(DUAL3D::dualCreated);
  return(DUAL3D::subCells.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getNumInterfaces() const
{
  assert(DUAL3D::dualCreated);
  return(DUAL3D::interfaces.size());
}



//_________________________________________________________________________________________________
// GET NUMBERS CONNECT
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getNumSubCellToElement() const
{
  assert(DUAL3D::connectCreated);
  return(DUAL3D::subCellToElement.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getNumSubCellToNode() const
{
  assert(DUAL3D::connectCreated);
  return(DUAL3D::subCellToNode.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getNumInterfaceToEdge() const
{
  assert(DUAL3D::connectCreated);
  return(DUAL3D::interfaceToEdge.size());
}



//_________________________________________________________________________________________________
// GET DUAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getDualPointL(const UInt & lid) const
{
  assert(DUAL3D::dualCreated);
  assert(lid >= 1);
  assert(lid <= DUAL3D::dualPoints.size());
  return(DUAL3D::dualPoints.getL(lid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getDualPointG(const UInt & gid) const
{
  assert(DUAL3D::dualCreated);
  assert(DUAL3D::dualPoints.isG(gid));
  return(DUAL3D::dualPoints.getG(gid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT3D &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellL(const UInt & lid) const
{
  assert(DUAL3D::dualCreated);
  assert(lid >= 1);
  assert(lid <= DUAL3D::subCells.size());
  return(DUAL3D::subCells.getItemL(lid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT3D &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellG(const UInt & gid) const
{
  assert(DUAL3D::dualCreated);
  assert(DUAL3D::subCells.isColG(gid));
  return(DUAL3D::subCells.getItemG(gid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT2D &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceL(const UInt & lid) const
{
  assert(DUAL3D::dualCreated);
  assert(lid >= 1);
  assert(lid <= DUAL3D::interfaces.size());
  return(DUAL3D::interfaces.getItemL(lid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT2D &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceG(const UInt & gid) const
{
  assert(DUAL3D::dualCreated);
  assert(DUAL3D::interfaces.isColG(gid));
  return(DUAL3D::interfaces.getItemG(gid));
}



//_________________________________________________________________________________________________
// GET DUAL CONNECTING FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellToElement(const UInt & i) const
{
  assert(DUAL3D::connectCreated);
  assert(i >= 1);
  assert(i <= DUAL3D::subCellToElement.size());
  return(DUAL3D::subCellToElement(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getSubcellToNode(const UInt & i) const
{
  assert(DUAL3D::connectCreated);
  assert(i >= 1);
  assert(i <= DUAL3D::subcellToNode.size());
  return(DUAL3D::subcellToNode(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceToEdge(const UInt & i) const
{
  assert(DUAL3D::connectCreated);
  assert(i >= 1);
  assert(i <= DUAL3D::interfaceToEdge.size());
  return(DUAL3D::interfaceToEdge(i));
}



//_________________________________________________________________________________________________
// GET FINITE VOLUME FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const Real &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getCellVolume(const UInt & i) const
{
  assert(DUAL3D::finiteVolumeCreated);
  assert(i >= 1);
  assert(i <= DUAL3D::cellVolume.size());
  return(DUAL3D::cellVolume(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeN(const UInt & i) const
{
  assert(DUAL3D::finiteVolumeCreated);
  assert(i >= 1);
  assert(i <= DUAL3D::edgeN.size());
  return(DUAL3D::edgeN(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const Real &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgeSurf(const UInt & i) const
{
  assert(DUAL3D::finiteVolumeCreated);
  assert(i >= 1);
  assert(i <= DUAL3D::edgeSurf.size());
  return(DUAL3D::edgeSurf(i));
}



//_________________________________________________________________________________________________
// DUMP INTERNAL ARRAYS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::NODESVECT &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getDualPoints() const
{
  assert(DUAL3D::dualCreated);
  return(DUAL3D::dualPoints);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH3D &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCells() const
{
  assert(DUAL3D::dualCreated);
  return(DUAL3D::subCells);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH2D &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaces() const
{
  assert(DUAL3D::dualCreated);
  return(DUAL3D::interfaces);
}



//_________________________________________________________________________________________________
// DUMP INTERNAL CONNECT ARRAYS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::SUBCELL_TO_ELEMENT &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellToElement() const
{
  assert(DUAL3D::connectCreated);
  return(DUAL3D::subCellToElement);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::SUBCELL_TO_NODE &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellToNode() const
{
  assert(DUAL3D::connectCreated);
  return(DUAL3D::subCellToNode);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::INTERFACE_TO_EDGE &
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceToEdge() const
{
  assert(DUAL3D::connectCreated);
  return(DUAL3D::interfaceToEdge);
}



//_________________________________________________________________________________________________
// DUAL MESHES PRINTOUT
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
typename dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::MESH3D
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellMesh()
{
  assert(DUAL3D::dualCreated);
  
  MESH3D outMesh;
  outMesh.setNodes(DUAL3D::dualPoints);
  outMesh.setElements(DUAL3D::subCells);
  outMesh.transferMap();
  return(outMesh);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
typename dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::MESH2D
dualMesh3d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceMesh()
{
  assert(DUAL3D::dualCreated);
  
  MESH2D outMesh;
  outMesh.setNodes(DUAL3D::dualPoints);
  outMesh.setElements(DUAL3D::interfaces);
  outMesh.transferMap();
  return(outMesh);
}

#endif
