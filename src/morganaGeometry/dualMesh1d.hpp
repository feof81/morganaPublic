/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef DUALMESH1D_HPP
#define DUALMESH1D_HPP

#include "mesh1d.hpp"
#include "dual1d.hpp"

using namespace std;


/*! Interface for the dual mesh 1d, see also \c dual1d */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class dualMesh1d : public dual1d<GEOSHAPE,ELMAP,NODEMAP>
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
    
    typedef dual1d<GEOSHAPE,ELMAP,NODEMAP> DUAL1D;
    //@}
    
    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 1;
    //@}

    /*! @name Constructors */ //@{ 
  public:
    dualMesh1d();
    dualMesh1d(const Teuchos::RCP<communicator> & CommDev);
    dualMesh1d(communicator & CommDev);
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
    const UInt & getNumInterfaceToElement() const;
    //@}
    
    /*! @name Get dual functions */ //@{
  public:
    const point3d      & getDualPointL(const UInt & lid) const;
    const point3d      & getDualPointG(const UInt & gid) const;
    const GEOELEMENT1D & getSubCellL(const UInt & lid) const;
    const GEOELEMENT1D & getSubCellG(const UInt & gid) const;
    const GEOELEMENT0D & getInterfaceL(const UInt & lid) const;
    const GEOELEMENT0D & getInterfaceG(const UInt & gid) const;
    //@}
    
    /*! @name Get dual connecting functions */ //@{
  public:
    const UInt & getSubCellToElement(const UInt & i) const;
    const UInt & getSubcellToNode(const UInt & i) const;  
    const UInt & getInterfaceToElement(const UInt & i) const;
    //@}
    
    /*! @name Get finite volume functions */ //@{
  public:
    const Real    & getSubCellVolume(const UInt & i) const;
    const point3d & getCellN(const UInt & i) const;
    //@}
    
    /*! @name Dump internal arrays */ //@{
  public:
    const NODESVECT & getDualPoints() const;
    const GRAPH1D   & getSubCells() const;
    const GRAPH0D   & getInterfaces() const;
    //@}
    
    /*! @name Dump internal connect arrays */ //@{
  public:
    const SUBCELL_TO_ELEMENT   & getSubCellToElement() const;
    const SUBCELL_TO_NODE      & getSubCellToNode() const;
    const INTERFACE_TO_ELEMENT & getInterfaceToElement() const;
    //@}
    
    /*! @name Dual meshes printout */ //@{
  public:
    MESH1D getSubCellMesh();
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
dualMesh1d() : dual1d<GEOSHAPE,ELMAP,NODEMAP>()
{ }

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
dualMesh1d(const Teuchos::RCP<communicator> & CommDev) : dual1d<GEOSHAPE,ELMAP,NODEMAP>(CommDev)
{ }

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
dualMesh1d(communicator & CommDev) : dual1d<GEOSHAPE,ELMAP,NODEMAP>(CommDev)
{ }


//_________________________________________________________________________________________________
// GET NUMBERS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumDualPoints() const
{
  assert(DUAL1D::connectCreated);
  return(DUAL1D::dualPoints.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumSubCells() const
{
  assert(DUAL1D::connectCreated);
  return(DUAL1D::subCells.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumInterfaces() const
{
  assert(DUAL1D::connectCreated);
  return(DUAL1D::interfaces.size());
}



//_________________________________________________________________________________________________
// GET NUMBERS CONNECT
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumSubCellToElement() const
{
  assert(DUAL1D::connectCreated);
  return(DUAL1D::subCellToElement.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumSubCellToNode() const
{
  assert(DUAL1D::connectCreated);
  return(DUAL1D::subCellToNode.size());
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getNumInterfaceToElement() const
{
  assert(DUAL1D::connectCreated);
  return(DUAL1D::interfaceToElement.size());
}



//_________________________________________________________________________________________________
// GET DUAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getDualPointL(const UInt & lid) const
{
  assert(DUAL1D::connectCreated);
  assert(lid >= 1);
  assert(lid <= DUAL1D::dualPoints.size());
  return(DUAL1D::dualPoints.getL(lid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getDualPointG(const UInt & gid) const
{
  assert(DUAL1D::connectCreated);
  assert(DUAL1D::dualPoints.isG(gid));
  return(DUAL1D::dualPoints.getG(gid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT1D &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellL(const UInt & lid) const
{
  assert(DUAL1D::connectCreated);
  assert(lid >= 1);
  assert(lid <= DUAL1D::subCells.size());
  return(DUAL1D::subCells.getItemL(lid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT1D &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellG(const UInt & gid) const
{
  assert(DUAL1D::connectCreated);
  assert(DUAL1D::subCells.isColG(gid));
  return(DUAL1D::subCells.getItemG(gid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT0D &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceL(const UInt & lid) const
{
  assert(DUAL1D::connectCreated);
  assert(lid >= 1);
  assert(lid <= DUAL1D::interfaces.size());
  return(DUAL1D::interfaces.getItemL(lid));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::GEOELEMENT0D &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceG(const UInt & gid) const
{
  assert(DUAL1D::connectCreated);
  assert(DUAL1D::interfaces.isColG(gid));
  return(DUAL1D::interfaces.getItemG(gid));
}



//_________________________________________________________________________________________________
// GET DUAL CONNECTING FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellToElement(const UInt & i) const
{
  assert(DUAL1D::connectCreated);
  assert(i >= 1);
  assert(i <= DUAL1D::subCellToElement.size());
  return(DUAL1D::subCellToElement(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getSubcellToNode(const UInt & i) const
{
  assert(DUAL1D::connectCreated);
  assert(i >= 1);
  assert(i <= DUAL1D::subcellToNode.size());
  return(DUAL1D::subcellToNode(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const UInt &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceToElement(const UInt & i) const
{
  assert(DUAL1D::connectCreated);
  assert(i >= 1);
  assert(i <= DUAL1D::interfaceToElement.size());
  return(DUAL1D::interfaceToElement(i));
}



//_________________________________________________________________________________________________
// GET FINITE VOLUME FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const Real &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellVolume(const UInt & i) const
{
  assert(DUAL1D::connectCreated);
  assert(DUAL1D::finiteVolumeCreated);
  assert(i >= 1);
  assert(i <= DUAL1D::subCellVolume.size());
  return(DUAL1D::subCellVolume(i));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const point3d &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getCellN(const UInt & i) const
{
  assert(DUAL1D::connectCreated);
  assert(DUAL1D::finiteVolumeCreated);
  assert(i >= 1);
  assert(i <= DUAL1D::subCellN.size());
  return(DUAL1D::subCellN(i));
}



//_________________________________________________________________________________________________
// DUMP INTERNAL ARRAYS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::NODESVECT &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getDualPoints() const
{
  assert(DUAL1D::connectCreated);
  return(DUAL1D::dualPoints);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH1D &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCells() const
{
  assert(DUAL1D::connectCreated);
  return(DUAL1D::subCells);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH0D &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaces() const
{
  assert(DUAL1D::connectCreated);
  return(DUAL1D::interfaces);
}



//_________________________________________________________________________________________________
// DUMP INTERNAL CONNECT ARRAYS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::SUBCELL_TO_ELEMENT &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellToElement() const
{
  assert(DUAL1D::connectCreated);
  return(DUAL1D::subCellToElement);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::SUBCELL_TO_NODE &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellToNode() const
{
  assert(DUAL1D::connectCreated);
  return(DUAL1D::subCellToNode);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
const typename dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::INTERFACE_TO_ELEMENT &
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getInterfaceToElement() const
{
  assert(DUAL1D::connectCreated);
  return(DUAL1D::interfaceToElement);
}



//_________________________________________________________________________________________________
// DUAL MESHES PRINTOUT
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
typename dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::MESH1D
dualMesh1d<GEOSHAPE,ELMAP,NODEMAP>::
getSubCellMesh()
{
  assert(DUAL1D::connectCreated);
  
  MESH1D outMesh;
  outMesh.setNodes(DUAL1D::dualPoints);
  outMesh.setElements(DUAL1D::subCells);
  outMesh.transferMap();
  return(outMesh);
}


#endif
