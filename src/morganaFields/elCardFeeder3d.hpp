/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef ELCARDFEEDER3D_HPP
#define ELCARDFEEDER3D_HPP

#include "Teuchos_RCP.hpp"

#include "mesh3d.hpp"
#include "connect3d.hpp"
#include "elCard3d.hpp"


/*! Transform the grid information into \c elCard3d */
template<typename GEOSHAPE, typename PMAPTYPE>
class elCardFeeder3d
{
    /*! @name Typedefs */ //@{
  public:    
    typedef mesh3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH3D;
    typedef connect3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT3D;
    //@}
   
    /*! @name Interal links */ //@{
  public:
    bool geometryLoaded;
    Teuchos::RCP<MESH3D>     grid3d;
    Teuchos::RCP<CONNECT3D>  connectGrid3d;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    elCardFeeder3d();
    elCardFeeder3d(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d);
    elCardFeeder3d(MESH3D & Mesh3d, CONNECT3D & ConnectGrid3d);
    void setGeometry(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d);
    void setGeometry(MESH3D & Mesh3d, CONNECT3D & ConnectGrid3d);
    elCard3d<GEOSHAPE,PMAPTYPE> getCardLocal(const UInt & lid) const;
    elCard3d<GEOSHAPE,PMAPTYPE> getCardGlobal(const UInt & gid) const;
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename PMAPTYPE>
elCardFeeder3d<GEOSHAPE,PMAPTYPE>::
elCardFeeder3d()
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  geometryLoaded = false;
}

template<typename GEOSHAPE, typename PMAPTYPE>
elCardFeeder3d<GEOSHAPE,PMAPTYPE>::
elCardFeeder3d(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d)
{
  geometryLoaded = true;
  grid3d         = Mesh3d;
  connectGrid3d  = ConnectGrid3d;
}

template<typename GEOSHAPE, typename PMAPTYPE>
elCardFeeder3d<GEOSHAPE,PMAPTYPE>::
elCardFeeder3d(MESH3D & Mesh3d, CONNECT3D & ConnectGrid3d)
{
  geometryLoaded = true;
  grid3d         = Teuchos::rcp(new MESH3D(Mesh3d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
}

template<typename GEOSHAPE, typename PMAPTYPE>
void
elCardFeeder3d<GEOSHAPE,PMAPTYPE>::
setGeometry(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d)
{
  geometryLoaded = true;
  grid3d         = Mesh3d;
  connectGrid3d  = ConnectGrid3d;
}

template<typename GEOSHAPE, typename PMAPTYPE>
void
elCardFeeder3d<GEOSHAPE,PMAPTYPE>::
setGeometry(MESH3D & Mesh3d, CONNECT3D & ConnectGrid3d)
{
  geometryLoaded = true;
  grid3d         = Teuchos::rcp(new MESH3D(Mesh3d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
}

template<typename GEOSHAPE, typename PMAPTYPE>
elCard3d<GEOSHAPE,PMAPTYPE>
elCardFeeder3d<GEOSHAPE,PMAPTYPE>::
getCardLocal(const UInt & lid) const
{
  //Asserts
  assert(geometryLoaded);
  
  //Allocate
  UInt nid;
  elCard3d<GEOSHAPE,PMAPTYPE> outCard;
  
  //Nodes
  for(UInt i=1; i <= grid3d->getElementL(lid).size(); ++i)
  {
    nid = grid3d->getElementL(lid).getCid(i);    
    outCard.setNode(i,grid3d->getNodeL(nid));
  }
  
  return(outCard);
}

template<typename GEOSHAPE, typename PMAPTYPE>
elCard3d<GEOSHAPE,PMAPTYPE>
elCardFeeder3d<GEOSHAPE,PMAPTYPE>::
getCardGlobal(const UInt & gid) const
{
  assert(geometryLoaded);
  UInt lid  = grid3d->getElements().getMapG(gid).getLid();
  return(getCardLocal(lid));
}



#endif
