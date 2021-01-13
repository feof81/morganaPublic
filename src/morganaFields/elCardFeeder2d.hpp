/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef ELCARDFEEDER2D_HPP
#define ELCARDFEEDER2D_HPP

#include "Teuchos_RCP.hpp"

#include "mesh2d.hpp"
#include "connect2d.hpp"
#include "elCard2d.hpp"



/*! Transform the grid information into \c elCard2d */
template<typename GEOSHAPE, typename PMAPTYPE>
class elCardFeeder2d
{
  /*! @name Typedefs */ //@{
  public:
    typedef typename   GEOSHAPE::GEOBSHAPE      GEOSHAPE1D;
    typedef geoElement<GEOSHAPE1D>              GEOELEMENT1D;
    
    typedef sVect<point3d>       NODES;
    typedef sVect<PMAPTYPE>      NODES_MAP;
    typedef sVect<GEOELEMENT1D>  EDGES;
    typedef sVect<PMAPTYPE>      EDGES_MAP;
    typedef sVect<bool>          EDGES_ORIENTATION;
    
    typedef mesh2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH2D;
    typedef connect2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT2D;
    //@}
    
    /*! @name Interal links */ //@{
  public:
    bool geometryLoaded;
    Teuchos::RCP<MESH2D>     grid2d;
    Teuchos::RCP<CONNECT2D>  connectGrid2d;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    elCardFeeder2d();
    elCardFeeder2d(const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    elCardFeeder2d(MESH2D & Mesh2d, CONNECT2D & ConnectGrid2d);
    void setGeometry(const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    void setGeometry(MESH2D & Mesh2d, CONNECT2D & ConnectGrid2d);
    elCard2d<GEOSHAPE,PMAPTYPE> getCardLocal(const UInt & lid) const;
    elCard2d<GEOSHAPE,PMAPTYPE> getCardGlobal(const UInt & gid) const;
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename PMAPTYPE>
elCardFeeder2d<GEOSHAPE,PMAPTYPE>::
elCardFeeder2d()
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  geometryLoaded = false;
}

template<typename GEOSHAPE, typename PMAPTYPE>
elCardFeeder2d<GEOSHAPE,PMAPTYPE>::
elCardFeeder2d(const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  grid2d         = Mesh2d;
  connectGrid2d  = ConnectGrid2d;
}

template<typename GEOSHAPE, typename PMAPTYPE>
elCardFeeder2d<GEOSHAPE,PMAPTYPE>::
elCardFeeder2d(MESH2D & Mesh2d, CONNECT2D & ConnectGrid2d)
{
  geometryLoaded = true;
  grid2d         = Teuchos::rcp(new MESH2D(Mesh2d));
  connectGrid2d  = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
}

template<typename GEOSHAPE, typename PMAPTYPE>
void
elCardFeeder2d<GEOSHAPE,PMAPTYPE>::
setGeometry(const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  grid2d         = Mesh2d;
  connectGrid2d  = ConnectGrid2d;
}

template<typename GEOSHAPE, typename PMAPTYPE>
void
elCardFeeder2d<GEOSHAPE,PMAPTYPE>::
setGeometry(MESH2D & Mesh2d, CONNECT2D & ConnectGrid2d)
{
  geometryLoaded = true;
  grid2d         = Teuchos::rcp(new MESH2D(Mesh2d));
  connectGrid2d  = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
}

template<typename GEOSHAPE, typename PMAPTYPE>
elCard2d<GEOSHAPE,PMAPTYPE>
elCardFeeder2d<GEOSHAPE,PMAPTYPE>::
getCardLocal(const UInt & lid) const
{
  //Asserts
  assert(geometryLoaded);
  assert(grid2d->getElements().colIsLocal());
  assert(grid2d->getEdges().colIsLocal());
  
  //Allocate
  UInt nid;
  elCard2d<GEOSHAPE,PMAPTYPE> outCard;
  
  //Nodes
  for(UInt i=1; i <= grid2d->getElementL(lid).size(); ++i)
  {
    nid = grid2d->getElementL(lid).getCid(i);    
    outCard.setNode(i,grid2d->getNodeL(nid));
  }
  
  return(outCard);
}

template<typename GEOSHAPE, typename PMAPTYPE>
elCard2d<GEOSHAPE,PMAPTYPE>
elCardFeeder2d<GEOSHAPE,PMAPTYPE>::
getCardGlobal(const UInt & gid) const
{
  assert(geometryLoaded);
  UInt lid  = grid2d->getElements().getMapG(gid).getLid();
  return(getCardLocal(lid));
}

#endif
