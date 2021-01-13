/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FEOSCPR3D_EXTERN_H
#define FEOSCPR3D_EXTERN_H

#include "feOscPr3d.hpp"
#include "feOsc3d_card.h"

#include "connect3d.hpp"

/*! Build the directions */
template<Int R, Int NW>
class feOscPr3d_dirGen
{ };

template<>
class feOscPr3d_dirGen<0,1>
{
  public:
    sVect<point3d> listH, listY;
  
  public:
    feOscPr3d_dirGen();
    const sVect<point3d> & getH() const;
    const sVect<point3d> & getY() const;
};

template<>
class feOscPr3d_dirGen<1,1>
{
  public:
    sVect<point3d> listH, listY;
  
  public:
    feOscPr3d_dirGen();
    const sVect<point3d> & getH() const;
    const sVect<point3d> & getY() const;
};


/*! Compute the \c feCards for the Oscillating Pr finite elements.  */
template<Int R, Int NW, typename PMAPTYPE>
class feOscPr3d_extern
{
    /*! @name Typedefs */ //@{
  public:
    typedef linearTetra                           GEOSHAPE;
    typedef feOsc3d_card                          FECARD;
    typedef mesh3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH3D;
    typedef connect3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT3D;
    typedef pVect<FECARD,PMAPTYPE>                FECARDS;
    //@}
    
    /*! @name Interal links and data */ //@{
  public:
    bool geometryLoaded;
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    Teuchos::RCP<MESH3D>       grid3d;
    Teuchos::RCP<CONNECT3D>    connectGrid3d;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    feOscPr3d_extern();
    feOscPr3d_extern(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d);
    feOscPr3d_extern(const MESH3D & Mesh3d, const CONNECT3D & ConnectGrid3d);
    void setGeometry(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d);
    void setGeometry(const MESH3D & Mesh3d, const CONNECT3D & ConnectGrid3d);
    void setCommDev(const Teuchos::RCP<communicator> & CommDev);
    void setCommDev(communicator & CommDev);
    //@}
    
    /*! @name Other functions */ //@{
  public:
    FECARDS buildFeCards(const Real & K) const;
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET
//-------------------------------------------------------------------------------------------------
template<Int R, Int NW, typename PMAPTYPE>
feOscPr3d_extern<R,NW,PMAPTYPE>::
feOscPr3d_extern()
{
  geometryLoaded = false;
  commDevLoaded  = false;
}

template<Int R, Int NW, typename PMAPTYPE>
feOscPr3d_extern<R,NW,PMAPTYPE>::
feOscPr3d_extern(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid3d         = Mesh3d;
  connectGrid3d  = ConnectGrid3d;
}

template<Int R, Int NW, typename PMAPTYPE>
feOscPr3d_extern<R,NW,PMAPTYPE>::
feOscPr3d_extern(const MESH3D & Mesh3d, const CONNECT3D & ConnectGrid3d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid3d         = Teuchos::rcp(new MESH3D(Mesh3d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
}

template<Int R, Int NW, typename PMAPTYPE>
void
feOscPr3d_extern<R,NW,PMAPTYPE>::
setGeometry(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid3d         = Mesh3d;
  connectGrid3d  = ConnectGrid3d;
}

template<Int R, Int NW, typename PMAPTYPE>
void
feOscPr3d_extern<R,NW,PMAPTYPE>::
setGeometry(const MESH3D & Mesh3d, const CONNECT3D & ConnectGrid3d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid3d         = Teuchos::rcp(new MESH3D(Mesh3d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
}

template<Int R, Int NW, typename PMAPTYPE>
void
feOscPr3d_extern<R,NW,PMAPTYPE>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<Int R, Int NW, typename PMAPTYPE>
void
feOscPr3d_extern<R,NW,PMAPTYPE>::
setCommDev(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}


//_________________________________________________________________________________________________
// OTHER
//-------------------------------------------------------------------------------------------------
template<Int R, Int NW, typename PMAPTYPE>
typename feOscPr3d_extern<R,NW,PMAPTYPE>::FECARDS
feOscPr3d_extern<R,NW,PMAPTYPE>::
buildFeCards(const Real & K) const
{
  //Asserts
  assert(geometryLoaded);
  assert(commDevLoaded);
  
  //Download vects
  feOscPr3d_dirGen<R,NW> vectGen;
  sVect<point3d> Href = vectGen.getH();
  sVect<point3d> Yref = vectGen.getY();
  
  for(UInt i=1; i <= Href.size(); ++i)
  { Href(i) *= K; }
  
  //Alloc
  point3d H, Y;
  tensor3d T;
  sVect<point3d> nodes;
  FECARD feCard(Href.size());
  
  //Gen the cards
  FECARDS cards(grid3d->getNumElements());
  cards.setMap(grid3d->getElements().getRowMap());
  
  for(UInt el=1; el <= grid3d->getNumElements(); ++el)
  {
    nodes = grid3d->getElementNodesL(el);
    
    T.setCol(1,nodes(2) - nodes(1));
    T.setCol(2,nodes(3) - nodes(1));
    T.setCol(3,nodes(4) - nodes(1));
    
    for(UInt i=1; i <= Href.size(); ++i)
    {
      Y = Yref(i);
      H = T.firstIndexSaturation(Href(i));
      
      feCard.set(Y,H,i);
    }
    
    cards(el) = feCard;
  }
  
  return(cards);
}

#endif
