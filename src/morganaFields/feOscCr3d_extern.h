/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FEOSCCR3D_EXTERN_H
#define FEOSCCR3D_EXTERN_H

#include "feOscCr3d.hpp"
#include "feOsc3d_card.h"

#include "connect3d.hpp"


/*! Build the directions */
template<Int R>
class feOscCr3d_dirGen
{ };

template<>
class feOscCr3d_dirGen<1>
{
  public:
    sVect<point3d> listH, listY;
  
  public:
    feOscCr3d_dirGen();
    const sVect<point3d> & getH() const;
    const sVect<point3d> & getY() const;
};


/*! Compute the \c feCards for the Oscillating Crouzeix Raviart finite elements.  */
template<Int R, typename PMAPTYPE>
class feOscCr3d_extern
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
    feOscCr3d_extern();
    feOscCr3d_extern(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d);
    feOscCr3d_extern(const MESH3D & Mesh3d, const CONNECT3D & ConnectGrid3d);
    void setGeometry(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d);
    void setGeometry(const MESH3D & Mesh3d, const CONNECT3D & ConnectGrid3d);
    void setCommDev(const Teuchos::RCP<communicator> & CommDev);
    void setCommDev(communicator & CommDev);
    //@}
    
    /*! @name Other functions */ //@{
  public:
    FECARDS buildFeCards(const Real & K) const;
    FECARDS buildFeCardsConj(const Real & K) const;
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET
//-------------------------------------------------------------------------------------------------
template<Int R, typename PMAPTYPE>
feOscCr3d_extern<R,PMAPTYPE>::
feOscCr3d_extern()
{
  geometryLoaded = false;
  commDevLoaded  = false;
}

template<Int R, typename PMAPTYPE>
feOscCr3d_extern<R,PMAPTYPE>::
feOscCr3d_extern(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid3d         = Mesh3d;
  connectGrid3d  = ConnectGrid3d;
}

template<Int R, typename PMAPTYPE>
feOscCr3d_extern<R,PMAPTYPE>::
feOscCr3d_extern(const MESH3D & Mesh3d, const CONNECT3D & ConnectGrid3d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid3d         = Teuchos::rcp(new MESH3D(Mesh3d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
}

template<Int R, typename PMAPTYPE>
void
feOscCr3d_extern<R,PMAPTYPE>::
setGeometry(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid3d         = Mesh3d;
  connectGrid3d  = ConnectGrid3d;
}

template<Int R, typename PMAPTYPE>
void
feOscCr3d_extern<R,PMAPTYPE>::
setGeometry(const MESH3D & Mesh3d, const CONNECT3D & ConnectGrid3d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid3d         = Teuchos::rcp(new MESH3D(Mesh3d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
}

template<Int R, typename PMAPTYPE>
void
feOscCr3d_extern<R,PMAPTYPE>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<Int R, typename PMAPTYPE>
void
feOscCr3d_extern<R,PMAPTYPE>::
setCommDev(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

//_________________________________________________________________________________________________
// OTHER
//-------------------------------------------------------------------------------------------------
template<Int R, typename PMAPTYPE>
typename feOscCr3d_extern<R,PMAPTYPE>::FECARDS
feOscCr3d_extern<R,PMAPTYPE>::
buildFeCards(const Real & K) const
{
  //Asserts
  assert(geometryLoaded);
  assert(commDevLoaded);
  
  //Download vects
  feOscCr3d_dirGen<1> vectGen;
  sVect<point3d> Href = vectGen.getH();
  sVect<point3d> Yref = vectGen.getY();
  
  //Alloc
  point3d H, Y;
  tensor3d T;
  Real alpha;
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
    T.computeInverse();
    
    for(UInt i=1; i <= Href.size(); ++i)
    {
      Y = Yref(i);
      H = T.firstIndexSaturation(Href(i));
      alpha = K / H.norm2();
      
      feCard.set(Y,Href(i) * alpha,i);
    }
    
    cards(el) = feCard;
  }
  
  return(cards);
}

template<Int R, typename PMAPTYPE>
typename feOscCr3d_extern<R,PMAPTYPE>::FECARDS
feOscCr3d_extern<R,PMAPTYPE>::
buildFeCardsConj(const Real & K) const
{
  //Asserts
  assert(geometryLoaded);
  assert(commDevLoaded);
  
  //Download vects
  feOscCr3d_dirGen<1> vectGen;
  sVect<point3d> Href = vectGen.getH();
  sVect<point3d> Yref = vectGen.getY();
  
  //Alloc
  point3d H, Y;
  tensor3d T;
  Real alpha;
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
    T.computeInverse();
    
    for(UInt i=1; i <= Href.size(); ++i)
    {
      Y = Yref(i);
      H = T.firstIndexSaturation(Href(i));
      alpha = K / H.norm2();
      
      feCard.set(Y,Href(i) * (-alpha),i);
    }
    
    cards(el) = feCard;
  }
  
  return(cards);
}

#endif
