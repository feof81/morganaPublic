/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FEOSCPR2D_EXTERN_H
#define FEOSCPR2D_EXTERN_H

#include "feOscPr2d.hpp"
#include "feOsc3d_card.h"

#include "connect2d.hpp"

/*! Build the directions */
template<Int R, Int NW>
class feOscPr2d_dirGen
{ };

template<>
class feOscPr2d_dirGen<0,0>
{
  public:
    sVect<point3d> listH, listY;
  
  public:
    feOscPr2d_dirGen();
    const sVect<point3d> & getH() const;
    const sVect<point3d> & getY() const;
};

template<>
class feOscPr2d_dirGen<1,0>
{
  public:
    sVect<point3d> listH, listY;
  
  public:
    feOscPr2d_dirGen();
    const sVect<point3d> & getH() const;
    const sVect<point3d> & getY() const;
};

template<>
class feOscPr2d_dirGen<0,1>
{
  public:
    sVect<point3d> listH, listY;
  
  public:
    feOscPr2d_dirGen();
    const sVect<point3d> & getH() const;
    const sVect<point3d> & getY() const;
};

template<>
class feOscPr2d_dirGen<1,1>
{
  public:
    sVect<point3d> listH, listY;
  
  public:
    feOscPr2d_dirGen();
    const sVect<point3d> & getH() const;
    const sVect<point3d> & getY() const;
};


/*! Compute the \c feCards for 3d Raviart Thomas finite elements.  */
template<Int R, Int NW, typename PMAPTYPE>
class feOscPr2d_extern
{
    /*! @name Typedefs */ //@{
  public:
    typedef linearTriangle                        GEOSHAPE;
    typedef feOsc3d_card                          FECARD;
    typedef mesh2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH2D;
    typedef connect2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT2D;
    typedef pVect<FECARD,PMAPTYPE>                FECARDS;
    //@}
    
    /*! @name Interal links and data */ //@{
  public:
    bool geometryLoaded;
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    Teuchos::RCP<MESH2D>       grid2d;
    Teuchos::RCP<CONNECT2D>    connectGrid2d;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    feOscPr2d_extern();
    feOscPr2d_extern(const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    feOscPr2d_extern(const MESH2D & Mesh2d, const CONNECT2D & ConnectGrid2d);
    void setGeometry(const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    void setGeometry(const MESH2D & Mesh2d, const CONNECT2D & ConnectGrid2d);
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
feOscPr2d_extern<R,NW,PMAPTYPE>::
feOscPr2d_extern()
{
  geometryLoaded = false;
  commDevLoaded  = false;
}

template<Int R, Int NW, typename PMAPTYPE>
feOscPr2d_extern<R,NW,PMAPTYPE>::
feOscPr2d_extern(const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Mesh2d;
  connectGrid2d  = ConnectGrid2d;
}

template<Int R, Int NW, typename PMAPTYPE>
feOscPr2d_extern<R,NW,PMAPTYPE>::
feOscPr2d_extern(const MESH2D & Mesh2d, const CONNECT2D & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Teuchos::rcp(new MESH2D(Mesh2d));
  connectGrid2d  = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
}

template<Int R, Int NW, typename PMAPTYPE>
void
feOscPr2d_extern<R,NW,PMAPTYPE>::
setGeometry(const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Mesh2d;
  connectGrid2d  = ConnectGrid2d;
}

template<Int R, Int NW, typename PMAPTYPE>
void
feOscPr2d_extern<R,NW,PMAPTYPE>::
setGeometry(const MESH2D & Mesh2d, const CONNECT2D & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Teuchos::rcp(new MESH2D(Mesh2d));
  connectGrid2d  = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
}

template<Int R, Int NW, typename PMAPTYPE>
void
feOscPr2d_extern<R,NW,PMAPTYPE>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<Int R, Int NW, typename PMAPTYPE>
void
feOscPr2d_extern<R,NW,PMAPTYPE>::
setCommDev(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}


//_________________________________________________________________________________________________
// OTHER
//-------------------------------------------------------------------------------------------------
template<Int R, Int NW, typename PMAPTYPE>
typename feOscPr2d_extern<R,NW,PMAPTYPE>::FECARDS
feOscPr2d_extern<R,NW,PMAPTYPE>::
buildFeCards(const Real & K) const
{
  //Asserts
  assert(geometryLoaded);
  assert(commDevLoaded);
  
  //Download vects
  feOscPr2d_dirGen<R,NW> vectGen;
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
  FECARDS cards(grid2d->getNumElements());
  cards.setMap(grid2d->getElements().getRowMap());
  
  for(UInt el=1; el <= grid2d->getNumElements(); ++el)
  {
    nodes = grid2d->getElementNodesL(el);
    
    T.setCol(1,nodes(2) - nodes(1));
    T.setCol(2,nodes(3) - nodes(1));
    T.completeThirdColoumn();
    
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
