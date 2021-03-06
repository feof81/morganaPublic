/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FEHY2D_EXTERN_HPP
#define FEHY2D_EXTERN_HPP

#include "Teuchos_RCP.hpp"

#include "typesInterface.hpp"
#include "feHY2d.hpp"


/*! Class for building the \c feCards for 2d hybrid finite elements */
template<typename PMAPTYPE>
class feHY2d_extern
{
  /*! @name Typedefs */ //@{
  public:
    typedef linearTriangle                          GEOSHAPE2D;    
    typedef feHY2d_card                             FECARD;
    typedef mesh2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>    MESH2D;
    typedef connect2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> CONNECT2D;    
    typedef pVect<FECARD,PMAPTYPE>                  FECARDS;
    //@}
    
    /*! @name Interal links */ //@{
  public:
    bool geometryLoaded;
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    Teuchos::RCP<MESH2D>       grid2d;
    Teuchos::RCP<CONNECT2D>    connectGrid2d;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    feHY2d_extern();
    feHY2d_extern(const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    feHY2d_extern(MESH2D & Mesh2d, CONNECT2D & ConnectGrid2d);
    void setGeometry(const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    void setGeometry(MESH2D & Mesh2d, CONNECT2D & ConnectGrid2d);
    void setCommDev(const Teuchos::RCP<communicator> & CommDev);
    void setCommDev(communicator & CommDev);
    //@}
    
    /*! @name Other functions */ //@{
  public:
    FECARDS buildFeCards() const;
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET
//-------------------------------------------------------------------------------------------------
template<typename PMAPTYPE>
feHY2d_extern<PMAPTYPE>::
feHY2d_extern()
{
  geometryLoaded = false;
  commDevLoaded  = false;
}

template<typename PMAPTYPE>
feHY2d_extern<PMAPTYPE>::
feHY2d_extern(const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Mesh2d;
  connectGrid2d  = ConnectGrid2d;
}

template<typename PMAPTYPE>
feHY2d_extern<PMAPTYPE>::
feHY2d_extern(MESH2D & Mesh2d, CONNECT2D & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Teuchos::rcp(new MESH2D(Mesh2d));
  connectGrid2d  = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
}

template<typename PMAPTYPE>
void
feHY2d_extern<PMAPTYPE>::
setGeometry(const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Mesh2d;
  connectGrid2d  = ConnectGrid2d;
}

template<typename PMAPTYPE>
void
feHY2d_extern<PMAPTYPE>::
setGeometry(MESH2D & Mesh2d, CONNECT2D & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Teuchos::rcp(new MESH2D(Mesh2d));
  connectGrid2d  = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
}

template<typename PMAPTYPE>
void
feHY2d_extern<PMAPTYPE>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename PMAPTYPE>
void
feHY2d_extern<PMAPTYPE>::
setCommDev(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}



//_________________________________________________________________________________________________
// OTHER
//-------------------------------------------------------------------------------------------------
template<typename PMAPTYPE>
typename feHY2d_extern<PMAPTYPE>::FECARDS
feHY2d_extern<PMAPTYPE>::
buildFeCards() const
{
  //Asserts
  assert(geometryLoaded);
  assert(commDevLoaded);
  
  //Gen the cards
  UInt flid;
  feHY2d_card card;
  
  FECARDS cards(grid2d->getNumElements());
  cards.setMap(grid2d->getElements().getRowMap());
  
  for(UInt elid=1; elid <= grid2d->getNumElements(); ++elid)
  {
    for(UInt j=1; j <= 3; ++j)
    {
      flid = connectGrid2d->getElementToEdge(elid,j);
      card.setActiveEdge(j, !(connectGrid2d->getEdgeIsBoundary(flid)) );
    }
    
    cards(elid) = card;
  }
  
  return(cards);
}


#endif
