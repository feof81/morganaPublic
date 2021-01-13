/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FEHY3D_EXTERN_HPP
#define FEHY3D_EXTERN_HPP

#include "feHY3d.hpp"
#include "connect2d.hpp"
#include "connect3d.hpp"


/*! Class for building the \c feCards for hybrid finite elements */
template<typename PMAPTYPE>
class feHY3d_extern
{
  /*! @name Typedefs */ //@{
  public:
    typedef linearTriangle                          GEOSHAPE2D;
    typedef linearTetra                             GEOSHAPE3D;
    typedef feHY3d_card                             FECARD;
    typedef mesh2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>    MESH2D;
    typedef mesh3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>    MESH3D;
    typedef connect2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> CONNECT2D;
    typedef connect3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> CONNECT3D;
    typedef pVect<FECARD,PMAPTYPE>                  FECARDS;
    //@}
    
    /*! @name Interal links */ //@{
  public:
    bool geometryLoaded;
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    Teuchos::RCP<MESH2D>       grid2d;
    Teuchos::RCP<MESH3D>       grid3d;
    Teuchos::RCP<CONNECT2D>    connectGrid2d;
    Teuchos::RCP<CONNECT3D>    connectGrid3d;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    feHY3d_extern();
    feHY3d_extern(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    feHY3d_extern(MESH3D & Mesh3d, CONNECT3D & ConnectGrid3d, MESH2D & Mesh2d, CONNECT2D & ConnectGrid2d);
    void setGeometry(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    void setGeometry(MESH3D & Mesh3d, CONNECT3D & ConnectGrid3d, MESH2D & Mesh2d, CONNECT2D & ConnectGrid2d);
    void setCommDev(const Teuchos::RCP<communicator> & CommDev);
    void setCommDev(communicator & CommDev);
    //@}
    
    /*! @name Other functions */ //@{
  public:
    FECARDS buildFeCards(const set<UInt> & geoIds3d, const set<UInt> & geoIds2d) const;
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET
//-------------------------------------------------------------------------------------------------
template<typename PMAPTYPE>
feHY3d_extern<PMAPTYPE>::
feHY3d_extern()
{
  geometryLoaded = false;
  commDevLoaded  = false;
}

template<typename PMAPTYPE>
feHY3d_extern<PMAPTYPE>::
feHY3d_extern(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Mesh2d;
  grid3d         = Mesh3d;
  connectGrid2d  = ConnectGrid2d;
  connectGrid3d  = ConnectGrid3d;
}

template<typename PMAPTYPE>
feHY3d_extern<PMAPTYPE>::
feHY3d_extern(MESH3D & Mesh3d, CONNECT3D & ConnectGrid3d, MESH2D & Mesh2d, CONNECT2D & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Teuchos::rcp(new MESH2D(Mesh2d));
  grid3d         = Teuchos::rcp(new MESH3D(Mesh3d));
  connectGrid2d  = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
}

template<typename PMAPTYPE>
void
feHY3d_extern<PMAPTYPE>::
setGeometry(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Mesh2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Mesh2d;
  grid3d         = Mesh3d;
  connectGrid2d  = ConnectGrid2d;
  connectGrid3d  = ConnectGrid3d;
}

template<typename PMAPTYPE>
void
feHY3d_extern<PMAPTYPE>::
setGeometry(MESH3D & Mesh3d, CONNECT3D & ConnectGrid3d, MESH2D & Mesh2d, CONNECT2D & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Teuchos::rcp(new MESH2D(Mesh2d));
  grid3d         = Teuchos::rcp(new MESH3D(Mesh3d));
  connectGrid2d  = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
}

template<typename PMAPTYPE>
void
feHY3d_extern<PMAPTYPE>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename PMAPTYPE>
void
feHY3d_extern<PMAPTYPE>::
setCommDev(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}



//_________________________________________________________________________________________________
// OTHER
//-------------------------------------------------------------------------------------------------

template<typename PMAPTYPE>
typename feHY3d_extern<PMAPTYPE>::FECARDS
feHY3d_extern<PMAPTYPE>::
buildFeCards(const set<UInt> & geoIds3d, const set<UInt> & geoIds2d) const
{
  //Asserts
  assert(geometryLoaded);
  assert(commDevLoaded);
  
  //Gen the cards
  UInt flid, geoId3d, geoId2d;
  feHY3d_card card;
  
  FECARDS cards(grid3d->getNumElements());
  cards.setMap(grid3d->getElements().getRowMap());
  
  for(UInt elid=1; elid <= grid3d->getNumElements(); ++elid)
  {
    geoId3d = grid3d->getElementL(elid).getGeoId();
    
    if(geoIds3d.count(geoId3d) == 1)
    {
      for(UInt j=1; j <= 4; ++j)
      {
        flid = connectGrid3d->getElementToFace(elid,j);
	
	if(connectGrid3d->getFaceIsBoundary(flid))
	{
	  geoId2d = grid2d->getElementL(connectGrid3d->getFaceBFace(flid)).getGeoId();
	  card.setActiveFace(j, !(geoIds2d.count(geoId2d) == 1));
	}
	else
	{
	  card.setActiveFace(j,true);
	}      
      }
    }
    else
    {
      for(UInt j=1; j <= 4; ++j)
      { card.setActiveFace(j,false); }
    }

    cards(elid) = card;
  }
  
  return(cards);
}


#endif
