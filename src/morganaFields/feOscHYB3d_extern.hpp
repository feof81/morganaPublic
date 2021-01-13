/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FEOSCHYB3D_EXTERN_HPP
#define FEOSCHYB3D_EXTERN_HPP

#include "feHYB3d.hpp"
#include "connect2d.hpp"
#include "connect3d.hpp"


/*! Class for building the \c feCards for hybrid finite elements */
template<typename PMAPTYPE>
class feOscHYB3d_extern
{
  /*! @name Typedefs */ //@{
  public:
    typedef linearTriangle                          GEOSHAPE2D;
    typedef linearTetra                             GEOSHAPE3D;
    typedef feOscHYB3d_card                         FECARD;
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
    feOscHYB3d_extern();
    feOscHYB3d_extern(const Teuchos::RCP<MESH3D>    & Mesh3d,
                      const Teuchos::RCP<CONNECT3D> & ConnectGrid3d,
                      const Teuchos::RCP<MESH2D>    & Mesh2d,
                      const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    
    feOscHYB3d_extern(MESH3D    & Mesh3d,
                      CONNECT3D & ConnectGrid3d,
                      MESH2D    & Mesh2d,
                      CONNECT2D & ConnectGrid2d);
    
    void setGeometry(const Teuchos::RCP<MESH3D>    & Mesh3d,
                     const Teuchos::RCP<CONNECT3D> & ConnectGrid3d,
                     const Teuchos::RCP<MESH2D>    & Mesh2d,
                     const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    
    void setGeometry(MESH3D    & Mesh3d,
                     CONNECT3D & ConnectGrid3d,
                     MESH2D    & Mesh2d,
                     CONNECT2D & ConnectGrid2d);
    
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
feOscHYB3d_extern<PMAPTYPE>::
feOscHYB3d_extern()
{
  geometryLoaded = false;
  commDevLoaded  = false;
}

template<typename PMAPTYPE>
feOscHYB3d_extern<PMAPTYPE>::
feOscHYB3d_extern(const Teuchos::RCP<MESH3D>    & Mesh3d,
	       const Teuchos::RCP<CONNECT3D> & ConnectGrid3d,
	       const Teuchos::RCP<MESH2D>    & Mesh2d,
	       const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Mesh2d;
  grid3d         = Mesh3d;
  connectGrid2d  = ConnectGrid2d;
  connectGrid3d  = ConnectGrid3d;
}

template<typename PMAPTYPE>
feOscHYB3d_extern<PMAPTYPE>::
feOscHYB3d_extern(MESH3D    & Mesh3d,
	       CONNECT3D & ConnectGrid3d,
	       MESH2D    & Mesh2d,
	       CONNECT2D & ConnectGrid2d)
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
feOscHYB3d_extern<PMAPTYPE>::
setGeometry(const Teuchos::RCP<MESH3D>    & Mesh3d,
	    const Teuchos::RCP<CONNECT3D> & ConnectGrid3d,
	    const Teuchos::RCP<MESH2D>    & Mesh2d,
	    const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
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
feOscHYB3d_extern<PMAPTYPE>::
setGeometry(MESH3D    & Mesh3d,
	    CONNECT3D & ConnectGrid3d,
	    MESH2D    & Mesh2d,
	    CONNECT2D & ConnectGrid2d)
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
feOscHYB3d_extern<PMAPTYPE>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename PMAPTYPE>
void
feOscHYB3d_extern<PMAPTYPE>::
setCommDev(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}



//_________________________________________________________________________________________________
// OTHER
//-------------------------------------------------------------------------------------------------
template<typename PMAPTYPE>
typename feOscHYB3d_extern<PMAPTYPE>::FECARDS
feOscHYB3d_extern<PMAPTYPE>::
buildFeCards(const set<UInt> & geoIds3d, const set<UInt> & geoIds2d) const
{
  //Asserts
  assert(geometryLoaded);
  assert(commDevLoaded);
  
   //Typedefs
  typedef typename CONNECT3D::ELEMENT_TO_FACE             ELEMENT_TO_FACE;
  typedef pGraphGlobalManip<pGraphItem,PMAPTYPE,PMAPTYPE> GRAPH_MANIP;
  
  //Inversion
  ELEMENT_TO_FACE elementToFace = connectGrid3d->getElementToFace();
  ELEMENT_TO_FACE faceToElement;
  
  GRAPH_MANIP manipulator(commDev);
  manipulator.inversion(elementToFace,faceToElement);
  
  faceToElement.pushToGlobal();
  
  //Gen the cards
  point3d H(0.0,0.0,0.0), Y(1.0/3.0, 1.0/3.0, 1.0/3.0);
  
  pGraphItem item;
  UInt flid, egid, fgid, emin, geoId3d, geoId2d;
  FECARD card(4);
  
  FECARDS cards(grid3d->getNumElements());
  cards.setMap(grid3d->getElements().getRowMap());
  
  for(UInt elid=1; elid <= grid3d->getNumElements(); ++elid)
  {
    geoId3d = grid3d->getElementL(elid).getGeoId();
    egid    = grid3d->getElements().getRowMapL(elid).getGid();
    
    if(geoIds3d.count(geoId3d) == 1)
    {
      for(UInt j=1; j <= 4; ++j)
      {
        flid = connectGrid3d->getElementToFace(elid,j);
	fgid = grid3d->getFaces().getMapL(flid).getGid();
	item = faceToElement.getItemG(fgid);
	
	emin = item.getCid(1);
        if(item.size() > 1) { emin = std::min(emin, item.getCid(2)); }
        assert(item.size() <= 2);
      
	card.setFaceOrientation(j,emin == egid);
	
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
      {
	flid = connectGrid3d->getElementToFace(elid,j);
	fgid = grid3d->getFaces().getMapL(flid).getGid();
	item = faceToElement.getItemG(fgid);
	
	emin = item.getCid(1);
        if(item.size() > 1) { emin = std::min(emin, item.getCid(2)); }
        assert(item.size() <= 2);
      
	card.setFaceOrientation(j,emin == egid);
	card.setActiveFace(j,false);
      }
    }
    
    card.set(Y,H,1);
    card.set(Y,H,2);
    card.set(Y,H,3);
    card.set(Y,H,4);
    
    cards(elid) = card;
  }
  
  return(cards);
}

#endif
