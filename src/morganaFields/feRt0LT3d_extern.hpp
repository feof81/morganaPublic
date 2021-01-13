/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FERT0LT3D_EXTERN_HPP
#define FERT0LT3D_EXTERN_HPP

#include "feRt0LT3d.hpp"
#include "feRt0LT3d_card.h"

#include "connect3d.hpp"


/*! Compute the \c feCards for 3d Raviart Thomas finite elements.  */
template<typename PMAPTYPE>
class feRt0LT3d_extern
{
    /*! @name Typedefs */ //@{
  public:
    typedef linearTetra                           GEOSHAPE;
    typedef feRt0LT3d_card                        FECARD;
    typedef mesh3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH3D;
    typedef connect3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT3D;
    typedef pVect<FECARD,PMAPTYPE>                FECARDS;
    //@}
    
    /*! @name Interal links */ //@{
  public:
    bool geometryLoaded;
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    Teuchos::RCP<MESH3D>       grid3d;
    Teuchos::RCP<CONNECT3D>    connectGrid3d;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    feRt0LT3d_extern();
    feRt0LT3d_extern(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d);
    feRt0LT3d_extern(const MESH3D & Mesh3d, const CONNECT3D & ConnectGrid3d);
    void setGeometry(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d);
    void setGeometry(const MESH3D & Mesh3d, const CONNECT3D & ConnectGrid3d);
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
feRt0LT3d_extern<PMAPTYPE>::
feRt0LT3d_extern()
{
  geometryLoaded = false;
  commDevLoaded  = false;
}

template<typename PMAPTYPE>
feRt0LT3d_extern<PMAPTYPE>::
feRt0LT3d_extern(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid3d         = Mesh3d;
  connectGrid3d  = ConnectGrid3d;
}

template<typename PMAPTYPE>
feRt0LT3d_extern<PMAPTYPE>::
feRt0LT3d_extern(const MESH3D & Mesh3d, const CONNECT3D & ConnectGrid3d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid3d         = Teuchos::rcp(new MESH3D(Mesh3d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
}

template<typename PMAPTYPE>
void
feRt0LT3d_extern<PMAPTYPE>::
setGeometry(const Teuchos::RCP<MESH3D> & Mesh3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid3d         = Mesh3d;
  connectGrid3d  = ConnectGrid3d;
}

template<typename PMAPTYPE>
void
feRt0LT3d_extern<PMAPTYPE>::
setGeometry(const MESH3D & Mesh3d, const CONNECT3D & ConnectGrid3d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid3d         = Teuchos::rcp(new MESH3D(Mesh3d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
}

template<typename PMAPTYPE>
void
feRt0LT3d_extern<PMAPTYPE>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename PMAPTYPE>
void
feRt0LT3d_extern<PMAPTYPE>::
setCommDev(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}



//_________________________________________________________________________________________________
// OTHER
//-------------------------------------------------------------------------------------------------
template<typename PMAPTYPE>
typename feRt0LT3d_extern<PMAPTYPE>::FECARDS
feRt0LT3d_extern<PMAPTYPE>::
buildFeCards() const
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
  FECARDS cards(grid3d->getNumElements());
  cards.setMap(grid3d->getElements().getRowMap());
  
  bool flag;
  pGraphItem item;
  UInt egid, flid, fgid, emin;
  feRt0LT3d_card card;
  
  for(UInt elid=1; elid <= grid3d->getNumElements(); ++elid)
  {
    egid = grid3d->getElements().getRowMapL(elid).getGid();
    
    for(UInt j=1; j <= 4; ++j)
    {
      flid = connectGrid3d->getElementToFace(elid,j);
      fgid = grid3d->getFaces().getMapL(flid).getGid();
      
      item = faceToElement.getItemG(fgid);
      
      emin = item.getCid(1);
      if(item.size() > 1) { emin = std::min(emin, item.getCid(2)); }
      assert(item.size() <= 2);
      
      flag = (emin == egid);
     
      card.setFaceOrientation(j,flag);      
    }
    
    cards(elid) = card;
  }
  
  return(cards);
}



#endif
