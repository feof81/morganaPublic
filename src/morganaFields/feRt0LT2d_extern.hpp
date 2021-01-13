/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FERT0LT2D_EXTERN_HPP
#define FERT0LT2D_EXTERN_HPP

#include "feRt0LT2d.hpp"
#include "feRt0LT2d_card.h"
#include "pMapItem.h"


/*! Produce the \c feCards for the 2d Raviarth - Thomas finite elements */
template<typename PMAPTYPE>
class feRt0LT2d_extern
{
    /*! @name Typedefs */ //@{
  public:
    typedef linearTriangle                        GEOSHAPE;
    typedef feRt0LT2d_card                        FECARD;
    typedef mesh2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH2D;
    typedef connect2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT2D;
    typedef pVect<FECARD,PMAPTYPE>                FECARDS;
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
    feRt0LT2d_extern();
    feRt0LT2d_extern(const Teuchos::RCP<MESH2D>    & Mesh2d,
                     const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    
    feRt0LT2d_extern(const MESH2D    & Mesh2d,
                     const CONNECT2D & ConnectGrid2d);
    
    void setGeometry(const Teuchos::RCP<MESH2D>    & Mesh2d,
                     const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    
    void setGeometry(const MESH2D    & Mesh2d,
                     const CONNECT2D & ConnectGrid2d);
    
    void setCommDev(const Teuchos::RCP<communicator> & CommDev);
    void setCommDev(communicator & CommDev);
    //@}
    
    /*! @name Other functions */ //@{
  public:
    FECARDS buildFeCards() const;
    FECARDS buildFeCards(const std::set<UInt> & activeGeoIds) const;
    FECARDS buildFeCards(const std::set<UInt> & activeGeoIds,
                         const std::set<UInt> & majorGeoIds) const;
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET
//-------------------------------------------------------------------------------------------------
template<typename PMAPTYPE>
feRt0LT2d_extern<PMAPTYPE>::
feRt0LT2d_extern()
{
  geometryLoaded = false;
  commDevLoaded  = false;
}

template<typename PMAPTYPE>
feRt0LT2d_extern<PMAPTYPE>::
feRt0LT2d_extern(const Teuchos::RCP<MESH2D>    & Mesh2d,
                 const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Mesh2d;
  connectGrid2d  = ConnectGrid2d;
}

template<typename PMAPTYPE>
feRt0LT2d_extern<PMAPTYPE>::
feRt0LT2d_extern(const MESH2D    & Mesh2d,
                 const CONNECT2D & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Teuchos::rcp(new MESH2D(Mesh2d));
  connectGrid2d  = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
}

template<typename PMAPTYPE>
void
feRt0LT2d_extern<PMAPTYPE>::
setGeometry(const Teuchos::RCP<MESH2D>    & Mesh2d,
            const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Mesh2d;
  connectGrid2d  = ConnectGrid2d;
}

template<typename PMAPTYPE>
void
feRt0LT2d_extern<PMAPTYPE>::
setGeometry(const MESH2D    & Mesh2d,
            const CONNECT2D & ConnectGrid2d)
{
  geometryLoaded = true;
  commDevLoaded  = false;
  
  grid2d         = Teuchos::rcp(new MESH2D(Mesh2d));
  connectGrid2d  = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
}

template<typename PMAPTYPE>
void
feRt0LT2d_extern<PMAPTYPE>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename PMAPTYPE>
void
feRt0LT2d_extern<PMAPTYPE>::
setCommDev(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}



//_________________________________________________________________________________________________
// OTHER
//-------------------------------------------------------------------------------------------------
template<typename PMAPTYPE>
typename feRt0LT2d_extern<PMAPTYPE>::FECARDS
feRt0LT2d_extern<PMAPTYPE>::
buildFeCards() const
{
  //Asserts----------------------------------------------------------
  assert(geometryLoaded);
  assert(commDevLoaded);
  
  //Typedefs---------------------------------------------------------
  typedef typename CONNECT2D::ELEMENT_TO_EDGE             ELEMENT_TO_EDGE;
  typedef pGraphGlobalManip<pGraphItem,PMAPTYPE,PMAPTYPE> GRAPH_MANIP;
  
  //Inversion--------------------------------------------------------
  ELEMENT_TO_EDGE elementToEdge = connectGrid2d->getElementToEdge();
  ELEMENT_TO_EDGE edgeToElement;
  
  GRAPH_MANIP manipulator(commDev);
  manipulator.inversion(elementToEdge,edgeToElement);
  
  edgeToElement.pushToGlobal();
  
  //Gen the cards----------------------------------------------------
  FECARDS cards(grid2d->getNumElements());
  cards.setMap(grid2d->getElements().getRowMap());
  
  bool flag;
  pGraphItem item;
  UInt egid, flid, fgid, emin;
  feRt0LT2d_card card;
  
  for(UInt elid=1; elid <= grid2d->getNumElements(); ++elid)
  {
    egid = grid2d->getElements().getRowMapL(elid).getGid();
    
    for(UInt j=1; j <= 3; ++j)
    {
      flid = connectGrid2d->getElementToEdge(elid,j);
      fgid = grid2d->getEdges().getMapL(flid).getGid();
      
      item = edgeToElement.getItemG(fgid);
      
      emin = item.getCid(1);
      if(item.size() > 1) { emin = std::min(emin, item.getCid(2)); }
      assert(item.size() <= 2);
      
      flag = (emin == egid);
     
      card.setEdgeOrientation(j,flag);      
    }
    
    cards(elid) = card;
  }
  
  return(cards);
}

template<typename PMAPTYPE>
typename feRt0LT2d_extern<PMAPTYPE>::FECARDS
feRt0LT2d_extern<PMAPTYPE>::
buildFeCards(const std::set<UInt> & activeGeoIds) const
{
  //Asserts----------------------------------------------------------
  assert(geometryLoaded);
  assert(commDevLoaded);
  
  //Typedefs---------------------------------------------------------
  typedef typename CONNECT2D::ELEMENT_TO_EDGE             ELEMENT_TO_EDGE;
  typedef pGraphGlobalManip<pGraphItem,PMAPTYPE,PMAPTYPE> GRAPH_MANIP;
  
  //Inversion
  ELEMENT_TO_EDGE elementToEdge = connectGrid2d->getElementToEdge();
  ELEMENT_TO_EDGE edgeToElement;
  
  GRAPH_MANIP manipulator(commDev);
  manipulator.inversion(elementToEdge,edgeToElement);
  
  edgeToElement.pushToGlobal();
  
  //GeoIds needed----------------------------------------------------
  typedef std::set<UInt> SET;
  typedef std::set<UInt>::iterator ITER;
  
  SET elGids;
  ITER iter;
  
  UInt lid = 1;
  pMapItem mapItem;
  pMap<pMapItem> newMap;
  pVect<UInt,pMapItem> vectGeoIds;
  pGraphItem item, tempItem;
  UInt egid, flid, fgid, emin, geoId;
  
  for(UInt i=1; i <= edgeToElement.size(); ++i)
  {
    tempItem = edgeToElement.getItemL(i);
    
    for(UInt k=1; k <= tempItem.size(); ++k)
    { elGids.insert(tempItem.getCid(k)); }
  }
  
  for(ITER iter = elGids.begin(); iter != elGids.end(); iter++)
  {
    newMap.push_back(pMapItem(lid,*iter,commDev->rank()));
    ++lid;
  }
  
  for(UInt i=1; i <= grid2d->getNumElements(); ++i)
  {
    mapItem.setPid(grid2d->getElements().getMapL(i).getPid());
    mapItem.setLid(grid2d->getElements().getMapL(i).getLid());
    mapItem.setGid(grid2d->getElements().getMapL(i).getGid());
    
    geoId = grid2d->getElementL(i).getGeoId();
    
    vectGeoIds.push_back(geoId,mapItem);
  }
  
  vectGeoIds.updateFinder();
  
  pVectGlobalManip<UInt,pMapItem> globalManip(commDev);
  globalManip.changeMap(vectGeoIds,newMap);
  
  //Gen the cards----------------------------------------------------
  FECARDS cards(grid2d->getNumElements());
  cards.setMap(grid2d->getElements().getRowMap());
  
  bool flag;
  feRt0LT2d_card card;
  
  for(UInt elid=1; elid <= grid2d->getNumElements(); ++elid)
  {
    egid  = grid2d->getElements().getRowMapL(elid).getGid();
    geoId = grid2d->getElementL(elid).getGeoId();
    
    if(activeGeoIds.count(geoId) > 0)
    {
      for(UInt j=1; j <= 3; ++j)
      {
        flid = connectGrid2d->getElementToEdge(elid,j);
        fgid = grid2d->getEdges().getMapL(flid).getGid();

        tempItem = edgeToElement.getItemG(fgid);
        item.clear();

        for(UInt k=1; k <= tempItem.size(); ++k)
        {
          geoId = vectGeoIds.getG(tempItem.getCid(k));
          if(activeGeoIds.count(geoId) > 0) { item.push_back(tempItem.getCid(k)); }
        }

        assert(item.size() >= 1); //At least one must be present
        assert(item.size() <= 2); //No more two adjactent elements are allowed

        emin = item.getCid(1);
        if(item.size() > 1) { emin = std::min(emin, item.getCid(2)); }

        flag = (emin == egid);
      
        card.setEdgeOrientation(j,flag);
      }
    }
    else
    {
      for(UInt j=1; j <= 3; ++j)
      { card.setEdgeOrientation(j,true); }
    }
    
    cards(elid) = card;
  }
  
  return(cards);
}

template<typename PMAPTYPE>
typename feRt0LT2d_extern<PMAPTYPE>::FECARDS
feRt0LT2d_extern<PMAPTYPE>::
buildFeCards(const std::set<UInt> & activeGeoIds,
             const std::set<UInt> & majorGeoIds) const
{
  //Asserts----------------------------------------------------------
  assert(geometryLoaded);
  assert(commDevLoaded);
  
  //Typedefs---------------------------------------------------------
  typedef typename CONNECT2D::ELEMENT_TO_EDGE             ELEMENT_TO_EDGE;
  typedef pGraphGlobalManip<pGraphItem,PMAPTYPE,PMAPTYPE> GRAPH_MANIP;
  
  //Inversion
  ELEMENT_TO_EDGE elementToEdge = connectGrid2d->getElementToEdge();
  ELEMENT_TO_EDGE edgeToElement;
  
  GRAPH_MANIP manipulator(commDev);
  manipulator.inversion(elementToEdge,edgeToElement);
  
  edgeToElement.pushToGlobal();
  
  //Vector of global geoIds------------------------------------------
  pVect<UInt,PMAPTYPE> geoIdVect;
  
  for(UInt el=1; el <= grid2d->getNumElements(); ++el)
  {
    geoIdVect.push_back(grid2d->getElements().getItemL(el).getGeoId(),
                        grid2d->getElements().getRowMapL(el));
  }
  
  geoIdVect.updateFinder();
  
  pVectGlobalManip<UInt,PMAPTYPE> geoIdManip(commDev);
  geoIdManip.changeMap(geoIdVect,
                       edgeToElement.getColMap());
  
  //GeoIds needed----------------------------------------------------
  typedef std::set<UInt> SET;
  typedef std::set<UInt>::iterator ITER;
  
  SET elGids;
  ITER iter;
  
  UInt lid = 1;
  pMapItem mapItem;
  pMap<pMapItem> newMap;
  pVect<UInt,pMapItem> vectGeoIds;
  pGraphItem item, tempItem;
  UInt egid, flid, fgid, emin, geoId, geoIdEl;
  
  for(UInt i=1; i <= edgeToElement.size(); ++i)
  {
    tempItem = edgeToElement.getItemL(i);
    
    for(UInt k=1; k <= tempItem.size(); ++k)
    { elGids.insert(tempItem.getCid(k)); }
  }
  
  for(ITER iter = elGids.begin(); iter != elGids.end(); iter++)
  {
    newMap.push_back(pMapItem(lid,*iter,commDev->rank()));
    ++lid;
  }
  
  for(UInt i=1; i <= grid2d->getNumElements(); ++i)
  {
    mapItem.setPid(grid2d->getElements().getMapL(i).getPid());
    mapItem.setLid(grid2d->getElements().getMapL(i).getLid());
    mapItem.setGid(grid2d->getElements().getMapL(i).getGid());
    
    geoId = grid2d->getElementL(i).getGeoId();
    
    vectGeoIds.push_back(geoId,mapItem);
  }
  
  vectGeoIds.updateFinder();
  
  pVectGlobalManip<UInt,pMapItem> globalManip(commDev);
  globalManip.changeMap(vectGeoIds,newMap);
  
  //Gen the cards----------------------------------------------------
  FECARDS cards(grid2d->getNumElements());
  cards.setMap(grid2d->getElements().getRowMap());
  
  bool flag;
  feRt0LT2d_card card;
  
  for(UInt elid=1; elid <= grid2d->getNumElements(); ++elid)
  {
    egid    = grid2d->getElements().getRowMapL(elid).getGid();
    geoIdEl = grid2d->getElementL(elid).getGeoId();
    
    if(activeGeoIds.count(geoIdEl) > 0)
    {
      for(UInt j=1; j <= 3; ++j)
      {
        flid = connectGrid2d->getElementToEdge(elid,j);
        fgid = grid2d->getEdges().getMapL(flid).getGid();

        tempItem = edgeToElement.getItemG(fgid);
        item.clear();

        for(UInt k=1; k <= tempItem.size(); ++k)
        {
          geoId = vectGeoIds.getG(tempItem.getCid(k));
          if(activeGeoIds.count(geoId) > 0) { item.push_back(tempItem.getCid(k)); }
        }
        
        if(item.size() <= 2)
        {
          emin = item.getCid(1);
  
          if(item.size() == 2) 
          { emin = std::min(emin, item.getCid(2)); }
  
          flag = (emin == egid);
        }
        else
        {
          flag = (majorGeoIds.count(geoIdEl) == 1);
  
          for(UInt k=1; k <= item.size(); ++k)
          {
            if((item.getCid(k) != egid) && flag)
            { assert(majorGeoIds.count(geoIdVect.getDataG(item.getCid(k))) != 1); }
          }
        }

        card.setEdgeOrientation(j,flag);
      }
    }
    else
    {
      for(UInt j=1; j <= 3; ++j)
      { card.setEdgeOrientation(j,true); }
    }
    
    cards(elid) = card;
  }
  
  return(cards);
}

#endif
