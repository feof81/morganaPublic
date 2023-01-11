/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PVECTMANIP_HPP
#define PVECTMANIP_HPP


#include <set>
#include <map>
#include <assert.h>
#include <iostream>

#include "typesInterface.hpp"
#include "pMap.hpp"
#include "pVect.hpp"
#include "pMapItemSendRecv.h"
#include "traitsSegmentationUtility.hpp"


/*! Serial manipulation class for \c pVect */
template<typename DATA, typename MAPITEM> class pVectManip
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMap<pMapItemSendRecv> SENDRECV;
    typedef pVect<DATA,MAPITEM>    PVECT;
    typedef sVect<DATA>            DATAVECT;
    typedef pMap<MAPITEM>          PMAP;
    //@}
  
    /*! @name Internal data - finder */ //@{
  public:
    Teuchos::RCP<PVECT>  vect;
    bool           vectLoaded;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    pVectManip();
    pVectManip(const Teuchos::RCP<PVECT> & Vect);
    pVectManip(PVECT & Vect);
    void setVect(const Teuchos::RCP<PVECT> & Vect);
    void setVect(PVECT & Vect);
    //@}
    
    /*! @name Set-Get functions */ //@{
  public:
    /*! Set \c vect to normal indexing */
    void setNormalIndexing();
    void setNormalIndexing(PVECT & inVect) const;
    void setNormalIndexing(Teuchos::RCP<PVECT> & inVect) const;
    
    /*! Re-order the position of the items in the vector using the \c lid data */
    void setIndexing();
    void setIndexing(PVECT & inVect) const;
    void setIndexing(Teuchos::RCP<PVECT> & inVect) const;
    
    /*! Returns the maximum \c gid of \c vect */
    UInt getMaxGid();
    UInt getMaxGid(const PVECT & inVect) const;
    UInt getMaxGid(const Teuchos::RCP<const PVECT> & inVect) const;
    //@}
    
    /*! @name Group functions */ //@{
  public:
    /*! Cut \c vect into \c maxSize segments */
    void segmentationSimple(sVect<PVECT> & targetVects,
                            const UInt   & maxSize) const;
    
    void segmentationSimple(sVect<PVECT> & targetVects,
                            const UInt   & maxSize,
                            const PVECT  & inVect) const;
    
    void segmentationSimple(sVect<PVECT> & targetVects,
                            const UInt   & maxSize,
                            const Teuchos::RCP<const PVECT> & inVect) const;
    
    /*! Cut \c vect into \c segments segments. Items are assigned to a specific segment using the \c gid key.
    The maximum of all the \c gid s on all the processors is passed with \c maxGid */
    void segmentationNormal(sVect<PVECT> & targetVects,
                            const UInt   & segments,
                            const UInt   & maxGid) const;
    
    void segmentationNormal(sVect<PVECT> & targetVects,
                            const UInt   & segments,
                            const UInt   & maxGid,
                            const PVECT  & inVect) const;
    
    void segmentationNormal(sVect<PVECT> & targetVects,
                            const UInt   & segments,
                            const UInt   & maxGid,
                            const Teuchos::RCP<const PVECT> & inVect) const;
    
    /*! Cut \c vect into \c maxPid segments using the \c pid key */
    void segmentationPid(sVect<PVECT> & targetVects,
                         const UInt   & maxPid) const;
    
    void segmentationPid(sVect<PVECT> & targetVects,
                         const UInt   & maxPid,
                         const PVECT  & inVect) const;
    
    void segmentationPid(sVect<PVECT> & targetVects,
                         const UInt   & maxPid,
                         const Teuchos::RCP<const PVECT> & inVect) const;
    
    /*! Cut \c vect into \c maxPid segments using the sending-id (sid) in \c mapSend as a key */
    void segmentationMap(sVect<PVECT>   & targetVects,
                         const SENDRECV & mapSend,
                         const UInt     & maxPid) const;
    
    void segmentationMap(sVect<PVECT>   & targetVects,
                         const SENDRECV & mapSend,
                         const UInt     & maxPid,
                         const PVECT    & inVect) const;
    
    void segmentationMap(sVect<PVECT>   & targetVects,
                         const SENDRECV & mapSend,
                         const UInt     & maxPid,
                         const Teuchos::RCP<const PVECT> & inVect) const;
    
    /*! Cut \c vect into \c maxPid segments using the sending-id (sid) in \c mapSend as a key,
     *  recursive, high performance, version \c targetVects should be already allocated and the map must not change */
    void segmentationMapR(sVect<PVECT>   & targetVects,
                          const SENDRECV & mapSend,
                          const UInt     & maxPid) const;
    
    void segmentationMapR(sVect<PVECT>   & targetVects,
                          const SENDRECV & mapSend,
                          const UInt     & maxPid,
                          const PVECT    & inVect) const;
    
    void segmentationMapR(sVect<PVECT>   & targetVects,
                          const SENDRECV & mapSend,
                          const UInt     & maxPid,
                          const Teuchos::RCP<const PVECT> & inVect) const;
    
    /*! Cut \c vect into \c maxPid segments using the \c data themselves as a key.
    The maximum and minimum data \c maxData and \c minData should be provided  */
    void segmentationData(sVect<PVECT> & targetVects,
                          const UInt   & maxPid,
                          const DATA   & minData,
                          const DATA   & maxData) const;
    
    void segmentationData(sVect<PVECT> & targetVects,
                          const UInt   & maxPid,
                          const DATA   & minData,
                          const DATA   & maxData,
                          const PVECT  & inVect) const;
    
    void segmentationData(sVect<PVECT> & targetVects,
                          const UInt   & maxPid,
                          const DATA   & minData,
                          const DATA   & maxData,
                          const Teuchos::RCP<const PVECT> & inVect) const;
    
    /*! Re-order the data using a multi-set, \c key = data */
    void orderData();
    void orderData(PVECT & inVect);
    void orderData(Teuchos::RCP<PVECT> & inVect);
    
    /*! Re-order the data using a multi-set, \c key = gid */
    void orderGid();
    void orderGid(PVECT & inVect);
    void orderGid(Teuchos::RCP<PVECT> & inVect);
    
    /*! Simple merge into \c vect of the \c sourceMaps. \c vect is not cleared */
    void mergeSimple(sVect<PVECT> & sourceVects);
    
    void mergeSimple(sVect<PVECT> & sourceVects,
                           PVECT  & outVect);
    
    void mergeSimple(sVect<PVECT>        & sourceVects,
                     Teuchos::RCP<PVECT> & outVect);
    
    /*! Simple merge into \c vect of the \c sourceMaps. \c vect is not cleared
        Recursive version the map is not changed and the outVect must have already been properly sized  */
    void mergeSimpleR(sVect<PVECT> & sourceVects);
    
    void mergeSimpleR(sVect<PVECT> & sourceVects,
                            PVECT  & outVect);
    
    void mergeSimpleR(sVect<PVECT>        & sourceVects,
                      Teuchos::RCP<PVECT> & outVect);
    
     /*! Merge into \c vect of the \c sourceMaps. The repeated \c gid s are eliminated so that there is only one entry for each \c gid.
    \c vect is cleared*/
    void unionExclusive(sVect<PVECT> & sourceVects);
    
    void unionExclusive(sVect<PVECT> & sourceVects,
                               PVECT & outVect);
    
    void unionExclusive(sVect<PVECT>        & sourceVects,
                        Teuchos::RCP<PVECT> & outVect);
    
    /*! Given a vector eliminates the elements with the same \c gid */
    void unionExclusive();
    void unionExclusive(PVECT & outVect);
    void unionExclusive(Teuchos::RCP<PVECT> & outVect);
    //@}
    
    /*! @name Other functions */ //@{
  public:
    size_t memSize() const;
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
pVectManip<DATA,MAPITEM>::
pVectManip() : vectLoaded(false)
{
}

template<typename DATA, typename MAPITEM>
pVectManip<DATA,MAPITEM>::
pVectManip(const Teuchos::RCP<PVECT> & Vect) : vectLoaded(true), vect(Vect)
{
}

template<typename DATA, typename MAPITEM>
pVectManip<DATA,MAPITEM>::
pVectManip(PVECT & Vect) : vectLoaded(true)
{
  vect = Teuchos::rcpFromRef(Vect);
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
setVect(const Teuchos::RCP<PVECT> & Vect)
{
  vectLoaded = true;
  vect       = Vect;
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
setVect(PVECT & Vect)
{
  vectLoaded = true;
  vect       = Teuchos::rcpFromRef(Vect);
}



//_________________________________________________________________________________________________
// SET-GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
setNormalIndexing()
{
  assert(vectLoaded);
  
  pMapManip<MAPITEM> manipulator;
  manipulator.setNormalIndexing(vect->getMapRcp());
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
setNormalIndexing(PVECT & inVect) const
{
  pMapManip<MAPITEM> manipulator;
  manipulator.setNormalIndexing(inVect.getMapRcp());
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
setNormalIndexing(Teuchos::RCP<PVECT> & inVect) const
{
  setNormalIndexing(*inVect);
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
setIndexing()
{
  assert(vectLoaded);
  
  //Vector copy
  PVECT vectCopy(*vect);
  
  UInt lid;
  for(UInt i=1; i <= vect->size(); ++i)
  {
    lid = vect->getMapL(i).getLid();
    assert(lid <= vect->size());
    
    vectCopy.getMapL(lid) = vect->getMapL(i);
    vectCopy.getL(lid)    = vect->getL(i);
  }
  
  *vect = vectCopy;
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
setIndexing(PVECT & inVect) const
{
  //Vector copy
  PVECT vectCopy(inVect);
  
  UInt lid;
  for(UInt i=1; i <= inVect.size(); ++i)
  {
    lid = inVect.getMapL(i).getLid();
    assert(lid <= inVect.size());
    
    vectCopy.getMapL(lid) = inVect.getMapL(i);
    vectCopy.getL(lid)    = inVect.getL(i);
  }
  
  inVect = vectCopy;
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
setIndexing(Teuchos::RCP<PVECT> & inVect) const
{
  setIndexing(*inVect);
}

template<typename DATA, typename MAPITEM>
UInt
pVectManip<DATA,MAPITEM>::
getMaxGid()
{
  assert(vectLoaded);
  
  pMapManip<MAPITEM> manipulator;
  return(manipulator.getMaxGid(vect->getMapRcp()));
}

template<typename DATA, typename MAPITEM>
UInt
pVectManip<DATA,MAPITEM>::
getMaxGid(const PVECT & inVect) const
{
  pMapManip<MAPITEM> manipulator;
  return(manipulator.getMaxGid(inVect.getMapRcp()));
}

template<typename DATA, typename MAPITEM>
UInt
pVectManip<DATA,MAPITEM>::
getMaxGid(const Teuchos::RCP<const PVECT> & inVect) const
{
  return(getMaxGid(*inVect));
}



//_________________________________________________________________________________________________
// GROUP FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationSimple(sVect<PVECT> & targetVects,
                   const UInt   & maxSize) const
{
  assert(vectLoaded);
   
  UInt j=1, k=1;
  targetVects.clear();
  targetVects.resize(1);
  
  for(UInt i=1; i <= vect->size(); ++i)
  {
    targetVects(k).push_back(vect->get(i),vect->getMapL(i));
    
    if(j > maxSize)
    {
      j = 1;
      ++k;
      targetVects.push_back(PVECT());
    }
    else
    {
      ++j;
    }
  }
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationSimple(sVect<PVECT> & targetVects,
                   const UInt   & maxSize,
                   const PVECT  & inVect) const
{
  UInt j=1, k=1;
  targetVects.clear();
  targetVects.resize(1);
  
  for(UInt i=1; i <= inVect.size(); ++i)
  {
    targetVects(k).push_back(inVect.get(i),inVect.getMapL(i));
    
    if(j > maxSize)
    {
      j = 1;
      ++k;
      targetVects.push_back(PVECT());
    }
    else
    {
      ++j;
    }
  }
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationSimple(sVect<PVECT> & targetVects,
                   const UInt   & maxSize,
                   const Teuchos::RCP<const PVECT> & inVect) const
{
  segmentationSimple(targetVects,maxSize,*inVect);
}
    
template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationNormal(sVect<PVECT> & targetVects,
                   const UInt   & segments,
                   const UInt   & maxGid) const
{
  assert(vectLoaded);
  
  //Check normal indexing
  pMapManip<MAPITEM> cheker;
  assert(cheker.checkNormalIndexing(vect->getMapRcp()));
  
  //Map segmentation
  sVect<PMAP> targetMaps;
  
  pMapManip<MAPITEM> manipulator;
  manipulator.segmentationNormal(targetMaps,segments,maxGid,vect->getMapRcp());
  
  //Data segmentation
  UInt    lid;
  DATA    data;
  MAPITEM mapItem;
  
  targetVects.clear();
  targetVects.resize(targetMaps.size());
  
  for(UInt k=1; k <= targetMaps.size(); ++k)
  {
    for(UInt i=1; i <= targetMaps(k).size(); ++i)
    {
      mapItem = targetMaps(k).get(i);
      lid     = mapItem.getLid();
      
      assert(lid <= vect->size());
      data = vect->getL(lid);
      
      targetVects(k).push_back(mapItem,data);
    }
  }
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationNormal(sVect<PVECT> & targetVects,
                   const UInt   & segments,
                   const UInt   & maxGid,
                   const PVECT  & inVect) const
{
  //Check normal indexing
  pMapManip<MAPITEM> cheker;
  assert(cheker.checkNormalIndexing(inVect.getMapRcp()));
  
  //Map segmentation
  sVect<PMAP> targetMaps;
  
  pMapManip<MAPITEM> manipulator;
  manipulator.segmentationNormal(targetMaps,segments,maxGid,inVect.getMapRcp());
  
  //Data segmentation
  UInt    lid;
  DATA    data;
  MAPITEM mapItem;
  
  targetVects.clear();
  targetVects.resize(targetMaps.size());
  
  for(UInt k=1; k <= targetMaps.size(); ++k)
  {
    for(UInt i=1; i <= targetMaps(k).size(); ++i)
    {
      mapItem = targetMaps(k).get(i);
      lid     = mapItem.getLid();
      
      assert(lid <= inVect.size());
      data = inVect.getL(lid);
      
      targetVects(k).push_back(mapItem,data);
    }
  }
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationNormal(sVect<PVECT> & targetVects,
                   const UInt   & segments,
                   const UInt   & maxGid,
                   const Teuchos::RCP<const PVECT> & inVect) const
{
  segmentationNormal(targetVects,segments,maxGid,*inVect);
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationPid(sVect<PVECT> & targetVects,
                const UInt   & maxPid) const
{
  assert(vectLoaded);
  
  //Check normal indexing
  pMapManip<MAPITEM> cheker;
  assert(cheker.checkNormalIndexing(vect->getMapRcp()));
  
  //Map segmentation
  sVect<PMAP> targetMaps;
  
  pMapManip<MAPITEM> manipulator;
  manipulator.segmentationPid(targetMaps,maxPid,vect->getMapRcp());
  
  //Data segmentation
  UInt    lid;
  DATA    data;
  MAPITEM mapItem;
  
  targetVects.clear();
  targetVects.resize(targetMaps.size());

  for(UInt k=1; k <= targetMaps.size(); ++k)
  {
    for(UInt i=1; i <= targetMaps(k).size(); ++i)
    {
      mapItem = targetMaps(k).get(i);
      lid     = mapItem.getLid();
      data    = vect->getL(lid);
      
      targetVects(k).push_back(mapItem,data);
    }
  }
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationPid(sVect<PVECT> & targetVects,
                const UInt   & maxPid,
                const PVECT  & inVect) const
{
  //Check normal indexing
  pMapManip<MAPITEM> cheker;
  assert(cheker.checkNormalIndexing(inVect.getMapRcp()));
  
  //Map segmentation
  sVect<PMAP> targetMaps;
  
  pMapManip<MAPITEM> manipulator;
  manipulator.segmentationPid(targetMaps,maxPid,inVect.getMapRcp());
  
  //Data segmentation
  UInt    lid;
  DATA    data;
  MAPITEM mapItem;
  
  targetVects.clear();
  targetVects.resize(targetMaps.size());

  for(UInt k=1; k <= targetMaps.size(); ++k)
  {
    for(UInt i=1; i <= targetMaps(k).size(); ++i)
    {
      mapItem = targetMaps(k).get(i);
      lid     = mapItem.getLid();
      data    = inVect.getL(lid);
      
      targetVects(k).push_back(mapItem,data);
    }
  }
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationPid(sVect<PVECT> & targetVects,
                const UInt   & maxPid,
                const Teuchos::RCP<const PVECT> & inVect) const
{
  segmentationPid(targetVects,maxPid,*inVect);
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationMap(      sVect<PVECT> & targetVects,
                const SENDRECV     & mapSend,
                const UInt         & maxPid) const
{
  assert(vectLoaded);
  
  //Clearing
  targetVects.clear();
  targetVects.resize(maxPid);
  
  //Segmentation
  UInt rid, gid;
  
  for(UInt i=1; i <= mapSend.size(); ++i)
  {
    rid = mapSend(i).getRid();
    gid = mapSend(i).getGid();
    
    assert(rid < maxPid);
    assert(gid == vect->getMapG(gid).getGid());
    targetVects[rid].push_back(vect->getMapG(gid), vect->getG(gid));
  }
  
  //Updating finders
  for(UInt i=1; i <= targetVects.size(); ++i)
  {
    targetVects(i).updateFinder();
  }
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationMap(      sVect<PVECT> & targetVects,
                const SENDRECV     & mapSend,
                const UInt         & maxPid,
                const PVECT        & inVect) const
{
  //Clearing
  targetVects.clear();
  targetVects.resize(maxPid);
  
  //Segmentation
  UInt rid, gid;
  
  for(UInt i=1; i <= mapSend.size(); ++i)
  {
    rid = mapSend(i).getRid();
    gid = mapSend(i).getGid();
    
    assert(rid < maxPid);
    assert(gid == inVect.getMapG(gid).getGid());
    targetVects[rid].push_back(inVect.getMapG(gid), inVect.getG(gid));
  }
  
  //Updating finders
  for(UInt i=1; i <= targetVects.size(); ++i)
  {
    targetVects(i).updateFinder();
  }
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationMap(      sVect<PVECT> & targetVects,
                const SENDRECV     & mapSend,
                const UInt         & maxPid,
                const Teuchos::RCP<const PVECT> & inVect) const
{
  segmentationMap(targetVects,mapSend,maxPid,*inVect);
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationMapR(      sVect<PVECT> & targetVects,
                 const SENDRECV     & mapSend,
                 const UInt         & maxPid) const
{
  //Assert
  assert(vectLoaded);
  assert(targetVects.size() == maxPid);
  
  //Segmentation
  UInt rid, gid;
  sVect<UInt> indexes(targetVects.size());
  
  for(UInt i=1; i <= indexes.size(); ++i)
  { indexes(i) = 1;}
  
  for(UInt i=1; i <= mapSend.size(); ++i)
  {    
    rid = mapSend(i).getRid();
    gid = mapSend(i).getGid();
    
    assert(rid < maxPid);
    assert(gid == vect->getMapG(gid).getGid());
    
    targetVects[rid](indexes[rid]) = vect->getG(gid);
    indexes[rid]++;
  }
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationMapR(      sVect<PVECT> & targetVects,
                 const SENDRECV     & mapSend,
                 const UInt         & maxPid,
                 const PVECT        & inVect) const
{
  //Assert
  assert(targetVects.size() == maxPid);
  
  //Segmentation
  UInt rid, gid;
  sVect<UInt> indexes(targetVects.size());
  
  for(UInt i=1; i <= indexes.size(); ++i)
  { indexes(i) = 1;}
  
  for(UInt i=1; i <= mapSend.size(); ++i)
  {
    rid = mapSend(i).getRid();
    gid = mapSend(i).getGid();
    
    assert(rid < maxPid);
    assert(gid == inVect.getMapG(gid).getGid());
    
    targetVects[rid](indexes[rid]) = inVect.getG(gid);
    indexes[rid]++;
  }
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationMapR(      sVect<PVECT> & targetVects,
                 const SENDRECV     & mapSend,
                 const UInt         & maxPid,
                 const Teuchos::RCP<const PVECT> & inVect) const
{
  segmentationMapR(targetVects,mapSend,maxPid,*inVect);
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationData(sVect<PVECT> & targetVects,
                 const UInt   & maxPid,
                 const DATA   & minData,
                 const DATA   & maxData) const
{
  assert(vectLoaded);
  
  //Allocation 
  UInt k;
  
  //Vect array
  targetVects.clear();
  targetVects.resize(maxPid);
  
  //Segmenting utility
  dataSegmentationUtility<DATA> utility(maxPid,minData,maxData);
  
  //Segmenting
  for(UInt i=1; i <= vect->size(); ++i)
  {
    k = utility.index(vect->getL(i));
    targetVects(k).push_back(vect->getMapL(i),vect->getL(i));
  }
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationData(sVect<PVECT> & targetVects,
                 const UInt   & maxPid,
                 const DATA   & minData,
                 const DATA   & maxData,
                 const PVECT  & inVect) const
{
  //Allocation 
  UInt k;
  
  //Vect array
  targetVects.clear();
  targetVects.resize(maxPid);
  
  //Segmenting utility
  dataSegmentationUtility<DATA> utility(maxPid,minData,maxData);
  
  //Segmenting
  for(UInt i=1; i <= inVect.size(); ++i)
  {
    k = utility.index(inVect.getL(i));
    targetVects(k).push_back(inVect.getMapL(i),inVect.getL(i));
  }
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
segmentationData(sVect<PVECT> & targetVects,
                 const UInt   & maxPid,
                 const DATA   & minData,
                 const DATA   & maxData,
                 const Teuchos::RCP<const PVECT> & inVect) const
{
  segmentationData(targetVects,maxPid,minData,maxData,*inVect);
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
orderData()
{
  assert(vectLoaded);
  
  typedef pair<DATA,MAPITEM>                        PAIR;
  typedef typename multimap<DATA,MAPITEM>::iterator ITERATOR;
  
  PAIR     pair;
  ITERATOR iter;
  PVECT    vectOrder;
  
  multimap<DATA,MAPITEM> trunk;
  
  //Uploading 
  for(UInt i=1; i <= vect->size(); ++i)
  {
    pair.first  = vect->getDataL(i);
    pair.second = vect->getMapL(i);
    
    trunk.insert(pair);
  }
  
  //Downloading
  vectOrder.reserve(vect->size());
  
  for(iter = trunk.begin(); iter != trunk.end(); ++iter)
  { vectOrder.push_back((*iter).first, (*iter).second); }
  
  vectOrder.updateFinder();
  
  //Copying
  *vect = vectOrder;
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
orderData(PVECT & inVect)
{
  typedef pair<DATA,MAPITEM>                        PAIR;
  typedef typename multimap<DATA,MAPITEM>::iterator ITERATOR;
  
  PAIR     pair;
  ITERATOR iter;
  PVECT    vectOrder;
  
  multimap<DATA,MAPITEM> trunk;
  
  //Uploading 
  for(UInt i=1; i <= inVect.size(); ++i)
  {
    pair.first  = inVect.getDataL(i);
    pair.second = inVect.getMapL(i);
    
    trunk.insert(pair);
  }
  
  //Downloading
  vectOrder.reserve(inVect.size());
  
  for(iter = trunk.begin(); iter != trunk.end(); ++iter)
  { vectOrder.push_back((*iter).first, (*iter).second); }
  
  vectOrder.updateFinder();
  
  //Copying
  inVect = vectOrder;
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
orderData(Teuchos::RCP<PVECT> & inVect)
{
  orderData(*inVect);
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
orderGid()
{
  assert(vectLoaded);
  
  typedef pair<MAPITEM,DATA>                        PAIR;
  typedef typename multimap<MAPITEM,DATA>::iterator ITERATOR;
  
  PAIR     pair;
  ITERATOR iter;
  PVECT    vectOrder;
  
  multimap<MAPITEM,DATA> trunk;
  
  //Uploading 
  for(UInt i=1; i <= vect->size(); ++i)
  {
    pair.first  = vect->getMapL(i);
    pair.second = vect->getDataL(i);
    
    trunk.insert(pair);
  }
  
  //Downloading
  vectOrder.reserve(vect->size());
  
  for(iter = trunk.begin(); iter != trunk.end(); ++iter)
  { vectOrder.push_back((*iter).second, (*iter).first); }
  
  vectOrder.updateFinder();
  
  //Copying
  *vect = vectOrder;
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
orderGid(PVECT & inVect)
{
  typedef pair<DATA,MAPITEM>                        PAIR;
  typedef typename multimap<DATA,MAPITEM>::iterator ITERATOR;
  
  PAIR     pair;
  ITERATOR iter;
  PVECT    vectOrder;
  
  multimap<DATA,MAPITEM> trunk;
  
  //Uploading 
  for(UInt i=1; i <= inVect.size(); ++i)
  {
    pair.first  = inVect.getDataL(i);
    pair.second = inVect.getMapL(i);
    
    trunk.insert(pair);
  }
  
  //Downloading
  vectOrder.reserve(inVect.size());
  
  for(iter = trunk.begin(); iter != trunk.end(); ++iter)
  { vectOrder.push_back((*iter).first, (*iter).second); }
  
  vectOrder.updateFinder();
  
  //Copying
  inVect = vectOrder;
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
orderGid(Teuchos::RCP<PVECT> & inVect)
{
  orderGid(*inVect);
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
mergeSimple(sVect<PVECT> & sourceVects)
{
  assert(vectLoaded);
  assert(vect->getMapRef().size() == vect->getDataRef().size());
  
  pMap<MAPITEM>  mapTot;
  sVect<DATA>   dataTot;
  
  for(UInt k=1; k <= sourceVects.size(); ++k)
  {
    for(UInt i=1; i <= sourceVects(k).size(); ++i)
    {
      mapTot.push_back(sourceVects(k).getMapL(i));
      dataTot.push_back(sourceVects(k).getL(i));
    }
  }
  
  vect->addData(mapTot,dataTot);
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
mergeSimple(sVect<PVECT> & sourceVects,
            PVECT        & outVect)
{
  assert(outVect.getMapRef().size() == outVect.getDataRef().size());
  
  pMap<MAPITEM>  mapTot;
  sVect<DATA>   dataTot;
  
  for(UInt k=1; k <= sourceVects.size(); ++k)
  {
    for(UInt i=1; i <= sourceVects(k).size(); ++i)
    {
      mapTot.push_back(sourceVects(k).getMapL(i));
      dataTot.push_back(sourceVects(k).getL(i));
    }
  }
  
  outVect.addData(mapTot,dataTot);
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
mergeSimple(sVect<PVECT>        & sourceVects,
            Teuchos::RCP<PVECT> & outVect)
{
  mergeSimple(sourceVects,*outVect);
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
mergeSimpleR(sVect<PVECT> & sourceVects)
{
  assert(vectLoaded);
  assert(vect->getMapRef().size() == vect->getDataRef().size());
  
  UInt index = 1;
  
  for(UInt k=1; k <= sourceVects.size(); ++k)
  {
    for(UInt i=1; i <= sourceVects(k).size(); ++i)
    {
      vect->operator()(index) = sourceVects(k).getL(i);
      index++;
    }
  }
  
  assert(vect->size() == (index-1));
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
mergeSimpleR(sVect<PVECT> & sourceVects,
             PVECT        & outVect)
{
  assert(outVect.getMapRef().size() == outVect.getDataRef().size());
  
  UInt index = 1;
  
  for(UInt k=1; k <= sourceVects.size(); ++k)
  {
    for(UInt i=1; i <= sourceVects(k).size(); ++i)
    {
      outVect(index) = sourceVects(k).getL(i);
      index++;
    }
  }
  
  assert(outVect.size() == (index-1));
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
mergeSimpleR(sVect<PVECT>        & sourceVects,
             Teuchos::RCP<PVECT> & outVect)
{
  mergeSimpleR(sourceVects,*outVect);
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
unionExclusive(sVect<PVECT> & sourceVects)
{
  assert(vectLoaded);
  assert(vect->getMapRef().size() == vect->getDataRef().size());
  
  //Adding the actual vector to the source list
  vect->clear();
  sourceVects.push_back(*vect);
  
  //Setting the pid to the corresponging sourceVect
  sVect<PMAP> sourceMaps(sourceVects.size());
  
  for(UInt k=1; k <= sourceVects.size(); ++k)
  {
    sourceMaps(k) = sourceVects(k).getMapRef();
    
    for(UInt i=1; i <= sourceVects(k).size(); ++i)
    {
      sourceMaps(k).get(i).setPid(k);
    }
  }

  //Creating a common map 
  pMap<MAPITEM>      newMap;
  pMapManip<MAPITEM> mapManipulator;
  mapManipulator.unionExclusive(sourceMaps,newMap);
  
  //Updating the pVect
  UInt pid, lid;
  DATA data;
  
  for(UInt i=1; i <= newMap.size(); ++i)
  {
    lid = newMap(i).getLid();
    pid = newMap(i).getPid();
    
    data = sourceVects(pid).get(lid);
    vect->push_back(data,newMap(i));
  }
  
  vect->updateFinder();
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
unionExclusive(sVect<PVECT> & sourceVects,
               PVECT        & outVect)
{
  assert(outVect.getMapRef().size() == outVect.getDataRef().size());
  
  //Adding the actual vector to the source list
  outVect.clear();
  sourceVects.push_back(outVect);
  
  //Setting the pid to the corresponging sourceVect
  sVect<PMAP> sourceMaps(sourceVects.size());
  
  for(UInt k=1; k <= sourceVects.size(); ++k)
  {
    sourceMaps(k) = sourceVects(k).getMapRef();
    
    for(UInt i=1; i <= sourceVects(k).size(); ++i)
    {
      sourceMaps(k).get(i).setPid(k);
    }
  }

  //Creating a common map 
  pMap<MAPITEM>      newMap;
  pMapManip<MAPITEM> mapManipulator;
  mapManipulator.unionExclusive(sourceMaps,newMap);
  
  //Updating the pVect
  UInt pid, lid;
  DATA data;
  
  for(UInt i=1; i <= newMap.size(); ++i)
  {
    lid = newMap(i).getLid();
    pid = newMap(i).getPid();
    
    data = sourceVects(pid).get(lid);
    outVect.push_back(data,newMap(i));
  }
  
  outVect.updateFinder();
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
unionExclusive(sVect<PVECT>        & sourceVects,
               Teuchos::RCP<PVECT> & outVect)
{
  unionExclusive(sourceVects,*outVect);
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
unionExclusive()
{
  typedef typename map<MAPITEM,DATA>::iterator ITER;
  
  //The stl map
  map<MAPITEM,DATA>  serialMap;
  pair<MAPITEM,DATA> pair;
  
  for(UInt i=1; i <= vect->size(); ++i)
  {
    pair.first  = vect->getMapL(i);
    pair.second = vect->getDataL(i);
    
    serialMap.insert(pair);
  }
  
  //De-serialization
  vect->clear();
  
  for(ITER it = serialMap.begin(); it != serialMap.end(); ++it)
  {
    vect->push_back(it->first,it->second);
  }
  
  //Updating
  vect->updateFinder();
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
unionExclusive(PVECT & outVect)
{
  typedef typename map<MAPITEM,DATA>::iterator ITER;
  
  //The stl map
  map<MAPITEM,DATA>  serialMap;
  pair<MAPITEM,DATA> pair;
  
  for(UInt i=1; i <= outVect.size(); ++i)
  {
    pair.first  = outVect.getMapL(i);
    pair.second = outVect.getDataL(i);
    
    serialMap.insert(pair);
  }
  
  //De-serialization
  outVect.clear();
  
  for(ITER it = serialMap.begin(); it != serialMap.end(); ++it)
  {
    outVect.push_back(it->first,it->second);
  }
  
  //Updating
  outVect.updateFinder();
}

template<typename DATA, typename MAPITEM>
void
pVectManip<DATA,MAPITEM>::
unionExclusive(Teuchos::RCP<PVECT> & outVect)
{
  unionExclusive(*outVect);
}


//_________________________________________________________________________________________________
// OTHER FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
size_t
pVectManip<DATA,MAPITEM>::
memSize() const
{
  return(sizeof(bool));
}

#endif
