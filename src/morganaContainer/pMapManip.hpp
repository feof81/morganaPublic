/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PMAPMANIP_HPP
#define PMAPMANIP_HPP

#include <set>
#include <assert.h>
#include <iostream>

#include "Teuchos_RCP.hpp"

#include "pMap.hpp"
#include "pMapItem.h"

using namespace std;


/*! Serial manipulation class for \c pMap */
template<typename ITEM> class pMapManip
{    
    /*! @name Typedefs */ //@{
  public:
    typedef pMap<ITEM> PMAP;
    //@}
  
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<PMAP> map;
     //@}
    
    /*! @name Internal data - finder */ //@{
  public:
    bool mapLoaded;
    bool finderOk;
    set<ITEM> container;
    //@}
    
    
    /*! @name Constructors and set functions */ //@{
  public:
    pMapManip();
    pMapManip(const Teuchos::RCP<PMAP> & Map);
    pMapManip(PMAP & Map);
    pMapManip(const pMapManip & Pmanip);
    pMapManip operator=(const pMapManip & Pmanip);
    void setMap(const Teuchos::RCP<PMAP> & Map);
    void setMap(PMAP & Map);
    //@}
    
    
    /*! @name Search functions */ //@{
  public:
    void buildFinder();
    void resetFinder();
    bool isItem(const ITEM & Item) const;
    UInt getLidItem(const ITEM & Item) const;
    //@}
    
    
    /*! @name Set-Get and indexing functions */ //@{
  public:
    /*! Set \c map to normal indexing */
    void setNormalIndexing();
    void setNormalIndexing(PMAP & inMap) const;
    void setNormalIndexing(Teuchos::RCP<PMAP> & inMap) const;
    
    /*! Check the normal indexing */
    bool checkNormalIndexing();
    bool checkNormalIndexing(const PMAP & inMap) const;
    bool checkNormalIndexing(const Teuchos::RCP<const PMAP> & inMap) const;
    
    /*! Re-order the position of the items in the vector using the \c lid data */
    void setIndexing();
    void setIndexing(PMAP & inMap) const;
    void setIndexing(Teuchos::RCP<PMAP> & inMap) const;
    
    /*! Returns the maximum \c gid of \c vect */
    UInt getMaxGid() const;
    UInt getMaxGid(const PMAP & inMap) const;
    UInt getMaxGid(const Teuchos::RCP<const PMAP> & inMap) const;
    
    /*! The graphs are equal? */
    bool isEqual(const PMAP & mapB);
    bool isEqual(const PMAP & mapA, const PMAP & mapB) const;
    bool isEqual(const Teuchos::RCP<const PMAP> & mapA, Teuchos::RCP<const PMAP> & mapB) const;
    //@}
    
    
    /*! @name Group functions */ //@{
  public:
    /*! Cut \c map into a number of segments such that none of them exceeds \c maxSize */
    void segmentationSimple(sVect<PMAP> & targetMaps, const UInt & maxSize) const;
    void segmentationSimple(sVect<PMAP> & targetMaps, const UInt & maxSize, const PMAP & inMap) const;
    void segmentationSimple(sVect<PMAP> & targetMaps, const UInt & maxSize, const Teuchos::RCP<const PMAP> & inMap) const;
    
    /*! Cut \c map into \c segments segments. Items are assigned to a specific segment using the \c gid key.
    The maximum of all the \c gid s on all the processors is passed with \c maxGid */
    void segmentationNormal(sVect<PMAP> & targetMaps, const UInt & segments, const UInt & maxGid) const;
    void segmentationNormal(sVect<PMAP> & targetMaps, const UInt & segments, const UInt & maxGid, const PMAP & inMap) const;
    void segmentationNormal(sVect<PMAP> & targetMaps, const UInt & segments, const UInt & maxGid, const Teuchos::RCP<const PMAP> & inMap) const;
    
    /*! Cut \c map into \c maxPid segments using the \c pid key */
    void segmentationPid(sVect<PMAP> & targetMaps, const UInt & maxPid) const;
    void segmentationPid(sVect<PMAP> & targetMaps, const UInt & maxPid, const PMAP & inMap) const;
    void segmentationPid(sVect<PMAP> & targetMaps, const UInt & maxPid, const Teuchos::RCP<const PMAP> & inMap) const;
    
    /*! Simple merge into \c map of the \c sourceMaps. \c map is not cleared */
    void mergeSimple(sVect<PMAP> & sourceMaps);
    void mergeSimple(sVect<PMAP> & sourceMaps, PMAP & inMap);
    void mergeSimple(sVect<PMAP> & sourceMaps, Teuchos::RCP<PMAP> & inMap);
    
    /*! Re-order the data using a multi-set using gid*/
    void orderData();
    void orderData(PMAP & inMap);
    void orderData(Teuchos::RCP<PMAP> & inMap);
    
    /*! Merge into \c map of the \c sourceMaps. The repeated \c gid s are eliminated so that there is only one entry for each \c gid.
    \c map is not cleared*/
    void unionExclusive(sVect<PMAP> & sourceMaps);
    void unionExclusive(sVect<PMAP> & sourceMaps, PMAP & inMap);
    void unionExclusive(sVect<PMAP> & sourceMaps, Teuchos::RCP<PMAP> & inMap);
    //@}
    
    /*! @name Memory functions */ //@{
  public:
    size_t memSize() const;
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename ITEM>
pMapManip<ITEM>::
pMapManip() : mapLoaded(false), finderOk(false)
{
}

template<typename ITEM>
pMapManip<ITEM>::
pMapManip(const Teuchos::RCP<PMAP> & Map) : map(Map)
{
  mapLoaded = true;
  finderOk  = false;
}

template<typename ITEM>
pMapManip<ITEM>::
pMapManip(PMAP & Map) : mapLoaded(true), finderOk(false)
{ 
  map = Teuchos::rcpFromRef(Map);
}

template<typename ITEM>
pMapManip<ITEM>::
pMapManip(const pMapManip & Pmanip)
{
  mapLoaded = Pmanip.mapLoaded;
  finderOk  = Pmanip.finderOk;
  container = Pmanip.container;
}

template<typename ITEM>
pMapManip<ITEM>
pMapManip<ITEM>::
operator=(const pMapManip & Pmanip)
{
  mapLoaded = Pmanip.mapLoaded;
  finderOk  = Pmanip.finderOk;
  container = Pmanip.container;
  
  return(*this);
}

template<typename ITEM>
void
pMapManip<ITEM>::
setMap(const Teuchos::RCP<PMAP> & Map)
{
  mapLoaded = true;
  map = Map;
}

template<typename ITEM>
void
pMapManip<ITEM>::
setMap(PMAP & Map)
{
  mapLoaded = true;
  
  map = Teuchos::rcpFromRef(Map);
}



//_________________________________________________________________________________________________
// FINDING FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename ITEM>
void
pMapManip<ITEM>::
buildFinder()
{
  assert(mapLoaded);
  finderOk = true;
  
  container.clear();
  
  for(UInt i=1; i<=map->size(); ++i)
  {
    container.insert(map->get(i));
  }
}

template<typename ITEM>
void
pMapManip<ITEM>::
resetFinder()
{
  finderOk = false;
  container.clear();
}

template<typename ITEM>
bool
pMapManip<ITEM>::
isItem(const ITEM & Item) const
{
  assert(mapLoaded);
  assert(finderOk);
  
  return((bool)container.count(Item));
}

template<typename ITEM>
UInt
pMapManip<ITEM>::
getLidItem(const ITEM & Item) const
{
  assert(mapLoaded);
  assert(finderOk);
  assert((bool)container.count(Item));
  
  typename set<ITEM>::iterator iteratore;
  pMapItem itemDelSet;
  
  iteratore  = container.find(Item);
  itemDelSet = *iteratore;
  
  return(itemDelSet.getLid());
}



//_________________________________________________________________________________________________
// SET-GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename ITEM>
void
pMapManip<ITEM>::
setNormalIndexing()
{
  assert(mapLoaded);
  
  for(UInt i=1; i <= map->size(); ++i)
  { map->get(i).setLid(i); }
}

template<typename ITEM>
void
pMapManip<ITEM>::
setNormalIndexing(PMAP & inMap) const
{
  for(UInt i=1; i <= inMap.size(); ++i)
  { inMap.get(i).setLid(i); }
}
 
template<typename ITEM>
void
pMapManip<ITEM>::
setNormalIndexing(Teuchos::RCP<PMAP> & inMap) const
{
  setNormalIndexing(*inMap);
}

template<typename ITEM>
bool
pMapManip<ITEM>::
checkNormalIndexing()
{
  assert(mapLoaded);
  
  bool flag = true;
  
  for(UInt i=1; i <= map->size(); ++i)
  { flag = flag & (map->get(i).getLid() == i); }
  
  return(flag);
}

template<typename ITEM>
bool
pMapManip<ITEM>::
checkNormalIndexing(const PMAP & inMap) const
{
  bool flag = true;
  
  for(UInt i=1; i <= inMap.size(); ++i)
  { flag = flag & (inMap.get(i).getLid() == i); }
  
  return(flag);
}

template<typename ITEM>
bool
pMapManip<ITEM>::
checkNormalIndexing(const Teuchos::RCP<const PMAP> & inMap) const
{
  return(checkNormalIndexing(*inMap));
}

template<typename ITEM>
void
pMapManip<ITEM>::
setIndexing()
{
  assert(mapLoaded);
  
  //Map copy
  PMAP mapCopy(*map);
  
  UInt lid;
  for(UInt i=1; i <= map->size(); ++i)
  {
    lid = map->get(i).getLid();
    assert(lid <= map->size());
    
    mapCopy.get(lid) = map->get(i);
  }
  
  *map = mapCopy;
}

template<typename ITEM>
void
pMapManip<ITEM>::
setIndexing(PMAP & inMap) const
{
  //Map copy
  PMAP mapCopy(inMap);
  
  UInt lid;
  for(UInt i=1; i <= inMap.size(); ++i)
  {
    lid = inMap.get(i).getLid();
    assert(lid <= inMap.size());
    
    mapCopy.get(lid) = inMap.get(i);
  }
  
  inMap = mapCopy;
}

template<typename ITEM>
void
pMapManip<ITEM>::
setIndexing(Teuchos::RCP<PMAP> & inMap) const
{
  setIndexing(*inMap);
}

template<typename ITEM>
UInt
pMapManip<ITEM>::
getMaxGid() const
{
  assert(mapLoaded);
  
  UInt maxSize = 0;
  
  for(UInt i=1; i <= map->size(); ++i)
  { maxSize = max(maxSize, map->get(i).getGid()); }
  
  return(maxSize);
}

template<typename ITEM>
UInt
pMapManip<ITEM>::
getMaxGid(const PMAP & inMap) const
{
  UInt maxSize = 0;
  
  for(UInt i=1; i <= inMap.size(); ++i)
  { maxSize = max(maxSize, inMap.get(i).getGid()); }
  
  return(maxSize);
}

template<typename ITEM>
UInt
pMapManip<ITEM>::
getMaxGid(const Teuchos::RCP<const PMAP> & inMap) const
{
  return(getMaxGid(*inMap));
}

template<typename ITEM>
bool
pMapManip<ITEM>::
isEqual(const PMAP & mapB)
{
  //Length check 
  if(map->size() != mapB.size())
  {
    return(false);
  }
  
  //Fine check
  bool flag = true;
  
  for(UInt i=1; i <= map->size(); ++i)
  {
    flag = flag & !(map->get(i) != mapB.get(i));
  }
  
  return(flag);
}

template<typename ITEM>
bool
pMapManip<ITEM>::
isEqual(const PMAP & mapA, const PMAP & mapB) const
{
  //Length check 
  if(mapA.size() != mapB.size())
  {
    return(false);
  }
  
  //Fine check
  bool flag = true;
  
  for(UInt i=1; i <= mapA.size(); ++i)
  {
    flag = flag & !(mapA.get(i) != mapB.get(i));
  }
  
  return(flag);
}

template<typename ITEM>
bool
pMapManip<ITEM>::
isEqual(const Teuchos::RCP<const PMAP> & mapA, Teuchos::RCP<const PMAP> & mapB) const
{
  return(isEqual(*mapA,*mapB));
}



//_________________________________________________________________________________________________
// GROUP FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename ITEM>
void
pMapManip<ITEM>::
segmentationSimple(sVect<PMAP> & targetMaps, const UInt & maxSize) const
{
  assert(mapLoaded);
  
  UInt j=1, k=1;
  targetMaps.clear();
  targetMaps.resize(1);
  
  //The vector is correcty sized
  if(map->size() < maxSize)
  {
    targetMaps(1) = *map;
  }
  //The vector is too long and should be divided
  else
  {
    for(UInt i=1; i <= map->size(); ++i)
    {
      targetMaps(k).push_back(map->get(i));
      
      if(j > maxSize)
      {
	j = 1;
        ++k;
        targetMaps.push_back(PMAP());
	targetMaps(k).reserve(maxSize);
      }
      else
      {
        ++j;
      }
    }
  }
}

template<typename ITEM>
void
pMapManip<ITEM>::
segmentationSimple(sVect<PMAP> & targetMaps, const UInt & maxSize, const PMAP & inMap) const
{
  UInt j=1, k=1;
  targetMaps.clear();
  targetMaps.resize(1);
  
  //The vector is correcty sized
  if(inMap.size() < maxSize)
  {
    targetMaps(1) = inMap;
  }
  //The vector is too long and should be divided
  else
  {
    for(UInt i=1; i <= inMap.size(); ++i)
    {
      targetMaps(k).push_back(inMap.get(i));
      
      if(j > maxSize)
      {
	j = 1;
        ++k;
        targetMaps.push_back(PMAP());
	targetMaps(k).reserve(maxSize);
      }
      else
      {
        ++j;
      }
    }
  }
}

template<typename ITEM>
void
pMapManip<ITEM>::
segmentationSimple(sVect<PMAP> & targetMaps, const UInt & maxSize, const Teuchos::RCP<const PMAP> & inMap) const
{
  segmentationSimple(targetMaps,maxSize,*inMap);
}

template<typename ITEM>
void
pMapManip<ITEM>::
segmentationNormal(sVect<PMAP> & targetMaps, const UInt & segments, const UInt & maxGid) const
{
  assert(mapLoaded);
  
  UInt k, gid;
  UInt size = ceil(Real(maxGid) / Real(segments));
  
  //Resizing and clearing
  targetMaps.clear();
  targetMaps.resize(segments);
  
  //Segmenting
  for(UInt i=1; i<=map->size(); ++i)
  {
    gid = map->get(i).getGid();
    k   = ceil(Real(gid) / Real(size));
    
    targetMaps(k).push_back(map->get(i));
  }
}

template<typename ITEM>
void
pMapManip<ITEM>::
segmentationNormal(sVect<PMAP> & targetMaps, const UInt & segments, const UInt & maxGid, const PMAP & inMap) const
{
  UInt k, gid;
  UInt size = ceil(Real(maxGid) / Real(segments));
  
  //Resizing and clearing
  targetMaps.clear();
  targetMaps.resize(segments);
  
  //Segmenting
  for(UInt i=1; i<= inMap.size(); ++i)
  {
    gid = inMap.get(i).getGid();
    k   = ceil(Real(gid) / Real(size));
    
    targetMaps(k).push_back(inMap.get(i));
  }
}

template<typename ITEM>
void
pMapManip<ITEM>::
segmentationNormal(sVect<PMAP> & targetMaps, const UInt & segments, const UInt & maxGid, const Teuchos::RCP<const PMAP> & inMap) const
{
  segmentationNormal(targetMaps,segments,maxGid,*inMap);
}

template<typename ITEM>
void
pMapManip<ITEM>::
segmentationPid(sVect<PMAP> & targetMaps, const UInt & maxPid) const
{
  assert(mapLoaded);
  
  //Alloc
  UInt pid;
  
  //Resizing and clearing
  targetMaps.clear();
  targetMaps.resize(maxPid);
  
  for(UInt i=1; i <= map->size(); ++i)
  {
    pid = map->get(i).getPid();
    assert(pid < maxPid);
    
    targetMaps[pid].push_back(map->get(i));
  }
}

template<typename ITEM>
void
pMapManip<ITEM>::
segmentationPid(sVect<PMAP> & targetMaps, const UInt & maxPid, const PMAP & inMap) const
{
  //Alloc
  UInt pid;
  
  //Resizing and clearing
  targetMaps.clear();
  targetMaps.resize(maxPid);
  
  for(UInt i=1; i <= inMap.size(); ++i)
  {
    pid = inMap.get(i).getPid();
    assert(pid < maxPid);
    
    targetMaps[pid].push_back(inMap.get(i));
  }
}

template<typename ITEM>
void
pMapManip<ITEM>::
segmentationPid(sVect<PMAP> & targetMaps, const UInt & maxPid, const Teuchos::RCP<const PMAP> & inMap) const
{
  segmentationPid(targetMaps,maxPid,*inMap);
}

template<typename ITEM>
void
pMapManip<ITEM>::
mergeSimple(sVect<PMAP> & sourceMaps)
{
  assert(mapLoaded);
  
  for(UInt k=1; k <= sourceMaps.size(); ++k)
  {
    for(UInt j=1; j <= sourceMaps(k).size(); ++j)
    {
      map->push_back(sourceMaps(k).get(j));
    }
  }
}

template<typename ITEM>
void
pMapManip<ITEM>::
mergeSimple(sVect<PMAP> & sourceMaps, PMAP & inMap)
{
  for(UInt k=1; k <= sourceMaps.size(); ++k)
  {
    for(UInt j=1; j <= sourceMaps(k).size(); ++j)
    {
      inMap.push_back(sourceMaps(k).get(j));
    }
  }
}

template<typename ITEM>
void
pMapManip<ITEM>::
mergeSimple(sVect<PMAP> & sourceMaps, Teuchos::RCP<PMAP> & inMap)
{
  mergeSimple(sourceMaps,*inMap);
}

template<typename ITEM>
void
pMapManip<ITEM>::
orderData()
{
  //Data serialization
  typedef typename multiset<ITEM>::iterator ITER;
  multiset<ITEM> tempMap;
  
  for(UInt i=1; i <= map->size(); ++i)
  {
    tempMap.insert(map(i));
  }
  
  //Data reload
  UInt size = map->size();
  
  map->clear();
  map->reserve(size);
  
  for(ITER iter = tempMap.begin(); iter != tempMap.end(); ++tempMap)
  {
    map->push_back(*iter);
  }
}

template<typename ITEM>
void
pMapManip<ITEM>::
orderData(PMAP & inMap)
{
  //Data serialization
  typedef typename multiset<ITEM>::iterator ITER;
  multiset<ITEM> tempMap;
  
  for(UInt i=1; i <= inMap.size(); ++i)
  {
    tempMap.insert(inMap(i));
  }
  
  //Data reload
  UInt size = inMap.size();
  
  inMap.clear();
  inMap.reserve(size);
  
  for(ITER iter = tempMap.begin(); iter != tempMap.end(); ++tempMap)
  {
    inMap.push_back(*iter);
  }
}

template<typename ITEM>
void
pMapManip<ITEM>::
orderData(Teuchos::RCP<PMAP> & inMap)
{
  orderData(*inMap);
}

template<typename ITEM>
void
pMapManip<ITEM>::
unionExclusive(sVect<PMAP> & sourceMaps)
{
  assert(mapLoaded);
  
  //Exclusive list creation
  set<ITEM> list;
  
  for(UInt k=1; k <= sourceMaps.size(); ++k)
  {
    for(UInt j=1; j <= sourceMaps(k).size(); ++j)
    {
      list.insert(sourceMaps(k).get(j));
    }
  }
  
  //Map bulding
  typedef typename set<ITEM>::iterator ITERATOR;
  ITERATOR iter    = list.begin();
  ITERATOR iterEnd = list.end();
  
  for(; iter != iterEnd; ++iter)
  {
    map->push_back(*iter);
  }
}

template<typename ITEM>
void
pMapManip<ITEM>::
unionExclusive(sVect<PMAP> & sourceMaps, PMAP & inMap)
{
  //Exclusive list creation
  set<ITEM> list;
  
  for(UInt k=1; k <= sourceMaps.size(); ++k)
  {
    for(UInt j=1; j <= sourceMaps(k).size(); ++j)
    {
      list.insert(sourceMaps(k).get(j));
    }
  }
  
  //Map bulding
  typedef typename set<ITEM>::iterator ITERATOR;
  ITERATOR iter    = list.begin();
  ITERATOR iterEnd = list.end();
  
  for(; iter != iterEnd; ++iter)
  {
    inMap.push_back(*iter);
  }
}

template<typename ITEM>
void
pMapManip<ITEM>::
unionExclusive(sVect<PMAP> & sourceMaps, Teuchos::RCP<PMAP> & inMap)
{
  unionExclusive(sourceMaps,*inMap);
}


//_________________________________________________________________________________________________
// MEMORY FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename ITEM>
size_t
pMapManip<ITEM>::
memSize() const
{
  if(container.size() == 0)
  { return(2 * sizeof(bool)); }
  else
  { 
    return( 2 * sizeof(bool)
            + (container.begin()->memSize()
             + 2 * sizeof(Real*) ) * container.size());
  }
}

#endif
