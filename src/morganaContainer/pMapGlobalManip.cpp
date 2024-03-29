/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "pMapGlobalManip.h"

using Teuchos::Array;
using Teuchos::ArrayView;


//_________________________________________________________________________________________________
// PMAPITEM
//-------------------------------------------------------------------------------------------------
pMapGlobalManip<pMapItem>::
pMapGlobalManip() : commDevLoaded(false)
{
}

pMapGlobalManip<pMapItem>::
pMapGlobalManip(const Teuchos::RCP<const communicator> & CommDev) : commDev(CommDev), commDevLoaded(true)
{
}

pMapGlobalManip<pMapItem>::
pMapGlobalManip(const communicator & CommDev) : commDevLoaded(true)
{
  commDev = Teuchos::rcpFromRef(CommDev);
}

void
pMapGlobalManip<pMapItem>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

void
pMapGlobalManip<pMapItem>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

UInt
pMapGlobalManip<pMapItem>::
sizeG(const Teuchos::RCP<const MAP> & Map) const
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  return(sizeG(*Map));
}

UInt
pMapGlobalManip<pMapItem>::
sizeG(const MAP & Map) const
{
  assert(commDevLoaded);
  
  //Massimo locale
  pMapManip<ITEM> doctor;
  UInt localMax = doctor.getMaxGid(Map);
  
  //Massimo globale
  sVect<UInt> values;
  all_gather(*commDev,localMax,values);
  
  for(UInt i=1; i <= values.size(); ++i)
  { localMax = max(localMax,values(i)); }
  
  return(localMax);
}

bool
pMapGlobalManip<pMapItem>::
check(const Teuchos::RCP<const MAP> & Map) const
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  return(check(*Map));
}

bool
pMapGlobalManip<pMapItem>::
check(const MAP & Map) const
{
  assert(commDevLoaded);
  
  //CopyTheMap
  MAP mapCopy(Map);
  
  //Normal Global indexing
  pMapComm<pMapItem> commSystem(commDev);
  commSystem.vectorNormal(mapCopy);
  
  //FinderStartup
  pMapManip<pMapItem> manipulator(mapCopy);
  manipulator.buildFinder();
  
  //Checking
  UInt pid    = commDev->rank();
  UInt pidMax = commDev->size();
  UInt maxGid = sizeG(mapCopy);
  UInt size   = ceil(Real(maxGid) / Real(pidMax));
  UInt gid_l  = pid * size + 1;
  UInt gid_u  = min((pid+1) * size, maxGid);
  
  bool flagGid = true;
  pMapItem item;
  
  for(UInt i=gid_l; i <= gid_u; ++i)
  {
    item.setGid(i);
    flagGid = flagGid & manipulator.isItem(item); 
  }
  
  if(!flagGid)
  { cout << "ERROR: Map with not contiguous gids. Pid: " << pid << endl; }
  
  //Check the pids
  bool flagPid = true;
  
  for(UInt i=1; i <= Map.size(); ++i)
  {
    flagPid = flagPid & (Map.get(i).getPid() == pid);
  }
  
  if(!flagPid)
  { cout << "ERROR: Map pid error. Pid: " << pid << endl; }
  
  //Global value reconstruction  
  UInt flagTemp = flagGid && flagPid;
  bool flag     = flagTemp;
  sVect<UInt> listTemp;
  
  all_gather(*commDev,flagTemp,listTemp);
  
  for(UInt i=1; i <= listTemp.size(); ++i)
  { flag = flag & bool(listTemp(i)); }
  
  return(flag);
}

void
pMapGlobalManip<pMapItem>::
destroyOverlap(Teuchos::RCP<MAP> & Map) const
{
  //Assert and copy
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  
  destroyOverlap(*Map);
}
    
void
pMapGlobalManip<pMapItem>::
destroyOverlap(MAP & Map) const
{
  //Assert and copy
  assert(check(Map));
  assert(commDevLoaded);
  MAP mapCopy(Map);
  
  //Normal communication
  UInt mapSize = sizeG(mapCopy);
  
  pMapComm<pMapItem> mapComm(commDev);
  mapComm.vectorNormal(mapCopy,mapSize);
  
  
  //Create the new map
  typedef typename set<ITEM>::iterator ITER;
  typedef pair<ITER,bool>              PAIR;

  MAP newMap;
  PAIR pair;
  ITER iter;
  set<ITEM> setMap;
  
  for(UInt i=1; i <= mapCopy.size(); ++i)
  {
    pair = setMap.insert(mapCopy(i));
    iter = pair.first;
    
    if( (!pair.second) && (mapCopy(i).getPid() < iter->getPid()) )
    {
      iter->setPid(mapCopy(i).getPid());
    }
  }
  
  newMap.reserve(setMap.size());
  
  for(ITER setIter = setMap.begin(); setIter != setMap.end(); ++setIter)
  {
    newMap.push_back(*setIter);
  }
    
  //Communicate back to the processors
  mapComm.vectorPid(newMap);
  
  pMapManip<pMapItem> serialManip;
  serialManip.setNormalIndexing(newMap);
  
  //Reset the vector
  Map = newMap;
}

void
pMapGlobalManip<pMapItem>::
exportEpetraMap(const Teuchos::RCP<const MAP> & Map,
                Teuchos::RCP<Epetra_Map>      & epetraMap,
                Epetra_MpiComm                & epetraComm,
                const UInt                    & base) const
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  exportEpetraMap(*Map,epetraMap,epetraComm,base);
}

void
pMapGlobalManip<pMapItem>::
exportEpetraMap(const MAP                & Map,
                Teuchos::RCP<Epetra_Map> & epetraMap,
                Epetra_MpiComm           & epetraComm,
                const UInt               & base) const
{
  assert(commDevLoaded);
  
  //Alloc
  int numMyElements = Map.size();
  int MyGlobalElements[numMyElements];
  
  //Renumbering
  for(int i=1; i <= numMyElements; ++i)
  {
    MyGlobalElements[i-1] = Map.get(i).getGid() -1 + base;
  }
  
  //Build map
  epetraMap = Teuchos::rcp(new Epetra_Map(-1,numMyElements,MyGlobalElements,base,epetraComm));
}

void
pMapGlobalManip<pMapItem>::
exportTpetraMap(const Teuchos::RCP<const MAP>  & Map,
                Teuchos::RCP<const TPETRA_MAP> & tpetraMap,
                const UInt                     & base) const
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  exportTpetraMap(*Map,tpetraMap,base);
}
    
void
pMapGlobalManip<pMapItem>::
exportTpetraMap(const MAP                      & Map,
                Teuchos::RCP<const TPETRA_MAP> & tpetraMap,
                const UInt                     & base) const
{
  assert(commDevLoaded);
  
  const Teuchos::RCP<const TPETRA_COMM> tpetraComm = Teuchos::rcp(new TPETRA_MPICOMM(*commDev));
  
  std::vector<TPETRA_GLOBAL_ORDINAL> tpetraGids(Map.size());
  for(UInt i=1; i <= Map.size(); ++i)
  { tpetraGids[i-1] = Map(i).getGid() - 1; }
  
  Teuchos::ArrayView<const TPETRA_GLOBAL_ORDINAL> teuchosArray(tpetraGids);
  
  
  TPETRA_GLOBAL_TYPE indexBase = base;
  TPETRA_GLOBAL_TYPE globSize = Teuchos::OrdinalTraits<TPETRA_GLOBAL_TYPE>::invalid();
  tpetraMap = Teuchos::RCP<TPETRA_MAP>(new TPETRA_MAP(globSize,
                                                      teuchosArray,
                                                      indexBase,
                                                      tpetraComm));
}

void
pMapGlobalManip<pMapItem>::
importEpetraMap(const Teuchos::RCP<const Epetra_Map> & EpetraMap,
                Teuchos::RCP<MAP>                    & Map)
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  
  importEpetraMap(*EpetraMap,Map);
}

void
pMapGlobalManip<pMapItem>::
importEpetraMap(const Epetra_Map  & EpetraMap,
                Teuchos::RCP<MAP> & Map)
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  
  //Alloc
  ITEM item;
  UInt  numMyElements    = EpetraMap.NumMyElements();
  int * MyGlobalElements = EpetraMap.MyGlobalElements();
  int   base             = EpetraMap.IndexBase();
  
  //Clearing
  Map->clear();
  Map->reserve(numMyElements);
  
  //Renumbering
  for(UInt i=1; i <= numMyElements; ++i)
  {
    item.setPid(commDev->rank());
    item.setLid(i);
    item.setGid(MyGlobalElements[i-1] +1 - base);
    
    Map->push_back(item);
  }
}

void 
pMapGlobalManip<pMapItem>::
importTpetraMap(const Teuchos::RCP<const TPETRA_MAP> & TpetraMap,
                Teuchos::RCP<MAP>                    & Map)
{
  assert(commDevLoaded);
  assert(Map.total_count()       > 0);
  assert(TpetraMap.total_count() > 0);
  
  //Clearing
  ITEM item;
  Map->clear();
  Map->reserve(TpetraMap->getNodeNumElements());
  
  //Downloading
  Teuchos::ArrayView<const TPETRA_GLOBAL_ORDINAL> newArray = TpetraMap->getNodeElementList();
  
  for(UInt i=1; i <= newArray.size(); ++i)
  {
    item.setPid(commDev->rank());
    item.setLid(i);
    item.setGid(newArray[i-1] + 1);
    
    Map->push_back(item);
  }
}

void
pMapGlobalManip<pMapItem>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MAP>          & OldMap,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MAP>                & NewMap)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldMap,
                    *NewCommDev,
                    *NewMap);
}

void
pMapGlobalManip<pMapItem>::
reduceCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MAP          & OldMap,
                   const communicator & NewCommDev,
                         MAP          & NewMap)
{
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() < OldCommDev.size());
    
    //Alloc
    pMapItem mItem;
    
    NewMap.clear();
    NewMap.resize(OldMap.size());
    
    //Main loop
    for(UInt i=1; i <= OldMap.size(); ++i)
    {
      mItem = OldMap(i);
      mItem.setPid(NewCommDev.rank());
      NewMap.set(i,mItem);
    }
  }
}

void
pMapGlobalManip<pMapItem>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MAP>          & OldMap,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MAP>                & NewMap)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldMap,
                    *NewCommDev,
                    *NewMap);
}
    
void
pMapGlobalManip<pMapItem>::
expandCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MAP          & OldMap,
                   const communicator & NewCommDev,
                         MAP          & NewMap)
{
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() > OldCommDev.size());
    
    //Alloc
    pMapItem mItem;
    
    NewMap.clear();
    NewMap.resize(OldMap.size());
    
    //Main loop
    for(UInt i=1; i <= OldMap.size(); ++i)
    {
      mItem = OldMap(i);
      mItem.setPid(NewCommDev.rank());
      NewMap.set(i,mItem);
    }
  }
}

size_t
pMapGlobalManip<pMapItem>::
memSize() const
{
  return(sizeof(double));
}


//_________________________________________________________________________________________________
// PMAPITEMSHARE
//-------------------------------------------------------------------------------------------------
pMapGlobalManip<pMapItemShare>::
pMapGlobalManip() : commDevLoaded(false)
{
}

pMapGlobalManip<pMapItemShare>::
pMapGlobalManip(const Teuchos::RCP<const communicator> & CommDev) : commDevLoaded(true), commDev(CommDev)
{
}

pMapGlobalManip<pMapItemShare>::
pMapGlobalManip(const communicator & CommDev) : commDevLoaded(true)
{
  commDev = Teuchos::rcpFromRef(CommDev);
}

void
pMapGlobalManip<pMapItemShare>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

void
pMapGlobalManip<pMapItemShare>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

UInt
pMapGlobalManip<pMapItemShare>::
sizeG(const Teuchos::RCP<const MAP> & Map) const
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  return(sizeG(*Map));
}

UInt
pMapGlobalManip<pMapItemShare>::
sizeG(const MAP & Map) const
{
  assert(commDevLoaded);

  //Massimo locale
  pMapManip<ITEM> doctor;
  UInt localMax = doctor.getMaxGid(Map);
  
  //Massimo globale
  sVect<UInt> values;
  all_gather(*commDev,localMax,values);
  
  for(UInt i=1; i <= values.size(); ++i)
  { localMax = max(localMax,values(i)); }
  
  return(localMax);
}

UInt
pMapGlobalManip<pMapItemShare>::
sharedL(const Teuchos::RCP<const MAP> & Map) const
{
  assert(Map.total_count() > 0);
  return(sharedL(*Map));
}

UInt
pMapGlobalManip<pMapItemShare>::
sharedL(const MAP & Map) const
{
  UInt tot = 0;
  
  for(UInt i=1; i <= Map.size(); ++i)
  {
    if(Map(i).getShared())
    { ++tot;}
  }
  
  return(tot);
}

UInt
pMapGlobalManip<pMapItemShare>::
sharedG(const Teuchos::RCP<const MAP> & Map) const
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  return(sharedG(*Map));
}

UInt
pMapGlobalManip<pMapItemShare>::
sharedG(const MAP & Map) const
{
  assert(commDevLoaded);
  
  //CopyTheMap
  MAP copyMap(Map);
  
  //Alloc
  UInt sum = 0;
  sVect<UInt>   values;
  set<pMapItemShare> container;
  
  //Normal Global indexing
  pMapComm<pMapItemShare> commSystem(commDev);
  commSystem.vectorNormal(copyMap);
  
  //Eliminazione doppioni
  for(UInt i=1; i <= copyMap.size(); ++i)
  {
    if(copyMap(i).getShared())
    { container.insert(copyMap(i)); }
  }
  
  //Compute the sum 
  UInt n = container.size();
  all_gather(*commDev,n,values);
  
  for(UInt i=1; i <= values.size(); ++i)
  { sum += values(i); }
  
  return(sum);
}

UInt
pMapGlobalManip<pMapItemShare>::
ownedL(const Teuchos::RCP<const MAP> & Map) const
{
  assert(Map.total_count() > 0);
  return(ownedL(*Map));
}

UInt
pMapGlobalManip<pMapItemShare>::
ownedL(const MAP & Map) const
{
  UInt tot = 0;
  
  for(UInt i=1; i <= Map.size(); ++i)
  {
    if(Map(i).getOwned())
    { ++tot;}
  }
  
  return(tot);
}

bool
pMapGlobalManip<pMapItemShare>::
check(const Teuchos::RCP<const MAP> & Map) const
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  return(check(*Map));
}

bool
pMapGlobalManip<pMapItemShare>::
check(const MAP & Map) const
{
  assert(commDevLoaded);
  
  //CopyTheMap
  MAP mapCopy(Map);
  
  //Normal Global indexing
  pMapComm<pMapItemShare> commSystem(commDev);
  commSystem.vectorNormal(mapCopy);
  
  //FinderStartup
  pMapManip<pMapItemShare> manipulator(mapCopy);
  manipulator.buildFinder();
  
  //Checking for map completeness------------------------------------
  UInt pid    = commDev->rank();
  UInt pidMax = commDev->size();
  UInt maxGid = sizeG(mapCopy);
  UInt size   = ceil(Real(maxGid) / Real(pidMax));
  UInt gid_l  = pid * size + 1;
  UInt gid_u  = min((pid+1) * size, maxGid);
  
  bool flagGid = true;
  pMapItemShare item;
  
  for(UInt i=gid_l; i <= gid_u; ++i)
  {
    item.setGid(i);
    flagGid = flagGid && manipulator.isItem(item);
  }
  
  if(!flagGid)
  { cout << "ERROR: Map with not contiguous gids. Pid: " << pid << endl; }
  
  //Check the pids---------------------------------------------------
  bool flagPid = true;
  
  for(UInt i=1; i <= Map.size(); ++i)
  {
    flagPid = flagPid & (Map.get(i).getPid() == pid);
  }
  
  if(!flagPid)
  { cout << "ERROR: Map pid error. Pid: " << pid << endl; }
  
  
  //Checking for single ownership------------------------------------
  typedef set<pMapItemShare>::iterator  ITERATOR;
  typedef pair<ITERATOR,bool>           PAIR;
  
  set<pMapItemShare> setList;
  PAIR               output;
  bool flagSSO = true;
  
  for(UInt i=1; i <= mapCopy.size(); ++i)
  {
    if(mapCopy(i).getOwned())
    {
      output  = setList.insert(mapCopy(i));
      flagSSO = flagSSO && output.second;
    }
  }
  
  if(!flagSSO)
  { cout << "ERROR: Map with multiple ownerships. Pid: " << pid << endl; }
  
  //Standard numbering-----------------------------------------------
  ITERATOR iter    = setList.begin();
  ITERATOR iterEnd = setList.end();
  
  UInt index = 1;
  sVect<bool> owningFlags;
  
  owningFlags.reserve(setList.size());
  
  for(; iter != iterEnd; ++iter)
  {
    iter->setLid(index);
    assert(iter->getLid() == index);
    
    owningFlags.push_back(iter->getShared());
    index++;
  }
  
  //Checking that not owned are shared-------------------------------
  sVect<bool> notOwningFlags(owningFlags.size());
  for(UInt i=1; i <= notOwningFlags.size(); ++i)
  {notOwningFlags(i) = false;}
  
  bool flagONS = true;
  
  UInt lid;
  for(UInt i=1; i <= mapCopy.size(); ++i)
  {
    if(!mapCopy(i).getOwned())
    {
      //NotOwned must be shared
      flagONS = flagONS && mapCopy(i).getShared();
      
      //NotOwning vector creation
      output = setList.insert(mapCopy(i));
      assert(!output.second);
      
      item = *(output.first);
      lid  = item.getLid();
      notOwningFlags(lid) = true;
    }
  }
  
  if(!flagONS)
  { cout << "ERROR: Map with not-owned items that are not-shared. Pid: " << pid << endl; }
  
  //Checking that all the owned-and-shared are associated to other notOwned-and-sharedG and viceversa
  assert(owningFlags.size() == notOwningFlags.size());
  
  bool flagS2S = true;
  
  for(UInt i=1; i <= owningFlags.size(); ++i)
  {
    flagS2S = flagS2S && (owningFlags(i) == notOwningFlags(i));
  }
  
  if(!flagS2S)
  { cout << "ERROR: Map with bad association. Shared-Owned items should be associated to other Shared-Not-Owned items. Pid: " << pid << endl; }
  
  //Global value reconstruction 
  UInt flagTemp = flagGid & flagPid & flagSSO & flagONS & flagS2S;
  bool flag     = flagTemp;
  sVect<UInt> listTemp;
  
  all_gather(*commDev,flagTemp,listTemp);
  
  for(UInt i=1; i <= listTemp.size(); ++i)
  { flag = flag & bool(listTemp(i)); }
  
  return(flag);
}

void
pMapGlobalManip<pMapItemShare>::
destroyOverlap(Teuchos::RCP<MAP> & Map) const
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  destroyOverlap(*Map);
}
    
void
pMapGlobalManip<pMapItemShare>::
destroyOverlap(MAP & Map) const
{
  //Assert and alloc
  assert(check(Map));
  assert(commDevLoaded);
  
  pMapItemShare item;
  MAP newMap;
  
  //New map 
  for(UInt i=1 ; i <= Map.size(); ++i)
  {
    if(Map(i).getOwned())
    {
      item = Map(i);
      item.setShared(false);
      
      newMap.push_back(item);
    }
  }
  
  //Reload
  Map = newMap;
}

void
pMapGlobalManip<pMapItemShare>::
updateOwningSharing(Teuchos::RCP<MAP> & Map) const
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  updateOwningSharing(*Map);
}
    
void
pMapGlobalManip<pMapItemShare>::
updateOwningSharing(MAP & Map) const
{
  assert(commDevLoaded);
  
  //Allocation
  UInt pid = commDev->rank();
  
  //Pid update 
  for(UInt i=1; i <= Map.size(); ++i)
  { Map.get(i).setPid(pid); }
  
  //Lid buffering
  Map.bufferLids();
  
  //CopyTheMap
  MAP mapCopy(Map);
  
  //Normal Global indexing
  pMapComm<pMapItemShare> commSystem(commDev);
  commSystem.vectorNormal(mapCopy);
    
  //Ordering
  typedef multiset<pMapItemShare>::iterator ITERATOR;
  
  UInt k=1;
  ITERATOR  iter;
  multiset<pMapItemShare> orderedMap;
  
  for(UInt i=1; i <= mapCopy.size(); ++i)
  { orderedMap.insert(mapCopy(i)); }
  
  //Fix ownership
  typedef set<UInt>::iterator ITERATOR_UINT;
  typedef pair<ITERATOR,ITERATOR> RANGE;
  
  UInt nItems, owningCount;
  pMapItemShare item;
  
  set<UInt> uniqueGids;
  ITERATOR_UINT iterUInt;
  RANGE range;
  
  for(iter = orderedMap.begin(); iter != orderedMap.end(); ++iter) //Compute a unique list of gids
  {
     uniqueGids.insert(iter->getGid());
  }
  
  for(iterUInt = uniqueGids.begin(); iterUInt != uniqueGids.end(); ++iterUInt)
  {
    item.setGid(*iterUInt);
    
    nItems = orderedMap.count(item);    
    assert(nItems >= 1);
    
    range = orderedMap.equal_range(item);
    
    owningCount = 0;
    
    for(iter = range.first; iter != range.second; ++iter)
    {
      owningCount += iter->getOwned();
    }
    
    if(owningCount == 0) //No ownership
    {
      iter = range.first;
      iter->setOwned(true);
    }
    
    if(owningCount >= 2) //Multiple ownerships
    {
      for(iter = range.first; iter != range.second; ++iter)
      {
        iter->setOwned(false);
      }
      
      iter = range.first;
      iter->setOwned(true);
    }
  }
  
  //Serialization onto mapCopy
  for(iter = orderedMap.begin(); iter != orderedMap.end(); ++iter)
  {
    mapCopy(k) = *iter;
    ++k;
  }
  
  //Fix Sharing
  bool shared = false;
  
  for(UInt i=1; i <= mapCopy.size(); ++i)
  {
    //Sharing determination
    shared = false;
    
    //Sharing determination
    if(i != mapCopy.size())
    {
      if( !(mapCopy(i) != mapCopy(i+1)) ) { shared = true; }
    }
    
    if(i != 1)
    {
      if( !(mapCopy(i) != mapCopy(i-1)) ) { shared = true; }
    }
    
    //Setting
    mapCopy(i).setShared(shared);
  }
  
  //Pid map communication
  commSystem.vectorPid(mapCopy);
  
  //Final download
  Map = mapCopy;
  Map.restoreLids();
  
  //Reordering
  pMapManip<pMapItemShare> reorder;
  reorder.setIndexing(Map);
}

void
pMapGlobalManip<pMapItemShare>::
updateSharing(Teuchos::RCP<MAP> & Map) const
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  updateSharing(*Map);
}
    
void
pMapGlobalManip<pMapItemShare>::
updateSharing(MAP & Map) const
{
  assert(commDevLoaded);
  
  //Allocation
  UInt pid = commDev->rank();
  
  //Pid update 
  for(UInt i=1; i <= Map.size(); ++i)
  { Map.get(i).setPid(pid); }
  
  //Lid buffering
  Map.bufferLids();
  
  //CopyTheMap
  MAP mapCopy(Map);
  
  //Normal Global indexing
  pMapComm<pMapItemShare> commSystem(commDev);
  commSystem.vectorNormal(mapCopy);
  
  //Ordering
  typedef multiset<pMapItemShare>::iterator ITERATOR;
  
  UInt k=1;
  ITERATOR  iter;
  multiset<pMapItemShare> orderedMap;
  
  for(UInt i=1; i <= mapCopy.size(); ++i)
  { orderedMap.insert(mapCopy(i)); }
  
  for(iter = orderedMap.begin(); iter != orderedMap.end(); ++iter)
  {
    mapCopy(k) = *iter;
    ++k;
  }
  
  //Sharing/Owning determination
  bool shared = false;
  
  for(UInt i=1; i <= mapCopy.size(); ++i)
  {
    shared = false;
    
    //Sharing determination
    if(i != mapCopy.size())
    {
      if( !(mapCopy(i) != mapCopy(i+1)) ) { shared = true; }
    }
    
    if(i != 1)
    {
      if( !(mapCopy(i) != mapCopy(i-1)) ) { shared = true; }
    }
    
    //Setting
    mapCopy(i).setShared(shared);
  }
  
  //Pid map communication
  commSystem.vectorPid(mapCopy);
  
  //Final download
  Map = mapCopy;
  Map.restoreLids();
  
  //Reordering
  pMapManip<pMapItemShare> reorder;
  reorder.setIndexing(Map);
}

void
pMapGlobalManip<pMapItemShare>::
createSendRecvMap(const Teuchos::RCP<const MAP> & Map,
                  Teuchos::RCP<SENDRECV>        & mapSend,
                  Teuchos::RCP<SENDRECV>        & mapRecv) const
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  assert(mapSend.total_count() > 0);
  assert(mapRecv.total_count() > 0);
  
  createSendRecvMap(*Map,*mapSend,*mapRecv);
}

void
pMapGlobalManip<pMapItemShare>::
createSendRecvMap(const MAP & Map,
                  SENDRECV  & mapSend,
                  SENDRECV  & mapRecv) const
{
  assert(commDevLoaded);
  
  //Clearing data
  mapSend.clear();
  mapRecv.clear();
  
  //Bulding reduced map of only shared objects-----------------------
  MAP mapShared;
  
  for(UInt i=1; i <= Map.size(); ++i)
  {
    if(Map.get(i).getShared())
    {
      mapShared.push_back(Map.get(i));
    }
  }
  
  //Pid designation--------------------------------------------------
  for(UInt i=1; i <= mapShared.size(); ++i)
  {
    mapShared(i).setPid(commDev->rank());
  }
  
  //Normal Global indexing
  pMapComm<pMapItemShare> commSystem(commDev);
  commSystem.vectorNormal(mapShared);
  
  //Owned list-------------------------------------------------------
  typedef set<pMapItemShare>::iterator  ITERATOR;
  typedef pair<ITERATOR,bool>           PAIR;
  
  set<pMapItemShare> setList;
  PAIR               output;
  
  for(UInt i=1; i <= mapShared.size(); ++i)
  {    
    if(mapShared(i).getOwned())
    {
      output = setList.insert(mapShared(i));
      assert(output.second);
    }
  }
  
  //Maps creation----------------------------------------------------
  pMapItemShare    owner,   shared;
  pMapItemSendRecv sender, reciver;
  
  for(UInt i=1; i <= mapShared.size(); ++i)
  {
    if(!mapShared(i).getOwned())
    {
      output = setList.insert(mapShared(i));
      assert(!output.second);
      
      owner  = *(output.first);
      shared = mapShared(i);
      
      assert(owner.getGid() == shared.getGid());
      
      sender.setPid(owner.getPid());
      sender.setLid(owner.getLid());
      sender.setGid(owner.getGid());
      sender.setSid(owner.getPid());
      sender.setRid(shared.getPid());
      mapSend.push_back(sender);
      
      reciver.setPid(shared.getPid());
      reciver.setLid(shared.getLid());
      reciver.setGid(shared.getGid());
      reciver.setSid(owner.getPid());
      reciver.setRid(shared.getPid());
      mapRecv.push_back(reciver);
    }
  }
  
  //Maps distribution------------------------------------------------
  pMapComm<pMapItemSendRecv> mapComm(commDev);
  
  mapComm.vectorPid(mapSend);
  mapComm.vectorPid(mapRecv);
}

void
pMapGlobalManip<pMapItemShare>::
exportEpetraMap(const Teuchos::RCP<const MAP> & Map,
                Teuchos::RCP<Epetra_Map>      & epetraMap,
                Epetra_MpiComm                & epetraComm,
                const UInt                    & base) const
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  exportEpetraMap(*Map,epetraMap,epetraComm,base);
}

void
pMapGlobalManip<pMapItemShare>::
exportEpetraMap(const MAP                & Map,
                Teuchos::RCP<Epetra_Map> & epetraMap,
                Epetra_MpiComm           & epetraComm,
                const UInt               & base) const
{
  assert(commDevLoaded);
  
  //Alloc
  int numMyElements      = Map.size();
  int MyGlobalElements[numMyElements];
  
  //Renumbering
  for(int i=1; i <= numMyElements; ++i)
  { MyGlobalElements[i-1] = Map.get(i).getGid() -1 + base; }
  
  //Build map
  epetraMap = Teuchos::rcp(new Epetra_Map(-1,numMyElements,MyGlobalElements,base,epetraComm));
}

void
pMapGlobalManip<pMapItemShare>::
exportTpetraMap(const Teuchos::RCP<const MAP>  & Map,
                Teuchos::RCP<const TPETRA_MAP> & tpetraMap,
                const UInt                     & base) const
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  exportTpetraMap(*Map,tpetraMap,base);
}
    
void
pMapGlobalManip<pMapItemShare>::
exportTpetraMap(const MAP                      & Map,
                Teuchos::RCP<const TPETRA_MAP> & tpetraMap,
                const UInt                     & base) const
{  
  assert(commDevLoaded);
  
  const Teuchos::RCP<const TPETRA_COMM> tpetraComm = Teuchos::rcp(new TPETRA_MPICOMM(*commDev));
  
  std::vector<TPETRA_GLOBAL_ORDINAL> tpetraGids(Map.size());
  for(UInt i=1; i <= Map.size(); ++i)
  { tpetraGids[i-1] = Map(i).getGid() - 1; }
  
  Teuchos::ArrayView<const TPETRA_GLOBAL_ORDINAL> teuchosArray(tpetraGids);
  
  
  TPETRA_GLOBAL_TYPE indexBase = base;
  TPETRA_GLOBAL_TYPE globSize = Teuchos::OrdinalTraits<TPETRA_GLOBAL_TYPE>::invalid();
  tpetraMap = Teuchos::RCP<TPETRA_MAP>(new TPETRA_MAP(globSize,
                                                      teuchosArray,
                                                      indexBase,
                                                      tpetraComm));
}

void
pMapGlobalManip<pMapItemShare>::
importEpetraMap(const Teuchos::RCP<const Epetra_Map> & EpetraMap,
                Teuchos::RCP<MAP>                    & Map)
{
  assert(commDevLoaded);
  assert(Map.total_count() > 0);
  importEpetraMap(*EpetraMap,Map);
}

void
pMapGlobalManip<pMapItemShare>::
importEpetraMap(const Epetra_Map  & EpetraMap,
                Teuchos::RCP<MAP> & Map)
{
  assert(commDevLoaded);
  
  //Alloc
  ITEM item;
  UInt  numMyElements    = EpetraMap.NumMyElements();
  int * MyGlobalElements = EpetraMap.MyGlobalElements();
  int   base             = EpetraMap.IndexBase();
  
  //Clearing
  Map->clear();
  Map->reserve(numMyElements);
  
  //Renumbering
  for(UInt i=1; i <= numMyElements; ++i)
  {
    item.setPid(commDev->rank());
    item.setLid(i);
    item.setGid(MyGlobalElements[i-1] +1 - base);
    
    Map->push_back(item);
  }
}

void 
pMapGlobalManip<pMapItemShare>::
importTpetraMap(const Teuchos::RCP<const TPETRA_MAP> & TpetraMap,
                Teuchos::RCP<MAP>                    & Map)
{
  assert(commDevLoaded);
  
  //Clearing
  ITEM item;
  Map->clear();
  Map->reserve(TpetraMap->getNodeNumElements());
  
  //Downloading
  Teuchos::ArrayView<const TPETRA_GLOBAL_ORDINAL> newArray = TpetraMap->getNodeElementList();
  
  for(UInt i=1; i <= newArray.size(); ++i)
  {
    item.setPid(commDev->rank());
    item.setLid(i);
    item.setGid(newArray[i-1] + 1);
    
    Map->push_back(item);
  }
}

void
pMapGlobalManip<pMapItemShare>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MAP>          & OldMap,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MAP>                & NewMap)
{
  reduceCommunicator(isActive,
                     OldCommDev,
                     OldMap,
                     NewCommDev,
                     NewMap);
}

void
pMapGlobalManip<pMapItemShare>::
reduceCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MAP          & OldMap,
                   const communicator & NewCommDev,
                         MAP          & NewMap)
{
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() < OldCommDev.size());
    
    //Alloc
    pMapItemShare mItem;
    
    NewMap.clear();
    NewMap.resize(OldMap.size());
    
    //Main loop
    for(UInt i=1; i <= OldMap.size(); ++i)
    {
      mItem = OldMap(i);
      mItem.setPid(NewCommDev.rank());
      NewMap.set(i,mItem);
    }
  }
}

void
pMapGlobalManip<pMapItemShare>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MAP>          & OldMap,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MAP>                & NewMap)
{
  expandCommunicator(isActive,
                     OldCommDev,
                     OldMap,
                     NewCommDev,
                     NewMap);
}
    
void
pMapGlobalManip<pMapItemShare>::
expandCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MAP          & OldMap,
                   const communicator & NewCommDev,
                         MAP          & NewMap)
{
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() > OldCommDev.size());
    
    //Alloc
    pMapItemShare mItem;
    
    NewMap.clear();
    NewMap.resize(OldMap.size());
    
    //Main loop
    for(UInt i=1; i <= OldMap.size(); ++i)
    {
      mItem = OldMap(i);
      mItem.setPid(NewCommDev.rank());
      NewMap.set(i,mItem);
    }
  }
}

size_t
pMapGlobalManip<pMapItemShare>::
memSize() const
{
  return(sizeof(double));
}


//_________________________________________________________________________________________________
// PMAPITEMSENDRECV
//-------------------------------------------------------------------------------------------------
pMapGlobalManip<pMapItemSendRecv>::
pMapGlobalManip() : commDevLoaded(false)
{
}

pMapGlobalManip<pMapItemSendRecv>::
pMapGlobalManip(const Teuchos::RCP<const communicator> & CommDev) : commDevLoaded(true), commDev(CommDev)
{
}

pMapGlobalManip<pMapItemSendRecv>::
pMapGlobalManip(const communicator & CommDev) : commDevLoaded(true)
{
  commDev = Teuchos::rcpFromRef(CommDev);
}

void
pMapGlobalManip<pMapItemSendRecv>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

void
pMapGlobalManip<pMapItemSendRecv>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

bool
pMapGlobalManip<pMapItemSendRecv>::
check(const Teuchos::RCP<const MAP> & mapSend, const Teuchos::RCP<const MAP> & mapRecv) const
{
  assert(commDevLoaded);
  assert(mapSend.total_count() > 0);
  assert(mapRecv.total_count() > 0);
  return(check(*mapSend,*mapRecv));
}

bool
pMapGlobalManip<pMapItemSendRecv>::
check(const MAP & mapSend, const MAP & mapRecv) const
{
  assert(commDevLoaded);
  
  //CopyTheMap
  MAP mapSendCopy(mapSend);
  MAP mapRecvCopy(mapRecv);
  
  //Normal Global indexing
  pMapComm<pMapItemSendRecv> commSystem(commDev);  
  commSystem.vectorNormal(mapSendCopy);
  commSystem.vectorNormal(mapRecvCopy);
  
  //Find matching
  typedef set<pMapItemSendRecv>::iterator  ITERATOR;
  typedef pair<ITERATOR,bool>              PAIR;
  
  bool  flag = true;
  PAIR  output;
  set<pMapItemSendRecv> setListSend;
  
  //Control
  if(mapSendCopy.size() == mapRecvCopy.size())
  {
    //Sending list
    for(UInt i=1; i <= mapSendCopy.size(); ++i)
    {
      output = setListSend.insert(mapSendCopy(i));
      assert(output.second);
    }
    
    //Matching
    for(UInt i=1; i <= mapRecvCopy.size(); ++i)
    {
      output = setListSend.insert(mapRecvCopy(i));
      flag   = flag & (!output.second);
    }
  }
  else
  { 
    flag = false;
  }
  
  
  //Global value reconstruction  
  UInt flagTemp = flag;
  sVect<UInt> listTemp;
  
  all_gather(*commDev,flagTemp,listTemp);
  
  for(UInt i=1; i <= listTemp.size(); ++i)
  { flag = flag & bool(listTemp(i)); }
  
  return(flag);
}

sVect<UInt>
pMapGlobalManip<pMapItemSendRecv>::
countSend(const Teuchos::RCP<const MAP> & mapSend) const
{
  assert(commDevLoaded);
  assert(mapSend.total_count() > 0);
  return(countSend(*mapSend));
}

sVect<UInt>
pMapGlobalManip<pMapItemSendRecv>::
countSend(const MAP & mapSend) const
{
  assert(commDevLoaded);
  
  sVect<UInt> inCount(commDev->size());
  UInt rid;
  
  for(UInt i=1; i <= mapSend.size(); ++i)
  {
    rid = mapSend(i).getRid();
    assert(rid < UInt(commDev->size()));
    
    inCount[rid]++;
  }
  
  return(inCount);
}
    
sVect<UInt>
pMapGlobalManip<pMapItemSendRecv>::
countRecv(const Teuchos::RCP<const MAP> & mapRecv) const
{
  assert(commDevLoaded);
  assert(mapRecv.total_count() > 0);
  return(countRecv(*mapRecv));
}
    
sVect<UInt>
pMapGlobalManip<pMapItemSendRecv>::
countRecv(const MAP & mapRecv) const
{
  assert(commDevLoaded);
  
  sVect<UInt> outCount(commDev->size());
  UInt sid;
  
  for(UInt i=1; i <= mapRecv.size(); ++i)
  {
    sid = mapRecv.get(i).getSid();
    assert(sid < UInt(commDev->size()));
    
    outCount[sid]++;
  }
  
  return(outCount);
}

void
pMapGlobalManip<pMapItemSendRecv>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MAP>          & OldMap,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MAP>                & NewMap)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldMap,
                    *NewCommDev,
                    *NewMap);
}

void
pMapGlobalManip<pMapItemSendRecv>::
reduceCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MAP          & OldMap,
                   const communicator & NewCommDev,
                         MAP          & NewMap)
{
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() < OldCommDev.size());
    
    //Alloc
    pMapItemSendRecv mItem;
    
    NewMap.clear();
    NewMap.resize(OldMap.size());
    
    //Main loop
    for(UInt i=1; i <= OldMap.size(); ++i)
    {
      mItem = OldMap(i);
      mItem.setPid(NewCommDev.rank());
      NewMap.set(i,mItem);
    }
  }
}

void
pMapGlobalManip<pMapItemSendRecv>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const MAP>          & OldMap,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<MAP>                & NewMap)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldMap,
                    *NewCommDev,
                    *NewMap);
}
    
void
pMapGlobalManip<pMapItemSendRecv>::
expandCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const MAP          & OldMap,
                   const communicator & NewCommDev,
                         MAP          & NewMap)
{
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() > OldCommDev.size());
    
    //Alloc
    pMapItemSendRecv mItem;
    
    NewMap.clear();
    NewMap.resize(OldMap.size());
    
    //Main loop
    for(UInt i=1; i <= OldMap.size(); ++i)
    {
      mItem = OldMap(i);
      mItem.setPid(NewCommDev.rank());
      NewMap.set(i,mItem);
    }
  }
}

size_t
pMapGlobalManip<pMapItemSendRecv>::
memSize() const
{
  return(sizeof(bool));
}
