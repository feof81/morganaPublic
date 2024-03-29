/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PVECTGLOBALMANIP_HPP
#define PVECTGLOBALMANIP_HPP

#include "Teuchos_RCP.hpp"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "pMapItem.h"
#include "pMapItemShare.h"
#include "pMapGlobalManip.h"
#include "pVectManip.hpp"
#include "pVectComm.hpp"

#include "traitsBasic.h"
#include "traitsMapItemFixer.hpp"


using namespace std;
using namespace boost::mpi;


//_________________________________________________________________________________________________
// UNSPECIALIZED
//-------------------------------------------------------------------------------------------------

/*! Unspecialized - empty */
template<typename DATA, typename MAPITEM> class pVectGlobalManip
{ };


//_________________________________________________________________________________________________
// PMAPITEM
//-------------------------------------------------------------------------------------------------

/*! Performs global manipulations, retrives informations and checks the \c pVect. Specialized version for \c pMapItem */
template<typename DATA> class pVectGlobalManip<DATA,pMapItem>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pVect<DATA,pMapItem>   PVECT;
    typedef sVect<DATA>            DATAVECT;
    typedef pMap<pMapItem>         PMAP;
    typedef pMap<pMapItemSendRecv> SENDRECV;
    //@}
    
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
    
    
    /*! @name Constructors and set functions */ //@{
  public:
    pVectGlobalManip();
    pVectGlobalManip(const Teuchos::RCP<const communicator> & CommDev);
    pVectGlobalManip(const communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    //@}
    
    
    /*! @name Data retriving functions */ //@{
  public:
    /*! Returns the maximum \c gid on all the processors */
    UInt sizeG(const Teuchos::RCP<const PVECT> & Vect) const;
    
    /*! Returns the maximum \c gid on all the processors */
    UInt sizeG(const PVECT & Vect) const;
    
    /*! Performs a checks on the \c pMap contained into \c Vect, see \c pMapGlobalManip */
    bool check(const Teuchos::RCP<const PVECT> & Vect) const;
    
    /*! Performs a checks on the \c pMap contained into \c Vect, see \c pMapGlobalManip */
    bool check(const PVECT & Vect) const;
    
    /*! Checks that items with the same \c gid have the same data */
    bool checkConsistency(const Teuchos::RCP<const PVECT> & Vect) const;
    
    /*! Checks that items with the same \c gid have the same data */
    bool checkConsistency(const PVECT & Vect) const;
    
    /*! Computes the minimum and maximum data */
    void dataMinMax(const Teuchos::RCP<const PVECT> & Vect,
                                               DATA & DataMin,
                                               DATA & DataMax);
    
    /*! Computes the minimum and maximum data */
    void dataMinMax(const PVECT & Vect,
                           DATA & DataMin,
                           DATA & DataMax);
    //@}
    
    
    /*! @name Global manipulations */ //@{
  public:
    /*! Assigns a unique \c gid to data. If two data are equal the same gid is assigned to both.
    Data should implement the inequality and less operator. */
    void buildGlobalNumbering(Teuchos::RCP<PVECT> & Vect) const;
    
    /*! Assigns a unique \c gid to data. If two data are equal the same gid is assigned to both.
    Data should implement the inequality and less operator. */
    void buildGlobalNumbering(PVECT & Vect) const;
    
    /*! A higher performance method of \c buildGlobalNumbering. Local elements are marked with 
     the \c isLocal vector and they are not broadcasted.*/
    void buildGlobalNumbering(      Teuchos::RCP<PVECT>              & Vect,
                              const Teuchos::RCP<const sVect<bool> > & isLocal) const;
    
    /*! A higher performance method of \c buildGlobalNumbering. Local elements are marked with 
     the \c isLocal vector and they are not broadcasted*/
    void buildGlobalNumbering(            PVECT & Vect,
                              const sVect<bool> & isLocal) const;
    
    /*! Creates a new vector whose new map is \c NewMap The data contained in \c Vect are communicated 
    to satisfy the new map*/
    void changeMap(      Teuchos::RCP<PVECT>      & Vect,
                   const Teuchos::RCP<const PMAP> & NewMap) const;
    
    /*! Creates a new vector whose new map is \c NewMap The data contained in \c Vect are communicated 
    to satisfy the new map*/
    void changeMap(     PVECT & Vect,
                   const PMAP & NewMap) const;
    
    /*! Distributes the vector using the \c pid logic. The repeated \c gid are eliminated */
    void vectorNormalExclusive(Teuchos::RCP<PVECT> & Vect);
    
    /*! Distributes the vector using the \c pid logic. The repeated \c gid are eliminated */
    void vectorNormalExclusive(PVECT & Vect);
    
    /*! Linear-epetra data distribution 
    \param Vect contains the data to send, then is cleared and loaded with the new data
    */
    Epetra_Map vectorLinear(Teuchos::RCP<PVECT> & Vect);
    
    /*! Linear-epetra data distribution 
    \param Vect contains the data to send, then is cleared and loaded with the new data
    */
    Epetra_Map vectorLinear(PVECT & Vect);
    
    /*! AllReduce - the repeated elements are reduced */
    template<typename OP>
    void allReduce(Teuchos::RCP<PVECT> & Vect, OP op) const;
    
    /*! AllReduce - the repeated elements are reduced */
    template<typename OP>
    void allReduce(PVECT & Vect, OP op) const;
    //@}
    
    
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copies the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const PVECT>        & OldVect,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<PVECT>              & NewVect);
    
    /*! Copies the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const PVECT        & OldVect,
                            const communicator & NewCommDev,
                                  PVECT        & NewVect);
    
    /*! Copies the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const PVECT>        & OldVect,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<PVECT>              & NewVect);
    
    /*! Copies the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const PVECT        & OldVect,
                            const communicator & NewCommDev,
                                  PVECT        & NewVect);
    //@}
    
    /*! @name Other functions */ //@{
  public:
    size_t memSize() const;
    //@}
};


template<typename DATA>
pVectGlobalManip<DATA,pMapItem>::
pVectGlobalManip() : commDevLoaded(false)
{
}

template<typename DATA>
pVectGlobalManip<DATA,pMapItem>::
pVectGlobalManip(const Teuchos::RCP<const communicator> & CommDev) : commDevLoaded(true), commDev(CommDev)
{
}

template<typename DATA>
pVectGlobalManip<DATA,pMapItem>::
pVectGlobalManip(const communicator & CommDev) : commDevLoaded(true)
{
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename DATA>
UInt
pVectGlobalManip<DATA,pMapItem>::
sizeG(const Teuchos::RCP<const PVECT> & Vect) const
{
  assert(commDevLoaded);
  pMapGlobalManip<pMapItem> manip(commDev);
  return(manip.sizeG(Vect->getMapRef()));
}

template<typename DATA>
UInt
pVectGlobalManip<DATA,pMapItem>::
sizeG(const PVECT & Vect) const
{
  assert(commDevLoaded);
  pMapGlobalManip<pMapItem> manip(commDev);
  return(manip.sizeG(Vect.getMapRef()));
}

template<typename DATA>
bool
pVectGlobalManip<DATA,pMapItem>::
check(const Teuchos::RCP<const PVECT> & Vect) const
{
  assert(commDevLoaded);
  pMapGlobalManip<pMapItem> manip(commDev);
  return(manip.check(Vect->getMapRef()));
}

template<typename DATA>
bool
pVectGlobalManip<DATA,pMapItem>::
check(const PVECT & Vect) const
{
  assert(commDevLoaded);
  pMapGlobalManip<pMapItem> manip(commDev);
  return(manip.check(Vect.getMapRef()));
}

template<typename DATA>
bool
pVectGlobalManip<DATA,pMapItem>::
checkConsistency(const Teuchos::RCP<const PVECT> & Vect) const
{
  //Assert
  assert(commDevLoaded);
  return(checkConsistency(*Vect));
}

template<typename DATA>
bool
pVectGlobalManip<DATA,pMapItem>::
checkConsistency(const PVECT & Vect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Vect copy
  PVECT vectCopy = Vect;
  
  //Data redistribution
  pVectComm<DATA,pMapItem> distributor(commDev);
  distributor.vectorNormal(vectCopy);
  
  //Ordering
  pVectManip<DATA,pMapItem> orderer;
  orderer.orderGid(vectCopy);
  
  //Checking
  bool flag = true;
  
  if(vectCopy.size() > 1)
  {
    for(UInt i=1; i < vectCopy.size(); ++i)
    {
      if(vectCopy.getMapL(i).getGid() == vectCopy.getMapL(i+1).getGid())
      {
	flag = flag & (!(vectCopy.getDataL(i) != vectCopy.getDataL(i+1)));
      }
    }
  }
  
  //Global value reconstruction  
  bool flagGlob = flag;
  sVect<UInt> listTemp;
  
  all_gather(*commDev,UInt(flag),listTemp);
  
  for(UInt i=1; i <= listTemp.size(); ++i)
  { flagGlob = flagGlob & bool(listTemp(i)); }
  
  return(flagGlob);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
dataMinMax(const Teuchos::RCP<const PVECT> & Vect,
                                      DATA & DataMin,
                                      DATA & DataMax)
{
  assert(commDevLoaded);
  dataMinMax(*Vect,DataMin,DataMax);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
dataMinMax(const PVECT & Vect,
                  DATA & DataMin,
                  DATA & DataMax)
{
  assert(commDevLoaded);
  
  //Allocations
  UInt numProc = commDev->size();
  
  //Local maximum and minimum data determination
  DATA dataMax, dataMin;
  dataMax = dataMin;
  
  UInt dataLoaded = (Vect.size() >= 1);
  
  if(dataLoaded)
  {
    dataMin = Vect.getL(1);
    dataMax = Vect.getL(1);
  }
  
  for(UInt i=1; i <= Vect.size(); ++i)
  {
    dataMin = std::min(dataMin,Vect.getL(i));
    dataMax = std::max(dataMax,Vect.getL(i));
  }  
  
  //Global maximum and minimum data determination
  sVect<UInt> allLoading(numProc);
  morgana_all_gather(*commDev,dataLoaded,allLoading);

  sVect<DATA> allMin(numProc);
  morgana_all_gather(*commDev,dataMin,allMin);
  
  sVect<DATA> allMax(numProc);
  morgana_all_gather(*commDev,dataMax,allMax);
  
  for(UInt i=1; i<=numProc; ++i)
  {
    if(allLoading(i))
    {
      dataMin = min(dataMin,allMin(i));
      dataMax = max(dataMax,allMax(i));
    }
  }
  
  //Data assignment
  DataMin = dataMin;
  DataMax = dataMax;
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
buildGlobalNumbering(Teuchos::RCP<PVECT> & Vect) const
{ 
  assert(commDevLoaded);
  buildGlobalNumbering(*Vect);  
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
buildGlobalNumbering(PVECT & Vect) const
{ 
  assert(commDevLoaded);
  
  //Allocations
  UInt numProc = commDev->size();
  UInt proc    = commDev->rank();
  
  //Vect copy
  Vect.bufferLids();
  PVECT vectCopy = Vect;
  
  //Pid buffering
  for(UInt i=1; i<=vectCopy.size(); ++i)
  {
    vectCopy.getMapL(i).setPid(proc);
  }

  //Data redistribution
  pVectComm<DATA,pMapItem> distributor(commDev);
  distributor.vectorData(vectCopy);  
  
  //Ordering
  pVectManip<DATA,pMapItem> orderer;
  orderer.orderData(vectCopy);
  
  //Local numbering
  UInt k = 1;
  
  if(vectCopy.size() == 0) { k = 0; }
  
  for(UInt i=1; i <= vectCopy.size(); ++i)
  {
    vectCopy.getMapL(i).setGid(k);
    
    if(i != vectCopy.size())
    {
      if(vectCopy(i) != vectCopy(i+1))
      { k++; }
    }
  }
  
  //Offset determination
  sVect<UInt> localSizes(numProc);
  all_gather(*commDev,k,localSizes);
  
  UInt offset = 0;
  for(UInt i=0; i < proc; ++i)
  {
    offset += localSizes[i];
  }
  
  //Global numbering
  for(UInt i=1; i <= vectCopy.size(); ++i)
  {
    k = vectCopy.getMapL(i).getGid();
    vectCopy.getMapL(i).setGid(k + offset);
  }
  
  //Pid map communication
  distributor.vectorPid(vectCopy);
  
  //Data download
  Vect = vectCopy;
  Vect.restoreLids();
  
  //Data lid ordering
  pVectManip<DATA,pMapItem> manipulator;
  manipulator.setIndexing(Vect);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
buildGlobalNumbering(      Teuchos::RCP<PVECT>              & Vect,
                     const Teuchos::RCP<const sVect<bool> > & isLocal) const
{
  assert(commDevLoaded);
  buildGlobalNumbering(*Vect,*isLocal);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
buildGlobalNumbering(            PVECT & Vect,
                     const sVect<bool> & isLocal) const
{  
  assert(commDevLoaded);
  
  //Assert
  assert(Vect.size() == isLocal.size());
  
  //Allocations
  UInt numProc = commDev->size();
  UInt proc    = commDev->rank();
  
  //Vector splitting
  sVect<UInt> lLids, gLids;
  PVECT lVect, gVect;
  
  for(UInt i=1; i <= Vect.size(); ++i)
  {
    if(isLocal(i))
    {
      lVect.push_back(Vect.getMapL(i), Vect.getDataL(i));
      lLids.push_back(i);
    }
    else
    {
      gVect.push_back(Vect.getMapL(i), Vect.getDataL(i));
      gLids.push_back(i);
    }
  }  
  
  lVect.updateFinder();
  gVect.updateFinder();  
  
  //Broadcast offsets
  sVect<UInt> localSizes(numProc);
  all_gather(*commDev,lVect.size(), localSizes);
  
  UInt offset    = 0;
  UInt totOffset = 0;
  
  for(UInt i=0; i < proc; ++i)
  { offset += localSizes[i]; }
  
  for(UInt i=0; i < numProc; ++i)
  { totOffset += localSizes[i]; }
  
  //Global vector temp numbering 
  buildGlobalNumbering(gVect);
  
  //Local vector offsetting 
  for(UInt i=1; i <= lVect.size(); ++i)
  {
    lVect.getMapL(i).setLid(lLids(i));
    lVect.getMapL(i).setGid(i + offset);
  }
  
  //Global vector offsetting
  UInt gid;
  
  for(UInt i=1; i <= gVect.size(); i++)
  {
    gid = gVect.getMapL(i).getGid();
    
    gVect.getMapL(i).setLid(gLids(i));
    gVect.getMapL(i).setGid(gid + totOffset);
    
    lVect.push_back(gVect.getMapL(i),gVect.getDataL(i));
  }
  
  //Lid reordering
  pVectManip<DATA,pMapItem> localManip;
  localManip.setIndexing(lVect);
  
  //Copy
  Vect = lVect;
  Vect.updateFinder();
  check(Vect);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
changeMap(      Teuchos::RCP<PVECT>      & Vect,
          const Teuchos::RCP<const PMAP> & NewMap) const
{
  assert(commDevLoaded);
  changeMap(*Vect,*NewMap);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
changeMap(     PVECT & Vect,
          const PMAP & NewMap) const
{
  assert(commDevLoaded);
  
  //Pre-checking
  pMapManip<pMapItem>       preLocal;
  pMapGlobalManip<pMapItem> preGlobal(commDev);
  
  assert(preLocal.checkNormalIndexing(NewMap));
  assert(preGlobal.check(NewMap));
  
  //Processor ids
  UInt proc = commDev->rank();
  
  //Allocations and vector copy
  Vect.bufferLids();
  PVECT vectCopy = Vect;
  PVECT newVect;
  
  //Pid buffering
  for(UInt i=1; i<=vectCopy.size(); ++i)
  { vectCopy.getMapL(i).setPid(proc); }
  
  //Maximum gid determination
  UInt maxGid = sizeG(vectCopy);
  
  //Local items of the vector
  UInt gid;
  PMAP retrivingMap;
  
  for(UInt i=1; i <= NewMap.size(); ++i)
  {
    gid = NewMap.get(i).getGid();
    
    if(Vect.isG(gid))
    { newVect.push_back(Vect.getDataG(gid),NewMap.get(i),false); }
    else
    { retrivingMap.push_back(NewMap.get(i)); }
  }
  
  retrivingMap.bufferLids();
  
  //Data normal indexing 
  pVectComm<DATA,pMapItem> vectDistributor(commDev);
  vectDistributor.vectorNormal(vectCopy);
  
  pMapComm<pMapItem> mapDistributor(commDev);
  mapDistributor.vectorNormal(retrivingMap,maxGid);
  
  retrivingMap.restoreLids();
  vectCopy.restoreLids();
  
  //Build commPattern
  SENDRECV sendMap(retrivingMap.size());
  pMapItemSendRecv sendItem;
  
  for(UInt i=1; i <= retrivingMap.size(); ++i)
  {
    gid = retrivingMap(i).getGid();
    assert(vectCopy.isG(gid));
    
    sendItem.setPid( vectCopy.getMapG(gid).getPid() );
    sendItem.setLid( vectCopy.getMapG(gid).getLid() );
    sendItem.setGid( gid );
    sendItem.setSid( vectCopy.getMapG(gid).getPid() );
    sendItem.setRid( retrivingMap(i).getPid() );
    
    sendMap(i) = sendItem;
  }
  
  //Pid indexing of send-recv
  pMapComm<pMapItemSendRecv> sendDistributor(commDev);
  sendDistributor.vectorPid(sendMap);
  
  //Sending of data-items
  PVECT retrivedVect;
  vectDistributor.sendRecv(sendMap,Vect,retrivedVect);  
  
  //VectMerging
  for(UInt i=1; i <= NewMap.size(); ++i)
  {
    gid = NewMap.get(i).getGid();
    
    if(retrivedVect.isG(gid))
    { newVect.push_back(retrivedVect.getDataG(gid),NewMap.get(i),false); }
  }
  
  //Data reordering  
  pVectManip<DATA,pMapItem> vectManipulator;
  vectManipulator.setIndexing(newVect);
  
  newVect.updateFinder();
  Vect = newVect;
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
vectorNormalExclusive(Teuchos::RCP<PVECT> & Vect)
{
  assert(commDevLoaded);
  vectorNormalExclusive(*Vect);
}
    
template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
vectorNormalExclusive(PVECT & Vect)
{
  assert(commDevLoaded);
  
  //Push to pid distribution
  pVectComm<DATA,pMapItem> distributor(commDev);
  distributor.vectorNormal(Vect);
  
  //Eliminate repeated gids
  pVectManip<DATA,pMapItem> serialManip;
  serialManip.unionExclusive(Vect);
  
  //Set pid
  UInt pid = commDev->rank();
  
  for(UInt i=1; i <= Vect.size(); ++i)
  {
    mapItemFixer(Vect.getMapL(i));
    Vect.getMapL(i).setPid(pid);
  }
  
  //Checking
  assert(check(Vect));
}

template<typename DATA>
Epetra_Map
pVectGlobalManip<DATA,pMapItem>::
vectorLinear(Teuchos::RCP<PVECT> & Vect)
{
  //Assert and pid
  assert(commDevLoaded);
  return(vectorLinear(*Vect));
}

template<typename DATA>
Epetra_Map
pVectGlobalManip<DATA,pMapItem>::
vectorLinear(PVECT & Vect)
{
  //Assert and pid
  assert(commDevLoaded);
  UInt pid = commDev->rank();
  
  //New map epetra
  UInt mapSize = sizeG(Vect);
  Epetra_MpiComm epetraComm(*commDev);
  Epetra_Map     newMap_epetra(int(mapSize), int(0), epetraComm);
  
  //New map morgana
  pMapGlobalManip<pMapItem> mapImporter(commDev);
  Teuchos::RCP<PMAP> newMap_morgana(new PMAP);
  mapImporter.importEpetraMap(newMap_epetra,newMap_morgana);
  
  for(UInt i=1; i <= newMap_morgana->size(); ++i)
  {
    mapItemFixer(newMap_morgana->get(i));
    newMap_morgana->get(i).setPid(pid);
  }
  
  //Redistribute data with the new map
  changeMap(Vect,*newMap_morgana);
  
  //Return
  return(newMap_epetra);
}

template<typename DATA>
template<typename OP>
void
pVectGlobalManip<DATA,pMapItem>::
allReduce(Teuchos::RCP<PVECT> & Vect, OP op) const
{
  allReduce(*Vect,op);
}

template<typename DATA>
template<typename OP>
void
pVectGlobalManip<DATA,pMapItem>::
allReduce(PVECT & Vect, OP op) const
{
  //Assert-----------------------------------------------------------
  assert(commDevLoaded);
  assert(check(Vect));
  
  //Typedef----------------------------------------------------------
  typedef traitsBasic<DATA> TRAITS;
  
  //Normal distribution----------------------------------------------
  pVectComm<DATA,pMapItem> pComm(commDev);
  pComm.vectorNormal(Vect);
  
  //Perform the reduction--------------------------------------------
  UInt gid;
  UInt minGid = std::numeric_limits<UInt>::max();
  UInt maxGid = 0;
  
  static const UInt numI = traitsBasic<DATA>::numI;
  static const UInt numJ = traitsBasic<DATA>::numJ;
  
  for(UInt i=1; i <= Vect.size(); ++i)
  {
    minGid = std::min(minGid, Vect.getMapL(i).getGid());
    maxGid = std::max(maxGid, Vect.getMapL(i).getGid());
  }
  
  sVect<DATA> redVect(maxGid -minGid +1);
  Real valueA, valueB, valueT;
  
  if(Vect.size() != 0)
  {
    gid = Vect.getMapL(1).getGid();
    redVect(gid -minGid +1) = Vect(1);
    
    for(UInt i=2; i <= Vect.size(); ++i)
    {
      gid = Vect.getMapL(i).getGid();
      
      for(UInt I=1; I <= numI; ++I)
      {
        for(UInt J=1; J <= numJ; J++)
        {
          valueA = TRAITS::getIJ(I,J,redVect(gid -minGid +1));
          valueB = TRAITS::getIJ(I,J,Vect(i));
          valueT = op.operator()(valueA,valueB);
          TRAITS::setIJ(I,J,valueT,redVect(gid -minGid +1));
        }
      }
    }
  }
  
  for(UInt i=1; i <= Vect.size(); ++i)
  {
    gid     = Vect.getMapL(i).getGid();
    Vect(i) = redVect(gid -minGid +1);
  }
  
  //Comm back--------------------------------------------------------
  pComm.vectorPid(Vect);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const PVECT>        & OldVect,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<PVECT>              & NewVect)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldVect,
                    *NewCommDev,
                    *NewVect);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
reduceCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const PVECT        & OldVect,
                   const communicator & NewCommDev,
                         PVECT        & NewVect)
{
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() < OldCommDev.size());
    
    //Alloc
    pMapItem mItem;
    
    NewVect.clear();
    NewVect.resize(OldVect.size());
    
    //Main loop
    for(UInt i=1; i <= OldVect.size(); ++i)
    {
      mItem = OldVect.getMapL(i);
      mItem.setPid(NewCommDev.rank());
      
      NewVect.getMapL(i)  = mItem;
      NewVect.getDataL(i) = OldVect.getDataL(i);
    }
    
    NewVect.updateFinder();
  }
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const PVECT>        & OldVect,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<PVECT>              & NewVect)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldVect,
                    *NewCommDev,
                    *NewVect);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItem>::
expandCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const PVECT        & OldVect,
                   const communicator & NewCommDev,
                         PVECT        & NewVect)
{
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() > OldCommDev.size());
    
    //Alloc
    pMapItem mItem;
    
    NewVect.clear();
    NewVect.resize(OldVect.size());
    
    //Main loop
    for(UInt i=1; i <= OldVect.size(); ++i)
    {
      mItem = OldVect.getMapL(i);
      mItem.setPid(NewCommDev.rank());
      
      NewVect.getMapL(i)  = mItem;
      NewVect.getDataL(i) = OldVect.getDataL(i);
    }
  }
}


//_________________________________________________________________________________________________
// OTHER FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DATA>
size_t
pVectGlobalManip<DATA,pMapItem>::
memSize() const
{
  return(sizeof(commDevLoaded));
}



//_________________________________________________________________________________________________
// PMAPITEMSHARE
//-------------------------------------------------------------------------------------------------

/*! Performs global manipulations, retrive informations and checks the \c pVect. Specialized version for \c pMapItemShare */
template<typename DATA> class pVectGlobalManip<DATA,pMapItemShare>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pVect<DATA,pMapItemShare> PVECT;
    typedef sVect<DATA>               DATAVECT;
    typedef pMap<pMapItemShare>       PMAP;
    typedef pMap<pMapItemSendRecv>    SENDRECV;
    typedef pVectUpdateRecursive<DATA,pMapItemShare>     PVUR;
    typedef typename pVectComm<DATA,pMapItemShare>::PVPS PVPS;
    //@}
    
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
    
    
    /*! @name Constructors and set functions */ //@{
  public:
    pVectGlobalManip();
    pVectGlobalManip(const Teuchos::RCP<const communicator> & CommDev);
    pVectGlobalManip(const communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    //@}
    
    
    /*! @name Data retriving functions */ //@{
  public:
    /*! Returns the maximum \c gid on all the processors */
    UInt sizeG(const Teuchos::RCP<const PVECT> & Vect) const;
    
    /*! Returns the maximum \c gid on all the processors */
    UInt sizeG(const PVECT & Vect) const;
    
    /*! Performs a check on the \c pMap contained into \c Vect, see \c pMapGlobalManip */
    bool check(const Teuchos::RCP<const PVECT> & Map) const;
    
    /*! Performs a check on the \c pMap contained into \c Vect, see \c pMapGlobalManip */
    bool check(const PVECT & Map) const;
    
     /*! Checks that items with the same \c gid have the same data */
    bool checkConsistency(const Teuchos::RCP<const PVECT> & Vect) const;
    
    /*! Checks that items with the same \c gid have the same data */
    bool checkConsistency(const PVECT & Vect) const;
    
    /*! Computes the minimum and maximum data */
    void dataMinMax(const Teuchos::RCP<const PVECT> & Vect,
                                               DATA & DataMin,
                                               DATA & DataMax);
    
    /*! Computes the minimum and maximum data */
    void dataMinMax(const PVECT & Vect,
                          DATA & DataMin,
                          DATA & DataMax);
    //@}
    
    
    /*! @name Global manipulations */ //@{
  public:
    /*! Assigns a unique \c gid to data. If two data are equal the same gid is assigned to both.
    The data should implement the inequality and less operator. The ownership structure is also created
    so that elements sharing the same gids are shared and only one of them is regarded as \c owned the others
    are intended as copy of the owned one. The owning structure is only corrected, if correct no fix is performed. */
    void buildGlobalNumbering(Teuchos::RCP<PVECT> & Vect) const;
    
    /*! Assigns a unique \c gid to data. If two data are equal the same gid is assigned to both.
    The data should implement the inequality and less operator. The ownership structure is also created
    so that elements sharing the same gids are shared and only one of them is regarded as \c owned the others
    are intended as copy of the owned one. The owning structure is only corrected, if correct no fix is performed. */
    void buildGlobalNumbering(PVECT & Vect) const;
    
    /*! A higher performance method of \c buildGlobalNumbering. Local elements are marked with 
     the \c isLocal vector and they are not broadcasted. The owning structure is only corrected, if correct no fix is performed.*/
    void buildGlobalNumbering(      Teuchos::RCP<PVECT>              & Vect,
			      const Teuchos::RCP<const sVect<bool> > & isLocal) const;
    
    /*! A higher performance method of \c buildGlobalNumbering. Local elements are marked with 
     the \c isLocal vector and they are not broadcasted. The owning structure is only corrected, if correct no fix is performed.*/
    void buildGlobalNumbering(            PVECT & Vect,
                              const sVect<bool> & isLocal) const;
    
    /*! Assigns a unique \c gid to data and constructs the sharing/owning data structure.
    Only items marked as \c shared = true are assumed to be shared.
    All items that are found to be not shared are marked as not shared.
    If two data are equal the same gid is assigned to both.
    Data should implement the inequality and less operator.*/
    void buildSharingOwnership(Teuchos::RCP<PVECT> & Vect) const;
    
    /*! Assigns a unique \c gid to data and constructs the sharing/owning data structure.
    Only items marked as \c shared = true are assumed to be shared.
    All items that are found to be not shared are marked as not shared.
    If two data are equal the same gid is assigned to both.
    Data should implement the inequality and less operator.*/
    void buildSharingOwnership(PVECT & Vect) const;
    
    /*! Sends and overwrites data of the shared-and-owned to all other shared-and-notOwned items.
    The communication map \c mapSend should be provided.
    See \c pMapGlobalManip to construct the commMap. */
    void updateData(      Teuchos::RCP<PVECT>          & Vect,
                    const Teuchos::RCP<const SENDRECV> & mapSend) const;
    
    /*! Sends and overwrites data of the shared-and-owned to all other shared-and-notOwned items.
    The communication map \c mapSend should be provided.
    See \c pMapGlobalManip to construct the commMap. */
    void updateData(         PVECT & Vect,
                    const SENDRECV & mapSend) const;
    
    /*! Sends and overwrites data of the shared-and-owned to all other shared-and-notOwned items.
    The communication map \c mapSend should be provided.
    See \c pMapGlobalManip to construct the commMap. */
    void updateData(      Teuchos::RCP<PVECT>          & Vect,
                    const Teuchos::RCP<const SENDRECV> & mapSend,
                    const Teuchos::RCP<const SENDRECV> & mapRecv) const;
    
    /*! Sends and overwrites data of the shared-and-owned to all other shared-and-notOwned items.
    The communication map \c mapSend should be provided.
    See \c pMapGlobalManip to construct the commMap. */
    void updateData(         PVECT & Vect,
                    const SENDRECV & mapSend,
                    const SENDRECV & mapRecv) const;
    
    /*! Creates a new vector whose new map is \c NewMap The data contained in \c Vect are communicated 
    to satisfy the new map. An assert chech is performed on the consistency of the new map. */
    void changeMap(      Teuchos::RCP<PVECT>      & Vect,
                   const Teuchos::RCP<const PMAP> & NewMap) const;
    
    /*! Creates a new vector whose new map is \c NewMap The data contained in \c Vect are communicated 
    to satisfy the new map. An assert chech is performed on the consistency of the new map. */
    void changeMap(     PVECT & Vect,
                   const PMAP & NewMap) const;
    
    /*! Distributes the vector using the \c pid logic. The repeated \c gid are eliminated */
    void vectorNormalExclusive(Teuchos::RCP<PVECT> & Vect);
    
    /*! Distributes the vector using the \c pid logic. The repeated \c gid are eliminated */
    void vectorNormalExclusive(PVECT & Vect);
    
    /*! Linear-epetra data distribution 
    \param Vect contains the data to send, then is cleared and loaded with the new data
    */
    Epetra_Map vectorLinear(Teuchos::RCP<PVECT> & Vect);
    
    /*! Linear-epetra data distribution 
    \param Vect contains the data to send, then is cleared and loaded with the new data
    */
    Epetra_Map vectorLinear(PVECT & Vect);
    
    /*! AllReduce - the repeated elements are reduced */
    template<typename OP>
    void allReduce(Teuchos::RCP<PVECT> & Vect, OP op) const;
    
    /*! AllReduce - the repeated elements are reduced */
    template<typename OP>
    void allReduce(PVECT & Vect, OP op) const;
    //@}
    
    /*! @name Non-pending update */ //@{
  public:
    /*! Sends and overwrites data of the shared-and-owned to all other shared-and-notOwned items.
    The communication map \c mapSend should be provided.
    See \c pMapGlobalManip to construct the commMap. */
    void updateDataI(      Teuchos::RCP<PVECT>          & Vect,
                     const Teuchos::RCP<const SENDRECV> & mapSend,
                           sVect<PVECT>                 & commSegments,
                           sVect<PVPS>                  & sendPvps,
                           sVect<PVPS>                  & recvPvps,
                     const UInt                         & channel = 2) const;

    /*! Sends and overwrites data of the shared-and-owned to all other shared-and-notOwned items.
    The communication map \c mapSend should be provided.
    See \c pMapGlobalManip to construct the commMap. */
    void updateDataO(      Teuchos::RCP<PVECT>          & Vect,
                     const Teuchos::RCP<const SENDRECV> & mapSend,
                           sVect<PVECT>                 & commSegments,
                           sVect<PVPS>                  & sendPvps,
                           sVect<PVPS>                  & recvPvps,
                     const UInt                         & channel = 2) const;
    
    /*! Sends and overwrites data of the shared-and-owned to all other shared-and-notOwned items.
    The communication map \c mapSend should be provided.
    See \c pMapGlobalManip to construct the commMap. */
    void updateDataI(      PVECT        & Vect,
                     const SENDRECV     & mapSend,
                           sVect<PVECT> & commSegments,
                           sVect<PVPS>  & sendPvps,
                           sVect<PVPS>  & recvPvps,
                     const UInt         & channel = 2) const;
    
    /*! Sends and overwrites data of the shared-and-owned to all other shared-and-notOwned items.
    The communication map \c mapSend should be provided.
    See \c pMapGlobalManip to construct the commMap. */
    void updateDataO(      PVECT        & Vect,
                     const SENDRECV     & mapSend,
                           sVect<PVECT> & commSegments,
                           sVect<PVPS>  & sendPvps,
                           sVect<PVPS>  & recvPvps,
                     const UInt         & channel = 2) const;
     
    /*! Sends and overwrites data of the shared-and-owned to all other shared-and-notOwned items.
    The communication map \c mapSend should be provided.
    See \c pMapGlobalManip to construct the commMap. */
    void updateDataI(      Teuchos::RCP<PVECT>          & Vect,
                     const Teuchos::RCP<const SENDRECV> & mapSend,
                     const Teuchos::RCP<const SENDRECV> & mapRecv,
                           sVect<PVECT>                 & commSegments,
                           sVect<PVPS>                  & sendPvps,
                           sVect<PVPS>                  & recvPvps,
                     const UInt                         & channel = 2) const;
  
    /*! Sends and overwrites data of the shared-and-owned to all other shared-and-notOwned items.
    The communication map \c mapSend should be provided.
    See \c pMapGlobalManip to construct the commMap. */
    void updateDataO(      Teuchos::RCP<PVECT>          & Vect,
                     const Teuchos::RCP<const SENDRECV> & mapSend,
                     const Teuchos::RCP<const SENDRECV> & mapRecv,
                           sVect<PVECT>                 & commSegments,
                           sVect<PVPS>                  & sendPvps,
                           sVect<PVPS>                  & recvPvps,
                     const UInt                         & channel = 2) const;
    
    /*! Sends and overwrites data of the shared-and-owned to all other shared-and-notOwned items.
    The communication map \c mapSend should be provided.
    See \c pMapGlobalManip to construct the commMap. */
    void updateDataI(      PVECT        & Vect,
                     const SENDRECV     & mapSend,
                     const SENDRECV     & mapRecv,
                           sVect<PVECT> & commSegments,
                           sVect<PVPS>  & sendPvps,
                           sVect<PVPS>  & recvPvps,
                     const UInt         & channel = 2) const;
   
    /*! Sends and overwrites data of the shared-and-owned to all other shared-and-notOwned items.
    The communication map \c mapSend should be provided.
    See \c pMapGlobalManip to construct the commMap. */
    void updateDataO(      PVECT        & Vect,
                     const SENDRECV     & mapSend,
                     const SENDRECV     & mapRecv,
                           sVect<PVECT> & commSegments,
                           sVect<PVPS>  & sendPvps,
                           sVect<PVPS>  & recvPvps,
                     const UInt         & channel = 2) const;  
    //@}
     

    /*! @name Recursive update */ //@{
  public:
    void updateDataRR(      Teuchos::RCP<PVECT>          & Vect,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
     
    void updateDataRR(      PVECT    & Vect,
                      const SENDRECV & mapSend,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
      
    void updateDataRI(      Teuchos::RCP<PVECT>          & Vect,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
     
    void updateDataRI(      PVECT    & Vect,
                      const SENDRECV & mapSend,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
      
    void updateDataRO(      Teuchos::RCP<PVECT>          & Vect,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
     
    void updateDataRO(      PVECT    & Vect,
                      const SENDRECV & mapSend,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
      
    void updateDataRR(      Teuchos::RCP<PVECT>          & Vect,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                      const Teuchos::RCP<const SENDRECV> & mapRecv,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
     
    void updateDataRR(      PVECT    & Vect,
                      const SENDRECV & mapSend,
                      const SENDRECV & mapRecv,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
      
    void updateDataRI(      Teuchos::RCP<PVECT>          & Vect,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                      const Teuchos::RCP<const SENDRECV> & mapRecv,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
     
    void updateDataRI(      PVECT    & Vect,
                      const SENDRECV & mapSend,
                      const SENDRECV & mapRecv,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
      
    void updateDataRO(      Teuchos::RCP<PVECT>          & Vect,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                      const Teuchos::RCP<const SENDRECV> & mapRecv,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
     
    void updateDataRO(      PVECT    & Vect,
                      const SENDRECV & mapSend,
                      const SENDRECV & mapRecv,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
    //@}
  
      
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copies the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const PVECT>        & OldVect,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<PVECT>              & NewVect);
    
    /*! Copies the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const PVECT        & OldVect,
                            const communicator & NewCommDev,
                                  PVECT        & NewVect);
    
    /*! Copies the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const PVECT>        & OldVect,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<PVECT>              & NewVect);
    
    /*! Copies the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const PVECT        & OldVect,
                            const communicator & NewCommDev,
                                  PVECT        & NewVect);
    //@}
    
    /*! @name Other functions */ //@{
  public:
    size_t memSize() const;
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DATA>
pVectGlobalManip<DATA,pMapItemShare>::
pVectGlobalManip() : commDevLoaded(false)
{
}

template<typename DATA>
pVectGlobalManip<DATA,pMapItemShare>::
pVectGlobalManip(const Teuchos::RCP<const communicator> & CommDev) : commDevLoaded(true), commDev(CommDev)
{
}

template<typename DATA>
pVectGlobalManip<DATA,pMapItemShare>::
pVectGlobalManip(const communicator & CommDev) : commDevLoaded(true)
{
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}


//_________________________________________________________________________________________________
// DATA RETRIVING FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DATA>
UInt
pVectGlobalManip<DATA,pMapItemShare>::
sizeG(const Teuchos::RCP<const PVECT> & Vect) const
{
  assert(commDevLoaded);
  
  pMapGlobalManip<pMapItemShare> manip(commDev);
  return(manip.sizeG(Vect->getMapRef()));
}

template<typename DATA>
UInt
pVectGlobalManip<DATA,pMapItemShare>::
sizeG(const PVECT & Vect) const
{
  assert(commDevLoaded);
  
  pMapGlobalManip<pMapItemShare> manip(commDev);
  return(manip.sizeG(Vect.getMapRef()));
}

template<typename DATA>
bool
pVectGlobalManip<DATA,pMapItemShare>::
check(const Teuchos::RCP<const PVECT> & Vect) const
{
  assert(commDevLoaded);
  
  pMapGlobalManip<pMapItemShare> manip(commDev);
  return(manip.check(Vect->getMapRef()));
}

template<typename DATA>
bool
pVectGlobalManip<DATA,pMapItemShare>::
check(const PVECT & Vect) const
{
  assert(commDevLoaded);
  
  pMapGlobalManip<pMapItemShare> manip(commDev);
  return(manip.check(Vect.getMapRef()));
}

template<typename DATA>
bool
pVectGlobalManip<DATA,pMapItemShare>::
checkConsistency(const Teuchos::RCP<const PVECT> & Vect) const
{
  //Assert
  assert(commDevLoaded);
  return(checkConsistency(*Vect));
}

template<typename DATA>
bool
pVectGlobalManip<DATA,pMapItemShare>::
checkConsistency(const PVECT & Vect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Vect copy
  PVECT vectCopy = Vect;
  
  //Data redistribution
  pVectComm<DATA,pMapItemShare> distributor(commDev);
  distributor.vectorNormal(vectCopy);
  
  //Ordering
  pVectManip<DATA,pMapItemShare> orderer;
  orderer.orderGid(vectCopy);
  
  //Checking
  bool flag = true;
  
  if(vectCopy.size() > 1)
  {
    for(UInt i=1; i < vectCopy.size(); ++i)
    {
      if(!(vectCopy.getMapL(i).getGid() != vectCopy.getMapL(i+1).getGid()))
      {
	flag = flag & (!(vectCopy.getDataL(i) != vectCopy.getDataL(i+1)));
      }
    }
  }
  
  //Global value reconstruction  
  bool flagGlob = flag;
  sVect<UInt> listTemp;
  
  all_gather(*commDev,UInt(flag),listTemp);
  
  for(UInt i=1; i <= listTemp.size(); ++i)
  { flagGlob = flagGlob & bool(listTemp(i)); }
  
  return(flagGlob);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
dataMinMax(const Teuchos::RCP<const PVECT> & Vect,
                                      DATA & DataMin,
                                      DATA & DataMax)
{
  assert(commDevLoaded);
  dataMinMax(*Vect,DataMin,DataMax);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
dataMinMax(const PVECT & Vect,
                  DATA & DataMin,
                  DATA & DataMax)
{
  assert(commDevLoaded);
  
  //Allocations
  UInt numProc = commDev->size();
  
  //Local maximum and minimum data determination
  DATA dataMax, dataMin;
  dataMax = dataMin;
  
  UInt dataLoaded = (Vect.size() >= 1);
  
  if(dataLoaded)
  {
    dataMin = Vect.getL(1);
    dataMax = Vect.getL(1);
  }
  
  for(UInt i=1; i <= Vect.size(); ++i)
  {
    dataMin = std::min(dataMin,Vect.getL(i));
    dataMax = std::max(dataMax,Vect.getL(i));
  }  
  
  //Global maximum and minimum data determination
  sVect<UInt> allLoading(numProc);
  morgana_all_gather(*commDev,dataLoaded,allLoading);

  sVect<DATA> allMin(numProc);
  morgana_all_gather(*commDev,dataMin,allMin);
  
  sVect<DATA> allMax(numProc);
  morgana_all_gather(*commDev,dataMax,allMax);
  
  for(UInt i=1; i<=numProc; ++i)
  {
    if(allLoading(i))
    {
      dataMin = min(dataMin,allMin(i));
      dataMax = max(dataMax,allMax(i));
    }
  }
  
  //Data assignment
  DataMin = dataMin;
  DataMax = dataMax;
}


//_________________________________________________________________________________________________
// GLOBAL MANIPULATION FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
buildGlobalNumbering(Teuchos::RCP<PVECT> & Vect) const
{
  assert(commDevLoaded);
  buildGlobalNumbering(*Vect);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
buildGlobalNumbering(PVECT & Vect) const
{
  assert(commDevLoaded);
  
  //Allocations
  UInt numProc = commDev->size();
  UInt proc    = commDev->rank();
  
  //Vect copy
  Vect.bufferLids();
  PVECT vectCopy = Vect;
  
  //Pid buffering
  for(UInt i=1; i<=vectCopy.size(); ++i)
  {
    vectCopy.getMapL(i).setPid(proc);
  }
  
  //Data redistribution
  pVectComm<DATA,pMapItemShare> distributor(commDev);
  distributor.vectorData(vectCopy);
  
  //Ordering
  pVectManip<DATA,pMapItemShare> orderer;
  orderer.orderData(vectCopy);
  
  //Counting the shared and owned items
  typedef pair<DATA,UInt> PAIR;
  typedef typename map<DATA,UInt>::iterator ITER;
  
  PAIR pair;
  ITER iter;
  std::pair<ITER,bool> ret;
  map<DATA,UInt> numShared, numOwned; 
  
  for(UInt i=1; i <= vectCopy.size(); ++i)
  {
    pair.first  = vectCopy(i);
    
    //Sharing count
    pair.second = 1;
    ret = numShared.insert(pair);
    
    if(!ret.second)
    {
      iter = ret.first;
      iter->second = iter->second + 1;
    }
    
    //Owning count
    pair.second = UInt(vectCopy.getMapL(i).getOwned());
    ret = numOwned.insert(pair);
    
    if(!ret.second)
    {
      iter = ret.first;
      iter->second = iter->second + UInt(vectCopy.getMapL(i).getOwned());
    }
  }
  
  //Local numbering
  UInt k = 1;
  bool owned  = true;
  bool shared = false;
  
  if(vectCopy.size() == 0) { k = 0; }
  
  for(UInt i=1; i <= vectCopy.size(); ++i)
  {
    //Sharing determination
    iter   = numShared.find(vectCopy(i));
    shared = (iter->second >= 2);    
    
    //Owning determination
    iter = numOwned.find(vectCopy(i));
    
    if(iter->second == 1)
    {
      owned = vectCopy.getMapL(i).getOwned();
    }
    else
    {
      if(iter->second == 0)
      {
	owned = true;
	iter->second = 1;
      }
      else
      {
	if(vectCopy.getMapL(i).getOwned()) 
	{
	  owned = false;
	  iter->second -= 1;
	}
	else
	{
	  owned = vectCopy.getMapL(i).getOwned();
	}
      }
    }
    
    //Settings
    vectCopy.getMapL(i).setGid(k);
    vectCopy.getMapL(i).setOwned(owned);
    vectCopy.getMapL(i).setShared(shared);
    
    //Se esite un elemento successivo
    if(i != vectCopy.size())
    {
      if(vectCopy(i) != vectCopy(i+1))
      { k++; }
    }
  }  
  
  //Offset determination
  sVect<UInt> localSizes(numProc);
  all_gather(*commDev,k,localSizes);
  
  UInt offset = 0;
  for(UInt i=0; i < proc; ++i)
  {
    offset += localSizes[i];
  }  
  
  //Global numbering
  for(UInt i=1; i <= vectCopy.size(); ++i)
  {
    k = vectCopy.getMapL(i).getGid();
    vectCopy.getMapL(i).setGid(k + offset);
  }
  
  //Pid map communication
  distributor.vectorPid(vectCopy);
  
  //Data download
  Vect = vectCopy;
  Vect.restoreLids();
  
  //Data lid ordering
  pVectManip<DATA,pMapItemShare> manipulator;
  manipulator.setIndexing(Vect);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
buildGlobalNumbering(      Teuchos::RCP<PVECT>              & Vect,
                     const Teuchos::RCP<const sVect<bool> > & isLocal) const
{
  //Assert
  assert(commDevLoaded);
  assert(Vect->size() == isLocal->size());
  
  buildGlobalNumbering(*Vect,*isLocal);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
buildGlobalNumbering(            PVECT & Vect,
                     const sVect<bool> & isLocal) const
{ 
  //Assert
  assert(commDevLoaded);
  assert(Vect.size() == isLocal.size());
  
  //Allocations
  UInt numProc = commDev->size();
  UInt proc    = commDev->rank();
  
  //Vector splitting
  sVect<UInt> lLids, gLids;
  PVECT lVect, gVect;
  
  for(UInt i=1; i <= Vect.size(); ++i)
  {
    if(isLocal(i))
    {
      Vect.getMapL(i).setShared(false);
      Vect.getMapL(i).setOwned(true);
      
      lVect.push_back(Vect.getMapL(i), Vect.getDataL(i));
      lLids.push_back(i);
    }
    else
    {
      gVect.push_back(Vect.getMapL(i), Vect.getDataL(i));
      gLids.push_back(i);
    }
  }
  
  lVect.updateFinder();
  gVect.updateFinder();  
  
  //Broadcast offsets
  sVect<UInt> localSizes(numProc);
  all_gather(*commDev, lVect.size(), localSizes);
  
  UInt offset    = 0;
  UInt totOffset = 0;
  
  for(UInt i=0; i < proc; ++i)
  { offset += localSizes[i]; }
  
  for(UInt i=0; i < numProc; ++i)
  { totOffset += localSizes[i]; }
  
  //Global vector temp numbering
  buildGlobalNumbering(gVect);
  
  //Local vector offsetting 
  for(UInt i=1; i <= lVect.size(); ++i)
  {
    lVect.getMapL(i).setLid(lLids(i));
    lVect.getMapL(i).setGid(i + offset);
  }   
   
  //Global vector offsetting
  UInt gid;
  
  for(UInt i=1; i <= gVect.size(); i++)
  {
    gid = gVect.getMapL(i).getGid();
    
    gVect.getMapL(i).setLid(gLids(i));
    gVect.getMapL(i).setGid(gid + totOffset);
    
    lVect.push_back(gVect.getMapL(i),gVect.getDataL(i));
  }
  
  //Lid reordering
  pVectManip<DATA,pMapItemShare> localManip;
  localManip.setIndexing(lVect);
  
  //Copy
  Vect = lVect;
  Vect.updateFinder();
  assert(check(Vect));
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
buildSharingOwnership(Teuchos::RCP<PVECT> & Vect) const
{
  assert(commDevLoaded);
  buildSharingOwnership(*Vect);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
buildSharingOwnership(PVECT & Vect) const
{
  assert(commDevLoaded);
  
  //Allocations
  UInt numProc = commDev->size();
  UInt proc    = commDev->rank();
  
  //Vect copy
  Vect.bufferLids();
  
  //Pid buffering
  for(UInt i=1; i<=Vect.size(); ++i)
  { Vect.getMapL(i).setPid(proc); }
  
  //Vect subdivision
  PVECT vectShared, vectNotShared;
  
  for(UInt i=1; i<=Vect.size(); ++i)
  {
    if(Vect.getMapL(i).getShared())
    { vectShared.push_back(Vect.getMapL(i),Vect.get(i)); }
    else
    { vectNotShared.push_back(Vect.getMapL(i),Vect.get(i)); }
  }
  
  vectShared.updateFinder();
  vectNotShared.updateFinder();
  
  //Data redistribution
  pVectComm<DATA,pMapItemShare> distributor(commDev);
  distributor.vectorData(vectShared);
  
  //Vector merging
  PVECT vectCopy;
  
  sVect<PVECT> sourceVects;
  sourceVects.push_back(vectShared);
  sourceVects.push_back(vectNotShared);
  
  pVectManip<DATA,pMapItemShare> manipulator;
  manipulator.mergeSimple(sourceVects,vectCopy);
  
  //Local numbering
  UInt k = 1;
  bool owned  = true;
  bool shared = false;
  
  if(vectCopy.size() == 0) { k = 0; }
  
  for(UInt i=1; i <= vectCopy.size(); ++i)
  {
    //Sharing determination
    shared = false;
    
    if(i != vectCopy.size())
    {
      if( !(vectCopy(i) != vectCopy(i+1)) )
      { shared = true; }
    }
    
    if(i != 1)
    {
      if( !(vectCopy(i) != vectCopy(i-1)) )
      { shared = true; }
    }
    
    //Settings
    vectCopy.getMapL(i).setGid(k);
    vectCopy.getMapL(i).setOwned(owned);
    vectCopy.getMapL(i).setShared(shared);
    
    //Se esite un elemento successivo
    if(i != vectCopy.size())
    {
      if(vectCopy(i) != vectCopy(i+1))
      { 
	//Incremento del gid
	k++;
	
	//Set owning flag
	owned  = true;
      }
      else
      {
	//Set owning flag
	owned  = false;
      }
    }
  }
  
  //Offset determination
  sVect<UInt> localSizes(numProc);
  all_gather(*commDev,k,localSizes);
  
  UInt offset = 0;
  for(UInt i=0; i < proc; ++i)
  {
    offset += localSizes[i];
  }
  
  //Global numbering
  for(UInt i=1; i <= vectCopy.size(); ++i)
  {
    k = vectCopy.getMapL(i).getGid();
    vectCopy.getMapL(i).setGid(k + offset);
  }
  
  //Pid map communication
  distributor.vectorPid(vectCopy);
  
  //Data download
  Vect = vectCopy;
  Vect.restoreLids();
  
  //Data lid ordering
  manipulator.setIndexing(Vect);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateData(      Teuchos::RCP<PVECT>          & Vect,
           const Teuchos::RCP<const SENDRECV> & mapSend) const
{
  assert(commDevLoaded);
  updateData(*Vect,*mapSend);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateData(         PVECT & Vect,
           const SENDRECV & mapSend) const
{
  //Assert
  assert(commDevLoaded);
  
  //Receiving vect
  PVECT rVect;
  
  //Sending
  pVectComm<DATA,pMapItemShare> sender(commDev);
  sender.sendRecv(mapSend,Vect,rVect);
  
  //Data update
  UInt gid;
  
  for(UInt i=1; i<=rVect.size(); ++i)
  {
    gid = rVect.getMapL(i).getGid();
    Vect.getG(gid) = rVect.get(i);
  }
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateData(      Teuchos::RCP<PVECT>          & Vect,
           const Teuchos::RCP<const SENDRECV> & mapSend,
           const Teuchos::RCP<const SENDRECV> & mapRecv) const
{
  assert(commDevLoaded);
  updateData(*Vect,*mapSend,*mapRecv);
}    

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateData(         PVECT & Vect,
           const SENDRECV & mapSend,
           const SENDRECV & mapRecv) const
{
  //Assert
  assert(commDevLoaded);
  
  //Receiving vect
  PVECT rVect;
  
  //Sending
  pVectComm<DATA,pMapItemShare> sender(commDev);
  sender.sendRecv(mapSend,mapRecv,Vect,rVect);
  
  //Data update
  UInt gid;
  
  for(UInt i=1; i<=rVect.size(); ++i)
  {
    gid = rVect.getMapL(i).getGid();
    Vect.getG(gid) = rVect.get(i);
  }
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
changeMap(      Teuchos::RCP<PVECT>      & Vect,
          const Teuchos::RCP<const PMAP> & NewMap) const
{
  assert(commDevLoaded);
  changeMap(*Vect,*NewMap);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
changeMap(     PVECT & Vect,
          const PMAP & NewMap) const
{
  assert(commDevLoaded);
  
  //Pre-checking
  pMapManip<pMapItemShare>       preLocal;
  pMapGlobalManip<pMapItemShare> preGlobal(commDev);
  
  assert(preLocal.checkNormalIndexing(NewMap));
  assert(preGlobal.check(NewMap));
  
  //Processor ids
  UInt proc = commDev->rank();
  
  //Allocations and vector copy
  Vect.bufferLids();
  PVECT vectCopy = Vect;
  PVECT newVect;
  
  //Maximum gid determination
  UInt maxGid = sizeG(vectCopy);
  
  //Pid buffering
  for(UInt i=1; i<=vectCopy.size(); ++i)
  { vectCopy.getMapL(i).setPid(proc); }
  
  //Local items of the vector
  UInt gid;
  PMAP retrivingMap;
  
  for(UInt i=1; i <= NewMap.size(); ++i)
  {
    gid = NewMap.get(i).getGid();
    
    if(Vect.isG(gid))
    { newVect.push_back(Vect.getDataG(gid),NewMap.get(i),false); }
    else
    { retrivingMap.push_back(NewMap.get(i)); }
  }
  
  retrivingMap.bufferLids();
  
  //Data normal indexing 
  pVectComm<DATA,pMapItemShare> vectDistributor(commDev);
  vectDistributor.vectorNormal(vectCopy);
  
  pMapComm<pMapItemShare> mapDistributor(commDev);
  mapDistributor.vectorNormal(retrivingMap,maxGid);
  
  retrivingMap.restoreLids();
  vectCopy.restoreLids();
  
  //Build commPattern
  SENDRECV sendMap(retrivingMap.size());
  pMapItemSendRecv sendItem;
  
  for(UInt i=1; i <= retrivingMap.size(); ++i)
  {
    gid = retrivingMap(i).getGid();
    assert(vectCopy.isG(gid));
    
    sendItem.setPid( vectCopy.getMapG(gid).getPid() );
    sendItem.setLid( vectCopy.getMapG(gid).getLid() );
    sendItem.setGid( gid );
    sendItem.setSid( vectCopy.getMapG(gid).getPid() );
    sendItem.setRid( retrivingMap(i).getPid() );
    
    sendMap(i) = sendItem;
  }
  
  //Pid indexing of send-recv
  pMapComm<pMapItemSendRecv> sendDistributor(commDev);
  sendDistributor.vectorPid(sendMap);
  
  //Sending of data-items
  PVECT retrivedVect;
  vectDistributor.sendRecv(sendMap,Vect,retrivedVect);
  
  //VectMerging
  for(UInt i=1; i <= NewMap.size(); ++i)
  {
    gid = NewMap.get(i).getGid();
    
    if(retrivedVect.isG(gid))
    { newVect.push_back(retrivedVect.getDataG(gid),NewMap.get(i),false); }
  }
  
  //Data reordering  
  pVectManip<DATA,pMapItemShare> vectManipulator;
  vectManipulator.setIndexing(newVect);
  
  newVect.updateFinder();
  Vect = newVect;
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
vectorNormalExclusive(Teuchos::RCP<PVECT> & Vect)
{
  assert(commDevLoaded);
  vectorNormalExclusive(*Vect);
}
    
template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
vectorNormalExclusive(PVECT & Vect)
{
  assert(commDevLoaded);
  
  //Push to pid distribution
  pVectComm<DATA,pMapItemShare> distributor(commDev);
  distributor.vectorNormal(Vect);
  
  //Eliminate repeated gids
  pVectManip<DATA,pMapItemShare> serialManip;
  serialManip.unionExclusive(Vect);
  
  //Set pid
  UInt pid = commDev->rank();
  
  for(UInt i=1; i <= Vect.size(); ++i)
  {
    mapItemFixer(Vect.getMapL(i));
    Vect.getMapL(i).setPid(pid);  
  }

  //Checking
  assert(check(Vect));
}

template<typename DATA>
Epetra_Map
pVectGlobalManip<DATA,pMapItemShare>::
vectorLinear(Teuchos::RCP<PVECT> & Vect)
{
  assert(commDevLoaded);
  
  return(vectorLinear(*Vect));
}

template<typename DATA>
Epetra_Map
pVectGlobalManip<DATA,pMapItemShare>::
vectorLinear(PVECT & Vect)
{
  //Assert and pid
  assert(commDevLoaded);
  UInt pid = commDev->rank();
  
  //New map epetra
  UInt mapSize = sizeG(Vect);
  Epetra_MpiComm epetraComm(*commDev);
  Epetra_Map     newMap_epetra(int(mapSize), int(0), epetraComm);
  
  //New map morgana
  pMapGlobalManip<pMapItemShare> mapImporter(commDev);
  Teuchos::RCP<PMAP> newMap_morgana(new PMAP);
  mapImporter.importEpetraMap(newMap_epetra,newMap_morgana);
  
  for(UInt i=1; i <= newMap_morgana->size(); ++i)
  {
    mapItemFixer(newMap_morgana->get(i));
    newMap_morgana->get(i).setPid(pid);
  }
  
  //Redistribute data with the new map
  changeMap(Vect,*newMap_morgana);
  
  //Return
  return(newMap_epetra);
}

template<typename DATA>
template<typename OP>
void
pVectGlobalManip<DATA,pMapItemShare>::
allReduce(Teuchos::RCP<PVECT> & Vect, OP op) const
{
  allReduce(*Vect,op);
}
    
template<typename DATA>
template<typename OP>
void
pVectGlobalManip<DATA,pMapItemShare>::
allReduce(PVECT & Vect, OP op) const
{
  //Assert-----------------------------------------------------------
  assert(commDevLoaded);
  assert(check(Vect));
  
  //Typedef----------------------------------------------------------
  typedef traitsBasic<DATA>          TRAITS;
  typedef std::map<UInt,DATA>        STD_MAP;
  typedef std::pair<UInt,DATA>       STD_PAIR;
  typedef typename STD_MAP::iterator STD_ITER;
  typedef std::pair<STD_ITER,bool>   STD_RET;
  
  //Create comm sub-vect---------------------------------------------
  PVECT commVect;
  
  for(UInt i=1; i <= Vect.size(); ++i)
  {
    if(Vect.getMapL(i).getShared())
    { 
      commVect.push_back(Vect.getMapL(i),
                         Vect.getDataL(i));
    }
  }
  
  commVect.updateFinder();
  
  //Normal distribution----------------------------------------------
  pVectComm<DATA,pMapItemShare> pComm(commDev);
  pComm.vectorNormal(commVect);
  
  //Perform the reduction--------------------------------------------
  STD_MAP  stdMap;
  STD_RET  stdRet;
  STD_PAIR stdPair;
  STD_ITER stdIter;
  Real valueA, valueB, valueT;
  UInt gid;
  
  static const UInt numI = traitsBasic<DATA>::numI;
  static const UInt numJ = traitsBasic<DATA>::numJ;
  
  for(UInt i=1; i <= commVect.size(); ++i)
  {
    stdPair.first  = commVect.getMapL(i).getGid();
    stdPair.second = commVect.getDataL(i);
    
    stdRet = stdMap.insert(stdPair);
    
    if(!stdRet.second)
    {
      stdIter = stdRet.first;
      
      for(UInt I=1; I <= numI; ++I)
      {
        for(UInt J=1; J <= numJ; J++)
        {
          valueA = TRAITS::getIJ(I,J,stdIter->second);
          valueB = TRAITS::getIJ(I,J,commVect(i));
          valueT = op.operator()(valueA,valueB);
          TRAITS::setIJ(I,J,valueT,stdIter->second);
        }
      }
    } 
  }
  
  for(UInt i=1; i <= commVect.size(); ++i)
  {
    gid = commVect.getMapL(i).getGid();
    assert(stdMap.count(gid));
    
    stdIter     = stdMap.find(gid);
    commVect(i) = stdIter->second;
  }
  
  //Comm back--------------------------------------------------------
  pComm.vectorPid(commVect);
  
  //Insert in the original vector------------------------------------
  for(UInt i=1; i <= commVect.size(); ++i)
  {
    gid = commVect.getMapL(i).getGid();
    
    assert(Vect.isG(gid));
    Vect.getDataG(gid) = commVect(i);
  }
}


//_________________________________________________________________________________________________
// NON PENDING UPDATE
//-------------------------------------------------------------------------------------------------
template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataI(      Teuchos::RCP<PVECT>          & Vect,
            const Teuchos::RCP<const SENDRECV> & mapSend,
                  sVect<PVECT>                 & commSegments,
                  sVect<PVPS>                  & sendPvps,
                  sVect<PVPS>                  & recvPvps,
            const UInt                         & channel) const
{
  assert(commDevLoaded);
  updateDataI(*Vect,
              *mapSend,
               commSegments,
               sendPvps,
               recvPvps,
               channel);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataO(      Teuchos::RCP<PVECT>          & Vect,
            const Teuchos::RCP<const SENDRECV> & mapSend,
                  sVect<PVECT>                 & commSegments,
                  sVect<PVPS>                  & sendPvps,
                  sVect<PVPS>                  & recvPvps,
            const UInt                         & channel) const
{
  assert(commDevLoaded);
  updateDataO(*Vect,
              *mapSend,
               commSegments,
               sendPvps,
               recvPvps,
               channel);
}
   
template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataI(      PVECT        & Vect,
            const SENDRECV     & mapSend,
                  sVect<PVECT> & commSegments,
                  sVect<PVPS>  & sendPvps,
                  sVect<PVPS>  & recvPvps,
            const UInt         & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sending
  PVECT bufVect;
  
  pVectComm<DATA,pMapItemShare> sender(commDev);
  sender.sendRecvI(mapSend,
                   Vect,
                   bufVect,
                   commSegments,
                   sendPvps,
                   recvPvps,
                   channel);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataO(      PVECT        & Vect,
            const SENDRECV     & mapSend,
                  sVect<PVECT> & commSegments,
                  sVect<PVPS>  & sendPvps,
                  sVect<PVPS>  & recvPvps,
            const UInt         & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sending
  PVECT bufVect;
  
  pVectComm<DATA,pMapItemShare> sender(commDev);
  sender.sendRecvO(mapSend,
                   Vect,
                   bufVect,
                   commSegments,
                   sendPvps,
                   recvPvps,
                   channel);
  
  //Data update
  UInt gid;
  
  for(UInt i=1; i <= bufVect.size(); ++i)
  {
    gid = bufVect.getMapL(i).getGid();
    Vect.getG(gid) = bufVect.get(i);
  }
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataI(      Teuchos::RCP<PVECT>          & Vect,
            const Teuchos::RCP<const SENDRECV> & mapSend,
            const Teuchos::RCP<const SENDRECV> & mapRecv,
                  sVect<PVECT>                 & commSegments,
                  sVect<PVPS>                  & sendPvps,
                  sVect<PVPS>                  & recvPvps,
            const UInt                         & channel) const
{
  assert(commDevLoaded);
  updateDataI(*Vect,
              *mapSend,
              *mapRecv,
               commSegments,
               sendPvps,
               recvPvps,
               channel);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataO(      Teuchos::RCP<PVECT>          & Vect,
            const Teuchos::RCP<const SENDRECV> & mapSend,
            const Teuchos::RCP<const SENDRECV> & mapRecv,
                  sVect<PVECT>                 & commSegments,
                  sVect<PVPS>                  & sendPvps,
                  sVect<PVPS>                  & recvPvps,
            const UInt                         & channel) const
{
  assert(commDevLoaded);
  updateDataO(*Vect,
              *mapSend,
              *mapRecv,
               commSegments,
               sendPvps,
               recvPvps,
               channel);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataI(      PVECT        & Vect,
            const SENDRECV     & mapSend,
            const SENDRECV     & mapRecv,
                  sVect<PVECT> & commSegments,
                  sVect<PVPS>  & sendPvps,
                  sVect<PVPS>  & recvPvps,
            const UInt         & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sending
  PVECT bufVect;
  
  pVectComm<DATA,pMapItemShare> sender(commDev);
  sender.sendRecvI(mapSend,
                   mapRecv,
                   Vect,
                   bufVect,
                   commSegments,
                   sendPvps,
                   recvPvps,
                   channel);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataO(      PVECT        & Vect,
            const SENDRECV     & mapSend,
            const SENDRECV     & mapRecv,
                  sVect<PVECT> & commSegments,
                  sVect<PVPS>  & sendPvps,
                  sVect<PVPS>  & recvPvps,
            const UInt         & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sending
  PVECT bufVect;
  
  pVectComm<DATA,pMapItemShare> sender(commDev);
  sender.sendRecvO(mapSend,
                   mapRecv,
                   Vect,
                   bufVect,
                   commSegments,
                   sendPvps,
                   recvPvps,
                   channel);
  
  //Data update
  UInt gid;
  
  for(UInt i=1; i <= bufVect.size(); ++i)
  {
    gid = bufVect.getMapL(i).getGid();
    Vect.getG(gid) = bufVect.get(i);
  }
}



//_________________________________________________________________________________________________
// RECURSIVE UPDATE
//-------------------------------------------------------------------------------------------------
template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataRR(      Teuchos::RCP<PVECT>          & Vect,
             const Teuchos::RCP<const SENDRECV> & mapSend,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRR(*Vect,
               *mapSend,
                commBuffer,
                channel);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataRR(      PVECT    & Vect,
             const SENDRECV & mapSend,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sending
  pVectComm<DATA,pMapItemShare> sender(commDev);
  sender.sendRecvRR(mapSend,
                    Vect,
                    commBuffer.rVect,
                    commBuffer.sendSegments,
                    commBuffer.recvSegments,
                    commBuffer.sendPvrs,
                    commBuffer.recvPvrs,
                    channel);
  
  //Data update
  UInt gid;
  
  for(UInt i=1; i <= commBuffer.rVect.size(); ++i)
  {
    gid = commBuffer.rVect.getMapL(i).getGid();
    Vect.getG(gid) = commBuffer.rVect.get(i);
  }
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataRI(      Teuchos::RCP<PVECT>          & Vect,
             const Teuchos::RCP<const SENDRECV> & mapSend,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRI(*Vect,
               *mapSend,
                commBuffer,
                channel);
}
     
template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataRI(      PVECT    & Vect,
             const SENDRECV & mapSend,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sending
  pVectComm<DATA,pMapItemShare> sender(commDev);
  sender.sendRecvRI(mapSend,
                    Vect,
                    commBuffer.rVect,
                    commBuffer.sendSegments,
                    commBuffer.recvSegments,
                    commBuffer.sendPvrs,
                    commBuffer.recvPvrs,
                    channel);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataRO(      Teuchos::RCP<PVECT>          & Vect,
             const Teuchos::RCP<const SENDRECV> & mapSend,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRO(*Vect,
               *mapSend,
                commBuffer,
                channel);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataRO(      PVECT    & Vect,
             const SENDRECV & mapSend,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sending
  pVectComm<DATA,pMapItemShare> sender(commDev);
  sender.sendRecvRO(mapSend,
                    Vect,
                    commBuffer.rVect,
                    commBuffer.sendSegments,
                    commBuffer.recvSegments,
                    commBuffer.sendPvrs,
                    commBuffer.recvPvrs,
                    channel);
  
  //Data update
  UInt gid;
  
  for(UInt i=1; i <= commBuffer.rVect.size(); ++i)
  {
    gid = commBuffer.rVect.getMapL(i).getGid();
    Vect.getG(gid) = commBuffer.rVect.get(i);
  }
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataRR(      Teuchos::RCP<PVECT>          & Vect,
             const Teuchos::RCP<const SENDRECV> & mapSend,
             const Teuchos::RCP<const SENDRECV> & mapRecv,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRR(*Vect,
               *mapSend,
               *mapRecv,
                commBuffer,
                channel);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataRR(      PVECT    & Vect,
             const SENDRECV & mapSend,
             const SENDRECV & mapRecv,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sending
  pVectComm<DATA,pMapItemShare> sender(commDev);
  sender.sendRecvRR(mapSend,
                    mapRecv,
                    Vect,
                    commBuffer.rVect,
                    commBuffer.sendSegments,
                    commBuffer.recvSegments,
                    commBuffer.sendPvrs,
                    commBuffer.recvPvrs,
                    channel);
  
  //Data update
  UInt gid;
  
  for(UInt i=1; i <= commBuffer.rVect.size(); ++i)
  {
    gid = commBuffer.rVect.getMapL(i).getGid();
    Vect.getG(gid) = commBuffer.rVect.get(i);
  }
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataRI(      Teuchos::RCP<PVECT>          & Vect,
             const Teuchos::RCP<const SENDRECV> & mapSend,
             const Teuchos::RCP<const SENDRECV> & mapRecv,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRI(*Vect,
               *mapSend,
               *mapRecv,
                commBuffer,
                channel);
}
     
template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataRI(      PVECT    & Vect,
             const SENDRECV & mapSend,
             const SENDRECV & mapRecv,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sending
  pVectComm<DATA,pMapItemShare> sender(commDev);
  sender.sendRecvRI(mapSend,
                    mapRecv,
                    Vect,
                    commBuffer.rVect,
                    commBuffer.sendSegments,
                    commBuffer.recvSegments,
                    commBuffer.sendPvrs,
                    commBuffer.recvPvrs,
                    channel);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataRO(      Teuchos::RCP<PVECT>          & Vect,
             const Teuchos::RCP<const SENDRECV> & mapSend,
             const Teuchos::RCP<const SENDRECV> & mapRecv,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRO(*Vect,
               *mapSend,
               *mapRecv,
                commBuffer,
                channel);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
updateDataRO(      PVECT    & Vect,
             const SENDRECV & mapSend,
             const SENDRECV & mapRecv,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sending
  pVectComm<DATA,pMapItemShare> sender(commDev);
  sender.sendRecvRO(mapSend,
                    mapRecv,
                    Vect,
                    commBuffer.rVect,
                    commBuffer.sendSegments,
                    commBuffer.recvSegments,
                    commBuffer.sendPvrs,
                    commBuffer.recvPvrs,
                    channel);
  
  //Data update
  UInt gid;
  
  for(UInt i=1; i <= commBuffer.rVect.size(); ++i)
  {
    gid = commBuffer.rVect.getMapL(i).getGid();
    Vect.getG(gid) = commBuffer.rVect.get(i);
  }
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const PVECT>        & OldVect,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<PVECT>              & NewVect)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldVect,
                    *NewCommDev,
                    *NewVect);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
reduceCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const PVECT        & OldVect,
                   const communicator & NewCommDev,
                         PVECT        & NewVect)
{
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() < OldCommDev.size());
    
    //Alloc
    pMapItemShare mItem;
    
    NewVect.clear();
    NewVect.resize(OldVect.size());
    
    //Main loop
    for(UInt i=1; i <= OldVect.size(); ++i)
    {
      mItem = OldVect.getMapL(i);
      mItem.setPid(NewCommDev.rank());
      
      NewVect.getMapL(i)  = mItem;
      NewVect.getDataL(i) = OldVect.getDataL(i);
    }
    
    NewVect.updateFinder();
  }
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const PVECT>        & OldVect,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<PVECT>              & NewVect)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldVect,
                    *NewCommDev,
                    *NewVect);
}

template<typename DATA>
void
pVectGlobalManip<DATA,pMapItemShare>::
expandCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const PVECT        & OldVect,
                   const communicator & NewCommDev,
                         PVECT        & NewVect)
{
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() > OldCommDev.size());
    
    //Alloc
    pMapItemShare mItem;
    
    NewVect.clear();
    NewVect.resize(OldVect.size());
    
    //Main loop
    for(UInt i=1; i <= OldVect.size(); ++i)
    {
      mItem = OldVect.getMapL(i);
      mItem.setPid(NewCommDev.rank());
      
      NewVect.getMapL(i)  = mItem;
      NewVect.getDataL(i) = OldVect.getDataL(i);
    }
  }
}


//_________________________________________________________________________________________________
// OTHER FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DATA>
size_t
pVectGlobalManip<DATA,pMapItemShare>::
memSize() const
{
  return(sizeof(bool));
}

#endif
