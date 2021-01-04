/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PMAPCOMM_HPP
#define PMAPCOMM_HPP

#include "Teuchos_RCP.hpp"
#include <Teuchos_RCPDecl.hpp> 


#include "pMpiOptimization.hpp"

#include "pMap.hpp"
#include "pMapItem.h"
#include "pMapManip.hpp"

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

using namespace std;
using namespace boost::mpi;


/*! Class for the communication of \c pMap */
template<typename ITEM> class pMapComm
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMap<ITEM> MAP;
    //@}
    
    
    /*! @name Internal data and comm-link */ //@{
  public:
    UInt maxSize;
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
    
    
    /*! @name Constructors */ //@{
  public:
    /*! Constructor */
    pMapComm();
    
    /*! Constructor */
    pMapComm(const Teuchos::RCP<const communicator> & CommDev);
    
    /*! Setting of the communication device */
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    //@}
    
    
    /*! @name Send and Recv - non-pending-comm */ //@{
  public:
    /*! Sending and Receiving
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    \param rMap the data received (rMap is cleared)
    */
    void sendRecv(const UInt & sid, const UInt & rid, const Teuchos::RCP<const MAP> & sMap, Teuchos::RCP<MAP> & rMap) const;
    
    /*! Sending and Receiving
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    \param rMap the data received (rMap is cleared)
    */
    void sendRecv(const UInt & sid, const UInt & rid, const MAP & sMap, MAP & rMap) const;
    
    /*! Sending
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    */
    void send(const UInt & sid, const UInt & rid, const Teuchos::RCP<const MAP> & sMap) const;
    
    /*! Sending
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    */
    void send(const UInt & sid, const UInt & rid, const MAP & sMap) const;
    
    /*! Not-sync sending
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    \param reqs the sending status vector
    */
    void isend(const UInt & sid, const UInt & rid, const Teuchos::RCP<const MAP> & sMap, sVect<request> & reqs) const;
    
    /*! Not-sync sending
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    \param reqs the sending status vector
    */
    void isend(const UInt & sid, const UInt & rid, const MAP & sMap, sVect<request> & reqs) const;
    
    /*! Receiving
    \param sid sending process
    \param rid receiving process
    \param rMap the data received (rMap is cleared)
    */
    void recv(const UInt & sid, const UInt & rid, Teuchos::RCP<MAP> & rMap) const;
    
    /*! Receiving
    \param sid sending process
    \param rid receiving process
    \param rMap the data received (rMap is cleared)
    */
    void recv(const UInt & sid, const UInt & rid, MAP & rMap) const;
    
    /*! Not-sync receiving
    \param sid sending process
    \param rid receiving process
    \param rMap the data received (rMap is cleared)
    */
    void irecv(const UInt & sid, const UInt & rid, Teuchos::RCP<MAP> & rMap) const;
    
    /*! Not-sync receiving
    \param sid sending process
    \param rid receiving process
    \param rMap the data received (rMap is cleared)
    */
    void irecv(const UInt & sid, const UInt & rid, MAP & rMap) const;
    
    /*! Merging the received data to \c rMap
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    \param rMap the data received (rMap is merged with the arriving data)
    */
    void merge(const UInt & sid, const UInt & rid, const MAP & sMap, MAP & rMap) const;
    
    /*! Merging the received data to \c rMap
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    \param rMap the data received (rMap is merged with the arriving data)
    */
    void merge(const UInt & sid, const UInt & rid, const Teuchos::RCP<const MAP> & sMap, Teuchos::RCP<MAP> & rMap) const;
    //@}
    
    
    /*! @name Pattern communications - non-pending-comm */ //@{
  public:
    /*! Normal data distribution (gid-based)
    \param map contains the data to send, then is cleared and loaded with the new data.
    The \c maxGid parameter can be used to force the vector to a different distribution
    on the communicator
    */
    void vectorNormal(Teuchos::RCP<MAP> & map, const UInt & maxGid = 0) const;
    
    /*! Normal data distribution (gid-based)
    \param map contains the data to send, then is cleared and loaded with the new data.
    The \c maxGid parameter can be used to force the vector to a different distribution
    on the communicator
    */
    void vectorNormal(MAP & map, const UInt & maxGid = 0) const;
    
    /*! Pid data distribution (pid-based)
    \param map contains the data to send, then is cleared and loaded with the new data
    */
    void vectorPid(Teuchos::RCP<MAP> & map) const;
    
    /*! Pid data distribution (pid-based)
    \param map contains the data to send, then is cleared and loaded with the new data
    */
    void vectorPid(MAP & map) const;
    //@}
};



//_______________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------------
template<typename ITEM>
pMapComm<ITEM>::
pMapComm()
{
  maxSize       = mpiMaxSize / sizeof(ITEM);
  commDevLoaded = false;
}

template<typename ITEM>
pMapComm<ITEM>::
pMapComm(const Teuchos::RCP<const communicator> & CommDev)
{
   maxSize       = mpiMaxSize / sizeof(ITEM);
   commDevLoaded = true;
   commDev       = CommDev;
}

template<typename ITEM>
void
pMapComm<ITEM>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}



//_______________________________________________________________________________________________________
// SEND AND RECV - NON-PENDING-COMM
//-------------------------------------------------------------------------------------------------------
template<typename ITEM>
void
pMapComm<ITEM>::
sendRecv(const UInt & sid, const UInt & rid, const Teuchos::RCP<const MAP> & sMap, Teuchos::RCP<MAP> & rMap) const
{
  assert(commDevLoaded);
  assert(sMap.total_count() > 0);
  assert(rMap.total_count() > 0);
  
  sendRecv(sid,rid,*sMap,*rMap);
}

#ifdef NOCOMMBUFFER
template<typename ITEM>
void
pMapComm<ITEM>::
sendRecv(const UInt & sid, const UInt & rid, const MAP & sMap, MAP & rMap) const
{
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    //Data transmission
    sVect<request> reqs(1);
    reqs(1) = commDev->isend(rid,1,sMap);
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Data transmission
    sVect<request> reqs(1);
    
    rMap.clear();
    reqs(1) = commDev->irecv(sid,1,rMap);
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
}
#else
template<typename ITEM>
void
pMapComm<ITEM>::
sendRecv(const UInt & sid, const UInt & rid, const MAP & sMap, MAP & rMap) const
{
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    //Segmentazione dati
    pMapManip<ITEM> manipulator;
    sVect<MAP> segmented;
    
    manipulator.segmentationSimple(segmented,maxSize,sMap);
    UInt n = segmented.size();
    assert(n <= mpiLayer);
    
    //Data transmission
    commDev->send(rid,1,n);
    
    //Inizio trasmissione
    sVect<request> reqs(n);
    
    for(UInt i=1; i<= n; ++i)
    { 
      reqs(i) = commDev->isend(rid,i,segmented(i));
    }
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Ricezione dimensioni ed allocazione mappa ricezione
    UInt n;
    commDev->recv(sid,1,n);
    
    sVect<MAP> segmented(n);
    
    //Data transmission
    sVect<request> reqs(n);
    
    for(UInt i=1; i<= n; ++i)
    { 
      reqs(i) = commDev->irecv(sid,i,segmented(i));
    }
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
    
    //Ricostruzione dati
    rMap.clear();
    pMapManip<ITEM> manipulator;
    manipulator.mergeSimple(segmented,rMap);
  }
}
#endif


template<typename ITEM>
void
pMapComm<ITEM>::
send(const UInt & sid, const UInt & rid, const Teuchos::RCP<const MAP> & sMap) const
{
  assert(commDevLoaded);
  assert(sMap.total_count() > 0);
  
  send(sid,rid,*sMap);
}

#ifdef NOCOMMBUFFER
template<typename ITEM>
void
pMapComm<ITEM>::
send(const UInt & sid, const UInt & rid, const MAP & sMap) const
{
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {    
    //Send data
    sVect<request> reqs(1);
    reqs(1) = commDev->isend(rid,1,sMap);
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
}
#else
template<typename ITEM>
void
pMapComm<ITEM>::
send(const UInt & sid, const UInt & rid, const MAP & sMap) const
{
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    //Segmentazione dati
    pMapManip<ITEM> manipulator;
    sVect<MAP> segmented;
    
    manipulator.segmentationSimple(segmented,maxSize,sMap);
    UInt n = segmented.size();
    assert(n <= mpiLayer);
    
    //Trasmissione dimensioni
    commDev->send(rid,1,n);
    
    //Inizio trasmissione
    sVect<request> reqs(n);
    
    for(UInt i=1; i<= n; ++i)
    { 
      reqs(i) = commDev->isend(rid,i,segmented(i));
    }
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
}
#endif


template<typename ITEM>
void
pMapComm<ITEM>::
isend(const UInt & sid, const UInt & rid, const Teuchos::RCP<const MAP> & sMap, sVect<request> & reqs) const
{
  assert(commDevLoaded);
  assert(sMap.total_count() > 0);
  
  isend(sid,rid,*sMap,reqs);
}

#ifdef NOCOMMBUFFER
template<typename ITEM>
void
pMapComm<ITEM>::
isend(const UInt & sid, const UInt & rid, const MAP & sMap, sVect<request> & reqs) const
{
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {    
    //Inizio trasmissione
    reqs.reserve(1);
    reqs.push_back( commDev->isend(rid,1,sMap) );
  }
}
#else
template<typename ITEM>
void
pMapComm<ITEM>::
isend(const UInt & sid, const UInt & rid, const MAP & sMap, sVect<request> & reqs) const
{
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    //Segmentazione dati
    pMapManip<ITEM> manipulator;
    sVect<MAP> segmented;
    
    manipulator.segmentationSimple(segmented,maxSize,sMap);
    UInt n = segmented.size();
    assert(n <= mpiLayer);
    
    //Trasmissione dimensioni
    commDev->send(rid,1,n);
    
    //Inizio trasmissione
    reqs.reserve(n);
    
    for(UInt i=1; i<= n; ++i)
    { 
      reqs.push_back( commDev->isend(rid,i,segmented(i)) );
    }
  }
}
#endif


template<typename ITEM>
void
pMapComm<ITEM>::
recv(const UInt & sid, const UInt & rid, Teuchos::RCP<MAP> & rMap) const
{
  assert(commDevLoaded);
  assert(rMap.total_count() > 0);
  
  recv(sid,rid,*rMap);
}

#ifdef NOCOMMBUFFER
template<typename ITEM>
void
pMapComm<ITEM>::
recv(const UInt & sid, const UInt & rid, MAP & rMap) const
{
  assert(commDevLoaded);
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Recv data
    rMap.clear();
    
    sVect<request> reqs(1);
    reqs(1) = commDev->irecv(sid,1,rMap);
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
}
#else
template<typename ITEM>
void
pMapComm<ITEM>::
recv(const UInt & sid, const UInt & rid, MAP & rMap) const
{
  assert(commDevLoaded);
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Ricezione dimensioni ed allocazione mappa ricezione
    UInt n;
    commDev->recv(sid,1,n);
    
    sVect<MAP> segmented(n);
    
    //Inizio trasmissione
    sVect<request> reqs(n);
    
    for(UInt i=1; i<= n; ++i)
    { 
      reqs(i) = commDev->irecv(sid,i,segmented(i));
    }
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
    
    //Ricostruzione dati
    rMap.clear();
    pMapManip<ITEM> manipulator;
    manipulator.mergeSimple(segmented,rMap);
  }
}
#endif


template<typename ITEM>
void
pMapComm<ITEM>::
irecv(const UInt & sid, const UInt & rid, Teuchos::RCP<MAP> & rMap) const
{
  assert(commDevLoaded);
  assert(rMap.total_count() > 0);
  
  irecv(sid,rid,*rMap);
}

#ifdef NOCOMMBUFFER
template<typename ITEM>
void
pMapComm<ITEM>::
irecv(const UInt & sid, const UInt & rid, MAP & rMap) const
{
  assert(commDevLoaded);
  
  //Reciver
  if(commDev->rank() == int(rid))
  {    
    //Inizio trasmissione
    rMap.clear();
    
    sVect<request> reqs(1);
    reqs(1) = commDev->irecv(sid,1,rMap);

    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
}
#else
template<typename ITEM>
void
pMapComm<ITEM>::
irecv(const UInt & sid, const UInt & rid, MAP & rMap) const
{
  assert(commDevLoaded);
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Ricezione dimensioni ed allocazione mappa ricezione
    UInt n;
    commDev->recv(sid,1,n);
    
    sVect<MAP> segmented(n);
    
    //Inizio trasmissione
    sVect<request> reqs(n);
    
    for(UInt i=1; i<= n; ++i)
    { 
      reqs(i) = commDev->irecv(sid,i,segmented(i));
    }
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
    
    //Ricostruzione dati
    rMap.clear();
    pMapManip<ITEM> manipulator;
    manipulator.mergeSimple(segmented,rMap);
  }
}
#endif


template<typename ITEM>
void
pMapComm<ITEM>::
merge(const UInt & sid, const UInt & rid, const Teuchos::RCP<const MAP> & sMap, Teuchos::RCP<MAP> & rMap) const
{
  assert(commDevLoaded);
  assert(sMap.total_count() > 0);
  assert(rMap.total_count() > 0);
  
  merge(sid,rid,*sMap,*rMap);
}

#ifdef NOCOMMBUFFER
template<typename ITEM>
void
pMapComm<ITEM>::
merge(const UInt & sid, const UInt & rid, const MAP & sMap, MAP & rMap) const
{
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {    
    //Inizio trasmissione
    sVect<request> reqs(1);
    reqs(1) = commDev->isend(rid,1,sMap);

    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    sVect<MAP> recMap(1);
    sVect<request> reqs(1);
    reqs(1) = commDev->irecv(sid,1,recMap(1));
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
    
    //Ricostruzione dati
    pMapManip<ITEM> manipulator;
    manipulator.mergeSimple(recMap,rMap);
  }
}
#else
template<typename ITEM>
void
pMapComm<ITEM>::
merge(const UInt & sid, const UInt & rid, const MAP & sMap, MAP & rMap) const
{
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    //Segmentazione dati
    pMapManip<ITEM> manipulator;
    sVect<MAP> segmented;
    
    manipulator.segmentationSimple(segmented,maxSize,sMap);
    UInt n = segmented.size();
    assert(n <= mpiLayer);
    
    //Trasmissione dimensioni
    commDev->send(rid,1,n);
    
    //Inizio trasmissione
    sVect<request> reqs(n);
    
    for(UInt i=1; i<= n; ++i)
    { 
      reqs(i) = commDev->isend(rid,i,segmented(i));
    }
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Ricezione dimensioni ed allocazione mappa ricezione
    UInt n;
    commDev->recv(sid,1,n);
    
    sVect<MAP> segmented(n);
    
    //Inizio trasmissione
    sVect<request> reqs(n);
    
    for(UInt i=1; i<= n; ++i)
    { 
      reqs(i) = commDev->irecv(sid,i,segmented(i));
    }
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
    
    //Ricostruzione dati
    pMapManip<ITEM> manipulator;
    manipulator.mergeSimple(segmented,rMap);
  }
}
#endif



//_______________________________________________________________________________________________________
// PATTERN COMMUNICATIONS - NON-PENDING-COMM
//-------------------------------------------------------------------------------------------------------
template<typename ITEM>
void
pMapComm<ITEM>::
vectorNormal(Teuchos::RCP<MAP> & map, const UInt & maxGid) const
{
  assert(commDevLoaded);
  assert(map.total_count() > 0);
  
  vectorNormal(*map,maxGid);
}

template<typename ITEM>
void
pMapComm<ITEM>::
vectorNormal(MAP & map, const UInt & maxGid) const
{
  assert(commDevLoaded);
  
  //Alloc
  UInt proc    = commDev->rank();
  UInt numProc = commDev->size();
  
  //Maximum gid determination
  pMapManip<ITEM>  manipulator;
  UInt globalMax = max(manipulator.getMaxGid(map), maxGid);
  
  sVect<UInt> values;
  all_gather(*commDev,globalMax,values);
  
  for(UInt i=1; i <= values.size(); ++i)
  {
    globalMax = max(globalMax,values(i));
  }
  
  //Segmenting
  sVect<MAP>  sendMaps;
  sVect<MAP>  recvMaps;
  MAP         tempMap;
  
  manipulator.segmentationNormal(sendMaps,numProc,globalMax,map);
  
  
  //Sending
  sVect<request> sendReqs;
  
  for(UInt k=0; k < numProc; ++k)
  {
    if(k != proc)
    {
      isend(proc,k,sendMaps[k],sendReqs);
    }
  }
  
  //Reciving 
  recvMaps.reserve(numProc);
  
  for(UInt k=0; k < numProc; ++k)
  {
    if(k != proc)
    {
      irecv(k,proc,tempMap);
      recvMaps.push_back(tempMap);
    }
  }
  
  //Waiting
  wait_all(sendReqs.begin(),sendReqs.end());
  
  //Merging
  recvMaps.push_back(sendMaps[proc]);
  
  map.clear();
  manipulator.mergeSimple(recvMaps,map);
}

template<typename ITEM>
void
pMapComm<ITEM>::
vectorPid(Teuchos::RCP<MAP> & map) const
{
  assert(commDevLoaded);
  assert(map.total_count() > 0);
  
  vectorPid(*map);
}

template<typename ITEM>
void
pMapComm<ITEM>::
vectorPid(MAP & map) const
{
  assert(commDevLoaded);
  
  //Alloc
  UInt proc    = commDev->rank();
  UInt numProc = commDev->size();
  sVect<MAP>  sendMaps;
  sVect<MAP>  recvMaps;
  MAP              tempMap;
  
  //Segmenting
  pMapManip<ITEM>  manipulator;
  manipulator.segmentationPid(sendMaps,numProc,map);
  
  //Sending
  sVect<request> sendReqs;
  
  for(UInt k=0; k < numProc; ++k)
  {
    if(k != proc)
    {
      isend(proc,k,sendMaps[k],sendReqs);
    }
  }
  
  //Reciving 
  recvMaps.reserve(numProc);
  
  for(UInt k=0; k < numProc; ++k)
  {
    if(k != proc)
    {
      irecv(k,proc,tempMap);
      recvMaps.push_back(tempMap);
    }
  }
  
  //Waiting
  wait_all(sendReqs.begin(),sendReqs.end());
  
  //Merging
  recvMaps.push_back(sendMaps[proc]);
  
  map.clear();
  manipulator.mergeSimple(recvMaps,map);
}

#endif
