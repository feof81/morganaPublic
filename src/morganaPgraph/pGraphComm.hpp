/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PGRAPHCOMM_HPP
#define PGRAPHCOMM_HPP

#include "Teuchos_RCP.hpp"
#include <Teuchos_RCPDecl.hpp> 

#include "pGraph.hpp"
#include "pVectManip.hpp"
#include "pGraphManip.hpp"
#include "pMapItemSendRecv.h"
#include "pMapGlobalManip.h"
#include "pVectGlobalManip.hpp"

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

using namespace std;
using namespace boost::mpi;


/*! Parallel graph communication */
template<typename ITEM, typename ROWMAP, typename COLMAP> class pGraphComm
{
  /*! @name Typedefs */ //@{
  public:
    typedef pGraph<ITEM,ROWMAP,COLMAP> PGRAPH;
    typedef pVect<ITEM,ROWMAP>        PVECT;
    typedef sVect<ITEM>               DATAVECT;
    typedef pMap<ROWMAP>              PMAP;
    typedef pMap<pMapItemSendRecv>    SENDRECV;
    //@}
  
    /*! @name Internal data and comm-link */ //@{
  public:
    UInt maxSize;
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
    
    
    /*! @name Constructors */ //@{
  public:
    pGraphComm();
    pGraphComm(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    //@}
    
    
    /*! @name Sending functions */ //@{
  public:
    /*! Sending and Receiving
    \param sid sending process
    \param rid receiving process
    \param sendGraph the data to send
    \param rGraph the data received (rGraph is cleared)
    */
    void sendRecv(const UInt & sid, const UInt & rid, const Teuchos::RCP<const PGRAPH> & sendGraph, Teuchos::RCP<PGRAPH> & rGraph) const;
    
    /*! Sending and Receiving
    \param sid sending process
    \param rid receiving process
    \param sendGraph the data to send
    \param rGraph the data received (rGraph is cleared)
    */
    void sendRecv(const UInt & sid, const UInt & rid, const PGRAPH & sendGraph, PGRAPH & rGraph) const;
    
    /*! Sending
    \param sid sending process
    \param rid receiving process
    \param sendGraph the data to send
    */
    void send(const UInt & sid, const UInt & rid, const Teuchos::RCP<const PGRAPH> & sendGraph) const;
    
    /*! Sending
    \param sid sending process
    \param rid receiving process
    \param sendGraph the data to send
    */
    void send(const UInt & sid, const UInt & rid, const PGRAPH & sendGraph) const;
    
    /*! Sending - nonBlocking
    \param sid sending process
    \param rid receiving process
    \param sendGraph the data to send
    */
    void isend(const UInt & sid, const UInt & rid, const Teuchos::RCP<const PGRAPH> & sendGraph, sVect<request> & reqs) const;
    
    /*! Sending - nonBlocking
    \param sid sending process
    \param rid receiving process
    \param sendGraph the data to send
    */
    void isend(const UInt & sid, const UInt & rid, const PGRAPH & sendGraph, sVect<request> & reqs) const;
    
    /*! Receiving
    \param sid sending process
    \param rid receiving process
    \param rGraph the data received (rGraph is cleared)
    */
    void recv(const UInt & sid, const UInt & rid, Teuchos::RCP<PGRAPH> & rGraph) const;
    
    /*! Receiving
    \param sid sending process
    \param rid receiving process
    \param rGraph the data received (rGraph is cleared)
    */
    void recv(const UInt & sid, const UInt & rid, PGRAPH & rGraph) const;
    
    /*! Receiving - pseudoNonBlocking
    \param sid sending process
    \param rid receiving process
    \param rGraph the data received (rGraph is cleared)
    */
    void irecv(const UInt & sid, const UInt & rid, Teuchos::RCP<PGRAPH> & rGraph) const;
    
    /*! Receiving - pseudoNonBlocking
    \param sid sending process
    \param rid receiving process
    \param rGraph the data received (rGraph is cleared)
    */
    void irecv(const UInt & sid, const UInt & rid, PGRAPH & rGraph) const;
    
    /*! Merging the received data to \c rMap
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    \param rGraph the data received (rGraph is merged with the arriving data)
    */
    void merge(const UInt & sid, const UInt & rid, const PGRAPH & sendGraph, PGRAPH & rGraph) const;
    
    /*! Merging the received data to \c rMap
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    \param rGraph the data received (rGraph is merged with the arriving data)
    */
    void merge(const UInt & sid, const UInt & rid, const Teuchos::RCP<const PGRAPH> & sendGraph, Teuchos::RCP<PGRAPH> & rGraph) const;
    //@}
    
    /*! @name Memory functions */ //@{
  public:
    size_t memSize() const;
    //@}
};



//_______________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
pGraphComm<ITEM,ROWMAP,COLMAP>::
pGraphComm()
{
  maxSize       = mpiMaxSize / (sizeof(ITEM) + sizeof(ROWMAP));
  commDevLoaded = false;
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
pGraphComm<ITEM,ROWMAP,COLMAP>::
pGraphComm(const Teuchos::RCP<const communicator> & CommDev)
{
   maxSize       = mpiMaxSize / (sizeof(ITEM) + sizeof(ROWMAP));
   commDevLoaded = true;
   commDev       = CommDev;
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphComm<ITEM,ROWMAP,COLMAP>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}



//_______________________________________________________________________________________________________
// SEND FUNCTIONS
//-------------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphComm<ITEM,ROWMAP,COLMAP>::
sendRecv(const UInt & sid, const UInt & rid, const Teuchos::RCP<const PGRAPH> & sendGraph, Teuchos::RCP<PGRAPH> & rGraph) const
{
  assert(commDevLoaded);
  sendRecv(sid,rid,*sendGraph,*rGraph);
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphComm<ITEM,ROWMAP,COLMAP>::
sendRecv(const UInt & sid, const UInt & rid, const PGRAPH & sendGraph, PGRAPH & rGraph) const
{
  assert(commDevLoaded);
  
  typedef pVectComm<ITEM,ROWMAP> PVECTCOMM;
  typedef pMapComm<COLMAP>       PMAPCOMM;
  
  //Row communication
  PVECTCOMM rowCommunicator(commDev);
  rowCommunicator.sendRecv(sid,rid,sendGraph,rGraph);
  
  //Col communication
  PMAPCOMM colCommunicator(commDev);
  pMap<COLMAP> tempColMap;
  
  colCommunicator.sendRecv(sid,rid,sendGraph.getColMap(),tempColMap);
  rGraph.setColMap(tempColMap);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphComm<ITEM,ROWMAP,COLMAP>::
send(const UInt & sid, const UInt & rid, const Teuchos::RCP<const PGRAPH> & sendGraph) const
{
  assert(commDevLoaded);
  send(sid,rid,*sendGraph);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphComm<ITEM,ROWMAP,COLMAP>::
send(const UInt & sid, const UInt & rid, const PGRAPH & sendGraph) const
{
  assert(commDevLoaded);
  
  typedef pVectComm<ITEM,ROWMAP> PVECTCOMM;
  typedef pMapComm<COLMAP>       PMAPCOMM;
  
  //Row communication
  PVECTCOMM rowCommunicator(commDev);
  rowCommunicator.send(sid,rid,sendGraph);
  
  //Col communication
  PMAPCOMM colCommunicator(commDev);  
  colCommunicator.send(sid,rid,sendGraph.getColMap());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphComm<ITEM,ROWMAP,COLMAP>::
isend(const UInt & sid, const UInt & rid, const Teuchos::RCP<const PGRAPH> & sendGraph, sVect<request> & reqs) const
{
  assert(commDevLoaded);
  isend(sid,rid,*sendGraph,reqs);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphComm<ITEM,ROWMAP,COLMAP>::
isend(const UInt & sid, const UInt & rid, const PGRAPH & sendGraph, sVect<request> & reqs) const
{
  assert(commDevLoaded);
  
  typedef pVectComm<ITEM,ROWMAP> PVECTCOMM;
  typedef pMapComm<COLMAP>       PMAPCOMM;
  
  //Row communication
  sVect<request> reqsRows;
  PVECTCOMM rowCommunicator(commDev);
  rowCommunicator.isend(sid,rid,sendGraph,reqsRows);
  
  //Col communication
  sVect<request> reqsCols;
  PMAPCOMM colCommunicator(commDev);  
  colCommunicator.isend(sid,rid,sendGraph.getColMap(),reqsCols);
  
  //Collecting reqs
  reqs.reserve(reqsRows.size() + reqsCols.size());
  
  for(UInt i=1; i <= reqsRows.size(); ++i)
  { reqs.push_back(reqsRows(i)); }
  
  for(UInt i=1; i <= reqsCols.size(); ++i)
  { reqs.push_back(reqsCols(i)); }
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphComm<ITEM,ROWMAP,COLMAP>::
recv(const UInt & sid, const UInt & rid, Teuchos::RCP<PGRAPH> & rGraph) const
{
  assert(commDevLoaded);
  recv(sid,rid,*rGraph);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphComm<ITEM,ROWMAP,COLMAP>::
recv(const UInt & sid, const UInt & rid, PGRAPH & rGraph) const
{
  assert(commDevLoaded);
  
  typedef pVectComm<ITEM,ROWMAP> PVECTCOMM;
  typedef pMapComm<COLMAP>       PMAPCOMM;
  
  //Row communication
  PVECTCOMM rowCommunicator(commDev);
  rowCommunicator.recv(sid,rid,rGraph);
  
  //Col communication
  PMAPCOMM colCommunicator(commDev);
  pMap<COLMAP> tempColMap;
  
  colCommunicator.recv(sid,rid,tempColMap);
  rGraph.setColMap(tempColMap);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphComm<ITEM,ROWMAP,COLMAP>::
irecv(const UInt & sid, const UInt & rid, Teuchos::RCP<PGRAPH> & rGraph) const
{
  assert(commDevLoaded);
  irecv(sid,rid,*rGraph);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphComm<ITEM,ROWMAP,COLMAP>::
irecv(const UInt & sid, const UInt & rid, PGRAPH & rGraph) const
{
  assert(commDevLoaded);
  
  typedef pVectComm<ITEM,ROWMAP> PVECTCOMM;
  typedef pMapComm<COLMAP>       PMAPCOMM;
  
  //Row communication
  PVECTCOMM rowCommunicator(commDev);
  rowCommunicator.irecv(sid,rid,rGraph);
  
  //Col communication
  PMAPCOMM colCommunicator(commDev);
  pMap<COLMAP> tempColMap;
  
  colCommunicator.irecv(sid,rid,tempColMap);
  rGraph.setColMap(tempColMap);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphComm<ITEM,ROWMAP,COLMAP>::
merge(const UInt & sid, const UInt & rid, const PGRAPH & sendGraph, PGRAPH & rGraph) const
{
  assert(commDevLoaded);
  merge(sid,rid,*sendGraph,*rGraph);
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphComm<ITEM,ROWMAP,COLMAP>::
merge(const UInt & sid, const UInt & rid, const Teuchos::RCP<const PGRAPH> & sendGraph, Teuchos::RCP<PGRAPH> & rGraph) const
{
  assert(commDevLoaded);
  
  typedef pVectComm<ITEM,ROWMAP> PVECTCOMM;
  typedef pMapComm<COLMAP>       PMAPCOMM;
  
  //Temporary graph
  PGRAPH tempGraph;
  
  //Row communication
  PVECTCOMM rowCommunicator(commDev);
  rowCommunicator.sendRecv(sid,rid,sendGraph,tempGraph);
  
  //Col communication
  PMAPCOMM colCommunicator(commDev);
  pMap<COLMAP> tempColMap;
  
  colCommunicator.sendRecv(sid,rid,sendGraph.getColMap(),tempColMap);
  tempGraph.setColMap(tempColMap);
  
  //Merging
  pGraphManip<ITEM,ROWMAP,COLMAP> merger;
  merger.mergeGraph(rGraph,tempGraph);
}


//_________________________________________________________________________________________________
// MEMORY FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
size_t
pGraphComm<ITEM,ROWMAP,COLMAP>::
memSize() const
{
  return(sizeof(UInt) + sizeof(bool));
}

#endif
