/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PVECTCOMM_HPP
#define PVECTCOMM_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_RCPDecl.hpp" 

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"

#include "simpleFormats.hpp"
#include "traitsMpiOptimization.hpp"

#include "morganaContainer.h"

#include "pVect.hpp"
#include "pVectManip.hpp"
#include "pMapItemSendRecv.h"
#include "pMapGlobalManip.h"

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

using namespace std;
using namespace boost::mpi;

//_________________________________________________________________________________________________
// COMM STRUCTS
//-------------------------------------------------------------------------------------------------

/*! Class for exchanging data for pending comm - Parallel Vector PendingS */
template<typename DATA, typename MAPITEM> 
class pVectPendings
{
  public:
    typedef pVect<DATA,MAPITEM> PVECT;
  
  public :
    sVect<request> reqs;
    sVect<PVECT>   bufVects;
};

/*! Class for exchanging data for recursive pending comm - Parallel Vector Recursives */
template<typename DATA, typename MAPITEM> 
class pVectRecursives
{
  public :
    sVect<request>      reqs;
    boost::mpi::content bufCont;
};

/*! Class for recursive pVect update */
template<typename DATA, typename MAPITEM> 
class pVectUpdateRecursive
{
  public:
    typedef pVect<DATA,MAPITEM>           PVECT;
    typedef pVectRecursives<DATA,MAPITEM> PVRS;
    
  public:
    PVECT        rVect;
    sVect<PVECT> sendSegments;
    sVect<PVECT> recvSegments;
    sVect<PVRS>  sendPvrs;
    sVect<PVRS>  recvPvrs;
};



//_________________________________________________________________________________________________
// COMM CLASS
//-------------------------------------------------------------------------------------------------

/*! Class for the communication of \c pVect */
template<typename DATA, typename MAPITEM> class pVectComm
{
    /*! @name Typedefs */ //@{
  public:
    typedef pVect<DATA,MAPITEM>           PVECT;
    typedef sVect<DATA>                   DATAVECT;
    typedef pMap<MAPITEM>                 PMAP;
    typedef pMap<pMapItemSendRecv>        SENDRECV;
    typedef pVectPendings<DATA,MAPITEM>   PVPS;
    typedef pVectRecursives<DATA,MAPITEM> PVRS;
    //@}
    
    
    /*! @name Internal data and comm-link */ //@{
  public:
    UInt maxSize;
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
    
    
    /*! @name Constructors */ //@{
  public:
    pVectComm();
    pVectComm(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    //@}
    
    
    /*! @name Send and Recv - non-pending-comm */ //@{
  public:
    /*! Sending and Receiving
    \param sid sending process
    \param rid receiving process
    \param sendVect the data to send
    \param rVect the data received (rVect is cleared)
    */
    void sendRecv(const UInt & sid,
                  const UInt & rid,
                  const Teuchos::RCP<const PVECT> & sendVect,
                        Teuchos::RCP<PVECT>       & rVect) const;
    
    /*! Sending and Receiving
    \param sid sending process
    \param rid receiving process
    \param sendVect the data to send
    \param rVect the data received (rVect is cleared)
    */
    void sendRecv(const UInt & sid,
                  const UInt & rid,
                  const PVECT & sendVect,
                        PVECT & rVect) const;
        
    /*! Sending
    \param sid sending process
    \param rid receiving process
    \param sendVect the data to send
    */
    void send(const UInt & sid,
              const UInt & rid,
              const Teuchos::RCP<const PVECT> & sendVect) const;
    
    /*! Sending
    \param sid sending process
    \param rid receiving process
    \param sendVect the data to send
    */
    void send(const UInt & sid,
              const UInt & rid,
              const PVECT & sendVect) const;
    
    /*! Not-sync sending
    \param sid sending process
    \param rid receiving process
    \param sendVect the data to send
    \param reqs the sending status vector
    */
    void isend(const UInt & sid,
               const UInt & rid,
               const Teuchos::RCP<const PVECT> & sendVect,
               sVect<request> & reqs) const;
    
    /*! Not-sync sending
    \param sid sending process
    \param rid receiving process
    \param sendVect the data to send
    \param reqs the sending status vector
    */
    void isend(const UInt & sid,
               const UInt & rid,
               const PVECT & sendVect,
               sVect<request> & reqs) const;
    
    /*! Receiving
    \param sid sending process
    \param rid receiving process
    \param rVect the data received (rVect is cleared)
    */
    void recv(const UInt & sid,
              const UInt & rid,
              Teuchos::RCP<PVECT> & rVect) const;
    
    /*! Receiving
    \param sid sending process
    \param rid receiving process
    \param rVect the data received (rVect is cleared)
    */
    void recv(const UInt & sid,
              const UInt & rid,
              PVECT & rVect) const;
    
    /*! Not-sync receiving - dummy
    \param sid sending process
    \param rid receiving process
    \param rVect the data received (rVect is cleared)
    */
    void irecv(const UInt & sid,
               const UInt & rid,
               Teuchos::RCP<PVECT> & rVect) const;
    
    /*! Not-sync receiving - dummy
    \param sid sending process
    \param rid receiving process
    \param rVect the data received (rVect is cleared)
    */
    void irecv(const UInt & sid,
               const UInt & rid,
               PVECT & rVect) const;
    
    /*! Merging the received data to \c rVect
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    \param rVect the data received (rVect is merged with the arriving data)
    */
    void merge(const UInt & sid,
               const UInt & rid,
               const PVECT & sendVect,
                     PVECT & rVect) const;
    
    /*! Merging the received data to \c rVect
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    \param rVect the data received (rVect is merged with the arriving data)
    */
    void merge(const UInt & sid,
               const UInt & rid,
               const Teuchos::RCP<const PVECT> & sendVect,
               Teuchos::RCP<PVECT>             & rVect) const;
    //@}
  

    /*! @name Send and Recv - pending-comm */ //@{
  public:
    /*! Not-sync sending
    \param sid sending process
    \param rid receiving process
    \param sendVect the data to send
    \param pvps pending commData
    */
    void sendI(const UInt                & sid,
               const UInt                & rid,
               const Teuchos::RCP<PVECT> & sendVect,
                     PVPS                & pvps,
               const UInt                & channel = 2) const;
   
    /*! Not-sync sending
    \param sid sending process
    \param rid receiving process
    \param sendVect the data to send
    \param pvps pending commData
    */
    void sendO(const UInt                & sid,
               const UInt                & rid,
               const Teuchos::RCP<PVECT> & sendVect,
                     PVPS                & pvps,
               const UInt                & channel = 2) const;
    
    /*! Not-sync sending
    \param sid sending process
    \param rid receiving process
    \param sendVect the data to send
    \param pvps pending commData
    */
    void sendI(const UInt  & sid,
               const UInt  & rid,
               const PVECT & sendVect,
                     PVPS  & pvps,
               const UInt  & channel = 2) const;
    
    /*! Not-sync sending
    \param sid sending process
    \param rid receiving process
    \param sendVect the data to send
    \param pvps pending commData
    */
    void sendO(const UInt  & sid,
               const UInt  & rid,
               const PVECT & sendVect,
                     PVPS  & pvps,
               const UInt  & channel = 2) const;
     
    /*! Not-sync receiving
    \param sid sending process
    \param rid receiving process
    \param rVect the data received (rVect is cleared)
    \param pvps pending commData
    */
    void recvI(const UInt                & sid,
               const UInt                & rid,
                     Teuchos::RCP<PVECT> & rVect,
                     PVPS                & pvps,
               const UInt                & channel = 2) const;
     
    /*! Not-sync receiving
    \param sid sending process
    \param rid receiving process
    \param rVect the data received (rVect is cleared)
    \param pvps pending commData
    */
    void recvO(const UInt                & sid,
               const UInt                & rid,
                     Teuchos::RCP<PVECT> & rVect,
                     PVPS                & pvps,
               const UInt                & channel = 2) const;
    
    /*! Not-sync receiving
    \param sid sending process
    \param rid receiving process
    \param rVect the data received (rVect is cleared)
    \param pvps pending commData
    */
    void recvI(const UInt & sid,
               const UInt & rid,
                    PVECT & rVect,
                     PVPS & pvps,
               const UInt & channel = 2) const;
     
    /*! Not-sync receiving
    \param sid sending process
    \param rid receiving process
    \param rVect the data received (rVect is cleared)
    \param pvps pending commData
    */
    void recvO(const UInt & sid,
               const UInt & rid,
                    PVECT & rVect,
                     PVPS & pvps,
               const UInt & channel = 2) const;
    
    /*! Sending and Receiving
    \param sid sending process
    \param rid receiving process
    \param sendVect the data to send
    \param rVect the data received (rVect is cleared)
    \param commSegments data vectors for handling the parts of the communication
    \param reqs wait requests
    */
    void sendRecvI(const UInt                      & sid,
                   const UInt                      & rid,
                   const Teuchos::RCP<const PVECT> & sendVect,
                         Teuchos::RCP<PVECT>       & rVect,
                         PVPS                      & pvps,
                   const UInt                      & channel = 2) const;

     /*! Sending and Receiving
    \param sid sending process
    \param rid receiving process
    \param sendVect the data to send
    \param rVect the data received (rVect is cleared)
    \param commSegments data vectors for handling the parts of the communication
    \param reqs wait requests
    */
    void sendRecvO(const UInt                      & sid,
                   const UInt                      & rid,
                   const Teuchos::RCP<const PVECT> & sendVect,
                         Teuchos::RCP<PVECT>       & rVect,
                         PVPS                      & pvps,
                   const UInt                      & channel = 2) const;
    
    /*! Sending and Receiving - start communications
    \param sid sending process
    \param rid receiving process
    \param sendVect the data to send
    \param rVect the data received (rVect is cleared)
    \param commSegments data vectors for handling the parts of the communication
    \param reqs wait requests
    */
    void sendRecvI(const UInt  & sid,
                   const UInt  & rid,
                   const PVECT & sendVect,
                         PVECT & rVect,
                         PVPS  & pvps,
                   const UInt  & channel = 2) const;
 
    /*! Sending and Receiving - stop communications
    \param sid sending process
    \param rid receiving process
    \param sendVect the data to send
    \param rVect the data received (rVect is cleared)
    \param commSegments data vectors for handling the parts of the communication
    \param reqs wait requests
    */
    void sendRecvO(const UInt  & sid,
                   const UInt  & rid,
                   const PVECT & sendVect,
                         PVECT & rVect,
                         PVPS  & pvps,
                   const UInt  & channel = 2) const;
       
    /*! Merging the received data to \c rVect
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    \param rVect the data received (rVect is merged with the arriving data)
    \param commSegments data vectors for handling the parts of the communication
    \param reqs wait requests
    */
    void mergeI(const UInt  & sid,
                const UInt  & rid,
                const PVECT & sendVect,
                      PVECT & rVect,
                      PVPS  & pvps,
                const UInt  & channel = 2) const;
     
    /*! Merging the received data to \c rVect
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    \param rVect the data received (rVect is merged with the arriving data)
    \param commSegments data vectors for handling the parts of the communication
    \param reqs wait requests
    */
    void mergeO(const UInt  & sid,
                const UInt  & rid,
                const PVECT & sendVect,
                      PVECT & rVect,
                      PVPS  & pvps,
                const UInt  & channel = 2) const;
    
    /*! Merging the received data to \c rVect
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    \param rVect the data received (rVect is merged with the arriving data)
    \param commSegments data vectors for handling the parts of the communication
    \param reqs wait requests
    */
    void mergeI(const UInt                      & sid,
                const UInt                      & rid,
                const Teuchos::RCP<const PVECT> & sendVect,
                      Teuchos::RCP<PVECT>       & rVect,
                      PVPS                      & pvps,
                const UInt                      & channel = 2) const;
       
    /*! Merging the received data to \c rVect
    \param sid sending process
    \param rid receiving process
    \param sMap the data to send
    \param rVect the data received (rVect is merged with the arriving data)
    \param commSegments data vectors for handling the parts of the communication
    \param reqs wait requests
    */
    void mergeO(const UInt                      & sid,
                const UInt                      & rid,
                const Teuchos::RCP<const PVECT> & sendVect,
                      Teuchos::RCP<PVECT>       & rVect,
                      PVPS                      & pvps,
                const UInt                      & channel = 2) const;
    //@}


    /*! @name Send and Recv - recursive-comm */ //@{
  public:
    void sendRR(const UInt                & sid,
                const UInt                & rid,
                const Teuchos::RCP<PVECT> & sendVect,
                      PVRS                & pvrs,
                      sVect<request>      & reqsRR,
                const UInt                & channel = 2) const;
       
    void sendRR(const UInt           & sid,
                const UInt           & rid,
                const PVECT          & sendVect,
                      PVRS           & pvrs,
                      sVect<request> & reqsRR,
                const UInt           & channel = 2) const;

    void recvRR(const UInt                & sid,
                const UInt                & rid,
                      Teuchos::RCP<PVECT> & rVect,
                      PVRS                & pvrs,
                const UInt                & channel = 2) const;

    void recvRR(const UInt & sid,
                const UInt & rid,
                     PVECT & rVect,
                      PVRS & pvrs,
                const UInt & channel = 2) const;
       
    void sendRI(const UInt                & sid,
                const UInt                & rid,
                const Teuchos::RCP<PVECT> & sendVect,
                      PVRS                & pvrs,
                const UInt                & channel = 2) const;
       
    void sendRI(const UInt  & sid,
                const UInt  & rid,
                const PVECT & sendVect,
                      PVRS  & pvrs,
                const UInt  & channel = 2) const;

    void sendRO(const UInt                & sid,
                const UInt                & rid,
                const Teuchos::RCP<PVECT> & sendVect,
                      PVRS                & pvrs,
                const UInt                & channel = 2) const;
       
    void sendRO(const UInt  & sid,
                const UInt  & rid,
                const PVECT & sendVect,
                      PVRS  & pvrs,
                const UInt  & channel = 2) const;
       
    void recvRI(const UInt                & sid,
                const UInt                & rid,
                      Teuchos::RCP<PVECT> & rVect,
                      PVRS                & pvrs,
                const UInt                & channel = 2) const;

    void recvRI(const UInt & sid,
                const UInt & rid,
                     PVECT & rVect,
                      PVRS & pvrs,
                const UInt & channel = 2) const;

    void recvRO(const UInt                & sid,
                const UInt                & rid,
                      Teuchos::RCP<PVECT> & rVect,
                      PVRS                & pvrs,
                const UInt                & channel = 2) const;

    void recvRO(const UInt & sid,
                const UInt & rid,
                     PVECT & rVect,
                      PVRS & pvrs,
                const UInt & channel = 2) const;
    //@}

       
    /*! @name Mapped communications - non-pending-comm */ //@{
  public:
    /*! Sending and Receiving
    \param mapSend map for sending items
    \param sendVect the data to send 
    \param rVect the data received (rVect is cleared)
    */
    void sendRecv(const Teuchos::RCP<SENDRECV>    & mapSend,
                  const Teuchos::RCP<const PVECT> & sendVect,
                        Teuchos::RCP<PVECT>       & rVect) const;
    
    /*! Sending and Receiving
    \param mapSend map for sending items
    \param sendVect the data to send 
    \param rVect the data received (rVect is cleared)
    */
    void sendRecv(const SENDRECV & mapSend,
                  const PVECT    & sendVect,
                        PVECT    & rVect) const;
    
    /*! Sending and Receiving
    \param mapSend map for sending items
    \param mapRecv map for receiving items
    \param sendVect the data to send 
    \param rVect the data received (rVect is cleared)
    */
    void sendRecv(const Teuchos::RCP<const SENDRECV> & mapSend,
                  const Teuchos::RCP<const SENDRECV> & mapRecv,
                  const Teuchos::RCP<const PVECT>    & sendVect,
                        Teuchos::RCP<PVECT>          & rVect) const;
    
    /*! Sending and Receiving
    \param mapSend map for sending items
    \param mapRecv map for receiving items
    \param sendVect the data to send 
    \param rVect the data received (rVect is cleared)
    */
    void sendRecv(const SENDRECV & mapSend,
                  const SENDRECV & mapRecv,
                  const PVECT    & sendVect,
                        PVECT    & rVect) const;
    //@}
    
 
    /*! @name Mapped communications - pending-comm */ //@{
  public:
    /*! Sending and Receiving
    \param mapSend map for sending items
    \param sendVect the data to send 
    \param rVect the data received (rVect is cleared)
    \param sendPvps send communication data
    \param recvPvps recv communication data
    */
    void sendRecvI(const Teuchos::RCP<SENDRECV>    & mapSend,
                   const Teuchos::RCP<const PVECT> & sendVect,
                         Teuchos::RCP<PVECT>       & rVect,
                         sVect<PVECT>              & commSegments,
                         sVect<PVPS>               & sendPvps,
                         sVect<PVPS>               & recvPvps,
                   const UInt                      & channel = 2) const;
 
    /*! Sending and Receiving
    \param mapSend map for sending items
    \param sendVect the data to send 
    \param rVect the data received (rVect is cleared)
    \param sendPvps send communication data
    \param recvPvps recv communication data
    */
    void sendRecvO(const Teuchos::RCP<SENDRECV>    & mapSend,
                   const Teuchos::RCP<const PVECT> & sendVect,
                         Teuchos::RCP<PVECT>       & rVect,
                         sVect<PVECT>              & commSegments,
                         sVect<PVPS>               & sendPvps,
                         sVect<PVPS>               & recvPvps,
                   const UInt                      & channel = 2) const;
    
    /*! Sending and Receiving
    \param mapSend map for sending items
    \param sendVect the data to send 
    \param rVect the data received (rVect is cleared)
    \param sendPvps send communication data
    \param recvPvps recv communication data
    */
    void sendRecvI(const SENDRECV     & mapSend,
                   const PVECT        & sendVect,
                         PVECT        & rVect,
                         sVect<PVECT> & commSegments,
                         sVect<PVPS>  & sendPvps,
                         sVect<PVPS>  & recvPvps,
                   const UInt         & channel = 2) const;

    /*! Sending and Receiving
    \param mapSend map for sending items
    \param sendVect the data to send 
    \param rVect the data received (rVect is cleared)
    \param sendPvps send communication data
    \param recvPvps recv communication data
    */
    void sendRecvO(const SENDRECV     & mapSend,
                   const PVECT        & sendVect,
                         PVECT        & rVect,
                         sVect<PVECT> & commSegments,
                         sVect<PVPS>  & sendPvps,
                         sVect<PVPS>  & recvPvps,
                   const UInt         & channel = 2) const;

    /*! Sending and Receiving
    \param mapSend map for sending items
    \param mapRecv map for receiving items
    \param sendVect the data to send 
    \param rVect the data received (rVect is cleared)
    \param commSegments data vectors for handling the parts of the communication
    \param reqs wait requests
    */
    void sendRecvI(const Teuchos::RCP<const SENDRECV> & mapSend,
                   const Teuchos::RCP<const SENDRECV> & mapRecv,
                   const Teuchos::RCP<const PVECT>    & sendVect,
                         Teuchos::RCP<PVECT>          & rVect,
                         sVect<PVECT>                 & commSegments,
                         sVect<PVPS>                  & sendPvps,
                         sVect<PVPS>                  & recvPvps,
                   const UInt                         & channel = 2) const;

    /*! Sending and Receiving
    \param mapSend map for sending items
    \param mapRecv map for receiving items
    \param sendVect the data to send 
    \param rVect the data received (rVect is cleared)
    \param commSegments data vectors for handling the parts of the communication
    \param reqs wait requests
    */
    void sendRecvO(const Teuchos::RCP<const SENDRECV> & mapSend,
                   const Teuchos::RCP<const SENDRECV> & mapRecv,
                   const Teuchos::RCP<const PVECT>    & sendVect,
                         Teuchos::RCP<PVECT>          & rVect,
                         sVect<PVECT>                 & commSegments,
                         sVect<PVPS>                  & sendPvps,
                         sVect<PVPS>                  & recvPvps,
                   const UInt                         & channel = 2) const;
    
    /*! Sending and Receiving
    \param mapSend map for sending items
    \param mapRecv map for receiving items
    \param sendVect the data to send 
    \param rVect the data received (rVect is cleared)
    \param commSegments data vectors for handling the parts of the communication
    \param reqs wait requests
    */
    void sendRecvI(const SENDRECV     & mapSend,
                   const SENDRECV     & mapRecv,
                   const PVECT        & sendVect,
                         PVECT        & rVect,
                         sVect<PVECT> & commSegments,
                         sVect<PVPS>  & sendPvps,
                         sVect<PVPS>  & recvPvps,
                   const UInt         & channel = 2) const;

    /*! Sending and Receiving
    \param mapSend map for sending items
    \param mapRecv map for receiving items
    \param sendVect the data to send 
    \param rVect the data received (rVect is cleared)
    \param commSegments data vectors for handling the parts of the communication
    \param reqs wait requests
    */
    void sendRecvO(const SENDRECV     & mapSend,
                   const SENDRECV     & mapRecv,
                   const PVECT        & sendVect,
                         PVECT        & rVect,
                         sVect<PVECT> & commSegments,
                         sVect<PVPS>  & sendPvps,
                         sVect<PVPS>  & recvPvps,
                   const UInt         & channel = 2) const;
    //@}

   
    /*! @name Mapped communications - recursive-comm */ //@{
  public:
    void sendRecvRR(const Teuchos::RCP<SENDRECV>    & mapSend,
                    const Teuchos::RCP<const PVECT> & sendVect,
                          Teuchos::RCP<PVECT>       & rVect,
                          sVect<PVECT>              & sendSegments,
                          sVect<PVECT>              & recvSegments,
                          sVect<PVRS>               & sendPvrs,
                          sVect<PVRS>               & recvPvrs,
                    const UInt                      & channel = 2) const;
    
    void sendRecvRR(const SENDRECV     & mapSend,
                    const PVECT        & sendVect,
                          PVECT        & rVect,
                          sVect<PVECT> & sendSegments,
                          sVect<PVECT> & recvSegments,
                          sVect<PVRS>  & sendPvrs,
                          sVect<PVRS>  & recvPvrs,
                    const UInt         & channel = 2) const;
 
    void sendRecvRI(const Teuchos::RCP<SENDRECV>    & mapSend,
                    const Teuchos::RCP<const PVECT> & sendVect,
                          Teuchos::RCP<PVECT>       & rVect,
                          sVect<PVECT>              & sendSegments,
                          sVect<PVECT>              & recvSegments,
                          sVect<PVRS>               & sendPvrs,
                          sVect<PVRS>               & recvPvrs,
                    const UInt                      & channel = 2) const;
    
    void sendRecvRI(const SENDRECV     & mapSend,
                    const PVECT        & sendVect,
                          PVECT        & rVect,
                          sVect<PVECT> & sendSegments,
                          sVect<PVECT> & recvSegments,
                          sVect<PVRS>  & sendPvrs,
                          sVect<PVRS>  & recvPvrs,
                    const UInt         & channel = 2) const;
    
    void sendRecvRO(const Teuchos::RCP<SENDRECV>    & mapSend,
                    const Teuchos::RCP<const PVECT> & sendVect,
                          Teuchos::RCP<PVECT>       & rVect,
                          sVect<PVECT>              & sendSegments,
                          sVect<PVECT>              & recvSegments,
                          sVect<PVRS>               & sendPvrs,
                          sVect<PVRS>               & recvPvrs,
                    const UInt                      & channel = 2) const;
    
    void sendRecvRO(const SENDRECV     & mapSend,
                    const PVECT        & sendVect,
                          PVECT        & rVect,
                          sVect<PVECT> & sendSegments,
                          sVect<PVECT> & recvSegments,
                          sVect<PVRS>  & sendPvrs,
                          sVect<PVRS>  & recvPvrs,
                    const UInt         & channel = 2) const;

    void sendRecvRR(const Teuchos::RCP<SENDRECV>    & mapSend,
                    const Teuchos::RCP<SENDRECV>    & mapRecv,
                    const Teuchos::RCP<const PVECT> & sendVect,
                          Teuchos::RCP<PVECT>       & rVect,
                          sVect<PVECT>              & sendSegments,
                          sVect<PVECT>              & recvSegments,
                          sVect<PVRS>               & sendPvrs,
                          sVect<PVRS>               & recvPvrs,
                    const UInt                      & channel = 2) const;
    
    void sendRecvRR(const SENDRECV     & mapSend,
                    const SENDRECV     & mapRecv,
                    const PVECT        & sendVect,
                          PVECT        & rVect,
                          sVect<PVECT> & sendSegments,
                          sVect<PVECT> & recvSegments,
                          sVect<PVRS>  & sendPvrs,
                          sVect<PVRS>  & recvPvrs,
                    const UInt         & channel = 2) const;
 
    void sendRecvRI(const Teuchos::RCP<SENDRECV>    & mapSend,
                    const Teuchos::RCP<SENDRECV>    & mapRecv,
                    const Teuchos::RCP<const PVECT> & sendVect,
                          Teuchos::RCP<PVECT>       & rVect,
                          sVect<PVECT>              & sendSegments,
                          sVect<PVECT>              & recvSegments,
                          sVect<PVRS>               & sendPvrs,
                          sVect<PVRS>               & recvPvrs,
                    const UInt                      & channel = 2) const;
    
    void sendRecvRI(const SENDRECV     & mapSend,
                    const SENDRECV     & mapRecv,
                    const PVECT        & sendVect,
                          PVECT        & rVect,
                          sVect<PVECT> & sendSegments,
                          sVect<PVECT> & recvSegments,
                          sVect<PVRS>  & sendPvrs,
                          sVect<PVRS>  & recvPvrs,
                    const UInt         & channel = 2) const;
    
    void sendRecvRO(const Teuchos::RCP<SENDRECV>    & mapSend,
                    const Teuchos::RCP<SENDRECV>    & mapRecv,
                    const Teuchos::RCP<const PVECT> & sendVect,
                          Teuchos::RCP<PVECT>       & rVect,
                          sVect<PVECT>              & sendSegments,
                          sVect<PVECT>              & recvSegments,
                          sVect<PVRS>               & sendPvrs,
                          sVect<PVRS>               & recvPvrs,
                    const UInt                      & channel = 2) const;
    
    void sendRecvRO(const SENDRECV     & mapSend,
                    const SENDRECV     & mapRecv,
                    const PVECT        & sendVect,
                          PVECT        & rVect,
                          sVect<PVECT> & sendSegments,
                          sVect<PVECT> & recvSegments,
                          sVect<PVRS>  & sendPvrs,
                          sVect<PVRS>  & recvPvrs,
                    const UInt         & channel = 2) const;
    //@}

    
    /*! @name Pattern communications - non-pending-comm */ //@{
  public:
     /*! Normal data distribution (gid-based)
    \param Vect contains the data to send, then is cleared and loaded with the new data.
    The \c maxGid parameter can be used to force the vector to a different distribution
    on the communciator
    */
    void vectorNormal(Teuchos::RCP<PVECT> & Vect,
                      const UInt & MaxGid = 0);
    
     /*! Normal data distribution (gid-based)
    \param Vect contains the data to send, then is cleared and loaded with the new data.
    The \c maxGid parameter can be used to force the vector to a different distribution
    on the communciator
    */
    void vectorNormal(PVECT & Vect,
                      const UInt & MaxGid = 0);
    
    /*! Pid data distribution (pid-based)
    \param Vect contains the data to send, then is cleared and loaded with the new data
    */
    void vectorPid(Teuchos::RCP<PVECT> & Vect);
    
    /*! Pid data distribution (pid-based)
    \param Vect contains the data to send, then is cleared and loaded with the new data
    */
    void vectorPid(PVECT & Vect);
    
    /*! Data distribution (data-based)
    \param Vect contains the data to send, then is cleared and loaded with the new data
    */
    void vectorData(Teuchos::RCP<PVECT> & Vect);
    
    /*! Data distribution (data-based)
    \param Vect contains the data to send, then is cleared and loaded with the new data
    */
    void vectorData(PVECT & Vect);
    
    /*! Data distribution (data-based)
    \param Vect contains the data to send, then is cleared and loaded with the new data
    \param dataMin the minimum datum 
    \param dataMax the maximum datum
    */
    void vectorData(Teuchos::RCP<PVECT> & Vect,
                    const DATA & dataMin,
                    const DATA & dataMax);
    
    /*! Data distribution (data-based)
    \param Vect contains the data to send, then is cleared and loaded with the new data
    \param dataMin the minimum datum 
    \param dataMax the maximum datum
    */
    void vectorData(PVECT & Vect,
                    const DATA & dataMin,
                    const DATA & dataMax);
    //@}
    
};



//_______________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
pVectComm<DATA,MAPITEM>::
pVectComm() : commDevLoaded(false)
{
  maxSize = mpiMaxSize / (sizeof(MAPITEM) + sizeof(DATA));
}

template<typename DATA, typename MAPITEM>
pVectComm<DATA,MAPITEM>::
pVectComm(const Teuchos::RCP<const communicator> & CommDev) : commDevLoaded(true), commDev(CommDev)
{
   maxSize = mpiMaxSize / (sizeof(MAPITEM) + sizeof(DATA));
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}



//_______________________________________________________________________________________________________
// SEND AND RECV- NON-PENDING-COMM
//-------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecv(const UInt & sid,
         const UInt & rid,
         const Teuchos::RCP<const PVECT> & sendVect,
               Teuchos::RCP<PVECT>       & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRecv(sid,rid,*sendVect,*rVect);
}


#ifdef NOCOMMBUFFER

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecv(const UInt & sid,
         const UInt & rid,
         const PVECT & sendVect,
               PVECT & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {    
    //Inizio trasmissione
    sVect<request> reqs(1);
    reqs(1) = commDev->isend(rid,1,sendVect);
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Inizio trasmissione
    rVect.clear();
    
    sVect<request> reqs(1);
    reqs(1) = commDev->irecv(sid,1,rVect);
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
}

#else

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecv(const UInt & sid,
         const UInt & rid,
         const PVECT & sendVect,
               PVECT & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
     //Segmentazione dati
    pVectManip<DATA,MAPITEM> manipulator;
    sVect<PVECT>               segmented;
    
    manipulator.segmentationSimple(segmented,maxSize,sendVect);
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
    
    assert(n <= mpiLayer);
    sVect<PVECT> segmented(n);
    
    //Inizio trasmissione
    sVect<request> reqs(n);
    
    for(UInt i=1; i<= n; ++i)
    { 
      reqs(i) = commDev->irecv(sid,i,segmented(i));
    }
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  
    //Ricostruzione dati
    rVect.clear();
    pVectManip<DATA,MAPITEM> manipulator;
    manipulator.mergeSimple(segmented,rVect);
  }
}

#endif


//-------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
send(const UInt & sid,
     const UInt & rid,
     const Teuchos::RCP<const PVECT> & sendVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  send(sid,rid,*sendVect);
}


#ifdef NOCOMMBUFFER

template<typename DATA, typename MAPITEM>    
void
pVectComm<DATA,MAPITEM>::
send(const UInt & sid,
     const UInt & rid,
     const PVECT & sendVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {    
    //Inizio trasmissione
    sVect<request> reqs(1);
    reqs(1) = commDev->isend(rid,1,sendVect);
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
}

#else

template<typename DATA, typename MAPITEM>    
void
pVectComm<DATA,MAPITEM>::
send(const UInt & sid,
     const UInt & rid,
     const PVECT & sendVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    //Segmentazione dati
    pVectManip<DATA,MAPITEM> manipulator;
    sVect<PVECT>               segmented;
    
    manipulator.segmentationSimple(segmented,maxSize,sendVect);
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


//-------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
isend(const UInt & sid,
      const UInt & rid,
      const Teuchos::RCP<const PVECT> & sendVect,
      sVect<request> & reqs) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  isend(sid,rid,*sendVect,reqs);
}


#ifdef NOCOMMBUFFER

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
isend(const UInt & sid,
      const UInt & rid,
      const PVECT & sendVect,
      sVect<request> & reqs) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {    
    //Inizio trasmissione   
    reqs.reserve(1);
    reqs.push_back( commDev->isend(rid,1,sendVect) );
  }
}

#else

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
isend(const UInt & sid,
      const UInt & rid,
      const PVECT & sendVect,
      sVect<request> & reqs) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    pVectManip<DATA,MAPITEM> manipulator;
    sVect<PVECT>             segmented;
      
    manipulator.segmentationSimple(segmented,maxSize,sendVect);
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


//-------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recv(const UInt & sid,
     const UInt & rid,
     Teuchos::RCP<PVECT> & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  recv(sid,rid,*rVect);
}

#ifdef NOCOMMBUFFER

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recv(const UInt & sid,
     const UInt & rid,
     PVECT & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Inizio trasmissione
    rVect.clear();
    
    sVect<request> reqs(1);
    reqs(1) = commDev->irecv(sid,1,rVect);
   
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
}

#else

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recv(const UInt & sid,
     const UInt & rid,
     PVECT & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Ricezione dimensioni ed allocazione mappa ricezione
    UInt n;
    commDev->recv(sid,1,n);
    
    sVect<PVECT> segmented(n);
    assert(n <= mpiLayer);
    
    //Inizio trasmissione
    sVect<request> reqs(n);
    
    for(UInt i=1; i<= n; ++i)
    { 
      reqs(i) = commDev->irecv(sid,i,segmented(i));
    }
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  
    //Ricostruzione dati
    rVect.clear();
    pVectManip<DATA,MAPITEM> manipulator;
    manipulator.mergeSimple(segmented,rVect);
  }
}

#endif


//-------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
irecv(const UInt & sid,
      const UInt & rid,
      Teuchos::RCP<PVECT> & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  irecv(sid,rid,*rVect);
}


#ifdef NOCOMMBUFFER

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
irecv(const UInt & sid,
      const UInt & rid,
      PVECT & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Reciver
  if(commDev->rank() == int(rid))
  {    
    //Inizio trasmissione
    rVect.clear();
    
    sVect<request> reqs(1);
    reqs(1) = commDev->irecv(sid,1,rVect);
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
}

#else

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
irecv(const UInt & sid,
      const UInt & rid,
      PVECT & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Ricezione dimensioni ed allocazione mappa ricezione
    UInt n;
    commDev->recv(sid,1,n);
    
    sVect<PVECT> segmented(n);
    assert(n <= mpiLayer);
    
    //Inizio trasmissione
    sVect<request> reqs(n);
    
    for(UInt i=1; i<= n; ++i)
    { 
      reqs(i) = commDev->irecv(sid,i,segmented(i));
    }
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
    
    //Ricostruzione dati
    rVect.clear();
    pVectManip<DATA,MAPITEM> manipulator;
    manipulator.mergeSimple(segmented,rVect);
  }
}

#endif


//-------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
merge(const UInt & sid,
      const UInt & rid,
      const Teuchos::RCP<const PVECT> & sendVect,
            Teuchos::RCP<PVECT>       & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  merge(sid,rid,*sendVect,*rVect);
}


#ifdef NOCOMMBUFFER

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
merge(const UInt & sid,
      const UInt & rid,
      const PVECT & sendVect,
            PVECT & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    //Inizio trasmissione
    sVect<request> reqs(1);
    reqs(1) = commDev->isend(rid,1,sendVect);
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
  }
  
  //Reciver
  if(commDev->rank() == int(rid))
  {    
    sVect<PVECT>   segmented(1);
    sVect<request> reqs(1);
    
    reqs(1) = commDev->irecv(sid,1,segmented(1));
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
    
    //Ricostruzione dati
    pVectManip<DATA,MAPITEM> manipulator;
    manipulator.mergeSimple(segmented,rVect);
  }
}

#else

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
merge(const UInt & sid,
      const UInt & rid,
      const PVECT & sendVect,
            PVECT & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    //Segmentazione dati
    pVectManip<DATA,MAPITEM> manipulator;
    sVect<PVECT>             segmented;
    
    manipulator.segmentationSimple(segmented,maxSize,sendVect);
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
    
    sVect<PVECT> segmented(n);
    assert(n <= mpiLayer);
    
    //Inizio trasmissione
    sVect<request> reqs(n);
    
    for(UInt i=1; i<= n; ++i)
    { 
      reqs(i) = commDev->irecv(sid,i,segmented(i));
    }
    
    //Waiting
    wait_all(reqs.begin(),reqs.end());
    
    //Ricostruzione dati
    pVectManip<DATA,MAPITEM> manipulator;
    manipulator.mergeSimple(segmented,rVect);
  }
}

#endif



//_______________________________________________________________________________________________________
// SEND AND RECV- PENDING-COMM
//-------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendI(const UInt                & sid,
      const UInt                & rid,
      const Teuchos::RCP<PVECT> & sendVect,
            PVPS                & pvps,
      const UInt                & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendI(sid,
        rid,
       *sendVect,
        pvps,
        channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendO(const UInt                & sid,
      const UInt                & rid,
      const Teuchos::RCP<PVECT> & sendVect,
            PVPS                & pvps,
      const UInt                & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendO(sid,
        rid,
       *sendVect,
        pvps,
        channel);
}


#ifdef NOCOMMBUFFER

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendI(const UInt  & sid,
      const UInt  & rid,
      const PVECT & sendVect,
            PVPS  & pvps,
      const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {    
    //Inizio trasmissione   
    pvps.reqs.reserve(1);
    pvps.reqs.push_back( commDev->isend(rid, (channel-1) * mpiLayer +1, sendVect) );
  }
}

#else

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendI(const UInt  & sid,
      const UInt  & rid,
      const PVECT & sendVect,
            PVPS  & pvps,
      const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    pVectManip<DATA,MAPITEM> manipulator;
    sVect<PVECT>             segmented;
      
    manipulator.segmentationSimple(segmented,maxSize,sendVect);
    UInt n = segmented.size();
    assert(n <= mpiLayer);
    
    //Trasmissione dimensioni
    commDev->send(rid, (channel-1) * mpiLayer +1, n);
    
    //Inizio trasmissione   
    pvps.reqs.reserve(n);
    
    for(UInt i=1; i<= n; ++i)
    { pvps.reqs.push_back( commDev->isend(rid, (channel-1) * mpiLayer +i, segmented(i)) ); }
  }
}

#endif

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendO(const UInt  & sid,
      const UInt  & rid,
      const PVECT & sendVect,
            PVPS  & pvps,
      const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Barrier
  if(commDev->rank() == int(sid))
  { wait_all(pvps.reqs.begin(),pvps.reqs.end()); }
}


//-------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recvI(const UInt                & sid,
      const UInt                & rid,
            Teuchos::RCP<PVECT> & rVect,
            PVPS                & pvps,
      const UInt                & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  recvI(sid,
        rid,
       *rVect,
        pvps,
        channel);
}
     
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recvO(const UInt                & sid,
      const UInt                & rid,
            Teuchos::RCP<PVECT> & rVect,
            PVPS                & pvps,
      const UInt                & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  recvO(sid,
        rid,
       *rVect,
        pvps,
        channel);
}     


#ifdef NOCOMMBUFFER

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recvI(const UInt  & sid,
      const UInt  & rid,
            PVECT & rVect,
            PVPS  & pvps,
      const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Reciver
  if(commDev->rank() == int(rid))
  {    
    //Inizio trasmissione
    rVect.clear();
    
    pvps.reqs.resize(1);
    pvps.reqs(1) = commDev->irecv(sid, (channel-1) * mpiLayer +1, rVect);
  }
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recvO(const UInt  & sid,
      const UInt  & rid,
           PVECT  & rVect,
            PVPS  & pvps,
      const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Waiting
  if(commDev->rank() == int(rid))
  { wait_all(pvps.reqs.begin(),pvps.reqs.end()); }
}

#else

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recvI(const UInt  & sid,
      const UInt  & rid,
            PVECT & rVect,
            PVPS  & pvps,
      const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Ricezione dimensioni ed allocazione mappa ricezione
    UInt n;
    commDev->recv(sid, (channel-1) * mpiLayer +1, n);
    
    pvps.bufVects.resize(n);
    
    //Inizio trasmissione
    pvps.reqs.resize(n);
    
    for(UInt i=1; i<= n; ++i)
    { pvps.reqs(i) = commDev->irecv(sid, (channel-1) * mpiLayer +i, segmented(i)); }
  }
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recvO(const UInt  & sid,
      const UInt  & rid,
            PVECT & rVect,
            PVPS  & pvps,
      const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Reciver
  if(commDev->rank() == int(rid))
  {   
    //Waiting
    wait_all(pvps.reqs.begin(),pvps.reqs.end());
    
    //Ricostruzione dati
    rVect.clear();
    pVectManip<DATA,MAPITEM> manipulator;
    manipulator.mergeSimple(pvps.bufVects,rVect);
  }
}

#endif


//-------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvI(const UInt                      & sid,
          const UInt                      & rid,
          const Teuchos::RCP<const PVECT> & sendVect,
                Teuchos::RCP<PVECT>       & rVect,
                PVPS                      & pvps,
          const UInt                      & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Send
  sendRecvI(sid,
            rid,
           *sendVect,
           *rVect,
            pvps,
            channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvO(const UInt                      & sid,
          const UInt                      & rid,
          const Teuchos::RCP<const PVECT> & sendVect,
                Teuchos::RCP<PVECT>       & rVect,
                PVPS                      & pvps,
          const UInt                      & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Send
  sendRecvO(sid,
            rid,
           *sendVect,
           *rVect,
            pvps,
            channel);
}


#ifdef NOCOMMBUFFER

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvI(const UInt  & sid,
          const UInt  & rid,
          const PVECT & sendVect,
                PVECT & rVect,
                PVPS  & pvps,
          const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {    
    //Inizio trasmissione
    pvps.reqs.resize(1);
    pvps.reqs(1) = commDev->isend(rid, (channel-1) * mpiLayer +1, sendVect);
  }
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Inizio trasmissione
    rVect.clear();
    
    pvps.reqs.resize(1);
    pvps.reqs(1) = commDev->irecv(sid, (channel-1) * mpiLayer +1, rVect);
  }
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvO(const UInt  & sid,
          const UInt  & rid,
          const PVECT & sendVect,
                PVECT & rVect,
                PVPS  & pvps,
          const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {    
    //Waiting
    wait_all(pvps.reqs.begin(),pvps.reqs.end());
  }
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Waiting
    wait_all(pvps.reqs.begin(),pvps.reqs.end());
  }
}
  
#else

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvI(const UInt  & sid,
          const UInt  & rid,
          const PVECT & sendVect,
                PVECT & rVect,
                PVPS  & pvps,
          const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
     //Segmentazione dati
    pVectManip<DATA,MAPITEM> manipulator;
    
    manipulator.segmentationSimple(pvps.bufVects,maxSize,sendVect);
    UInt n = pvps.bufVects.size();
    assert(n <= mpiLayer);
    
    //Trasmissione dimensioni
    commDev->send(rid, (channel-1) * mpiLayer +1, n);
    
    //Inizio trasmissione
    pvps.reqs.resize(n);
    
    for(UInt i=1; i<= n; ++i)
    { pvps.reqs(i) = commDev->isend(rid, (channel-1) * mpiLayer +i, pvps.bufVects(i)); }
  }
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Ricezione dimensioni ed allocazione mappa ricezione
    UInt n;
    commDev->recv(sid, (channel-1) * mpiLayer +1, n);
    
    pvps.bufVects.resize(n);
    
    //Inizio trasmissione
    pvps.reqs.resize(n);
    
    for(UInt i=1; i<= n; ++i)
    { pvps.reqs(i) = commDev->irecv(sid, (channel-1) * mpiLayer +i, pvps.bufVects(i)); }

  }
}
 
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvO(const UInt  & sid,
          const UInt  & rid,
          const PVECT & sendVect,
                PVECT & rVect,
                PVPS  & pvps,
          const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    //Waiting
    wait_all(pvps.reqs.begin(),pvps.reqs.end());
  }
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Waiting
    wait_all(pvps.reqs.begin(),pvps.reqs.end());
  
    //Ricostruzione dati
    rVect.clear();
    pVectManip<DATA,MAPITEM> manipulator;
    manipulator.mergeSimple(pvps.bufVects,rVect);
  }
}

#endif


//-------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
mergeI(const UInt                      & sid,
       const UInt                      & rid,
       const Teuchos::RCP<const PVECT> & sendVect,
             Teuchos::RCP<PVECT>       & rVect,
             PVPS                      & pvps,
       const UInt                      & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Send
  mergeI(sid,
         rid,
        *sendVect,
        *rVect,
         pvps,
         channel);
}
       
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
mergeO(const UInt                      & sid,
       const UInt                      & rid,
       const Teuchos::RCP<const PVECT> & sendVect,
             Teuchos::RCP<PVECT>       & rVect,
             PVPS                      & pvps,
       const UInt                      & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Send
  mergeO(sid,
         rid,
        *sendVect,
        *rVect,
         pvps,
         channel);
}


#ifdef NOCOMMBUFFER

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
mergeI(const UInt  & sid,
       const UInt  & rid,
       const PVECT & sendVect,
             PVECT & rVect,
             PVPS  & pvps,
       const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    //Inizio trasmissione
    pvps.reqs.resize(1);
    pvps.reqs(1) = commDev->isend(rid, (channel-1) * mpiLayer +1, sendVect);
  }
  
  //Reciver
  if(commDev->rank() == int(rid))
  {    
    pvps.bufVects.resize(1);
    pvps.reqs.resize(1);
    
    pvps.reqs(1) = commDev->irecv(sid, (channel-1) * mpiLayer +1, pvps.bufVects(1));
  }
}
    
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
mergeO(const UInt  & sid,
       const UInt  & rid,
       const PVECT & sendVect,
             PVECT & rVect,
             PVPS  & pvps,
       const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    //Waiting
    wait_all(pvps.reqs.begin(),pvps.reqs.end());
  }
  
  //Reciver
  if(commDev->rank() == int(rid))
  {    
    //Waiting
    wait_all(pvps.reqs.begin(),pvps.reqs.end());
    
    //Ricostruzione dati
    pVectManip<DATA,MAPITEM> manipulator;
    manipulator.mergeSimple(pvps.bufVects,rVect);
  }
}

#else

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
mergeI(const UInt  & sid,
       const UInt  & rid,
       const PVECT & sendVect,
             PVECT & rVect,
             PVPS  & pvps,
       const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    //Segmentazione dati
    pVectManip<DATA,MAPITEM> manipulator;
    
    manipulator.segmentationSimple(pvps.bufVects,maxSize,sendVect);
    UInt n = pvps.bufVects.size();
    assert(n <= mpiLayer);
    
    //Trasmissione dimensioni
    commDev->send(rid, (channel-1) * mpiLayer +1, n);
    
    //Inizio trasmissione
    pvps.reqs.resize(n);
    
    for(UInt i=1; i<= n; ++i)
    { pvps.reqs(i) = commDev->isend(rid, (channel-1) * mpiLayer +i, pvps.bufVects(i)); }
  }
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Ricezione dimensioni ed allocazione mappa ricezione
    UInt n;
    commDev->recv(sid, (channel-1) * mpiLayer +1, n);
    pvps.bufVects.resize(n);
    
    //Inizio trasmissione
    pvps.reqs.resize(n);
    
    for(UInt i=1; i<= n; ++i)
    { pvps.reqs(i) = commDev->irecv(sid, (channel-1) * mpiLayer +i, pvps.bufVects(i)); }
  }
}
    
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
mergeO(const UInt  & sid,
       const UInt  & rid,
       const PVECT & sendVect,
             PVECT & rVect,
             PVPS  & pvps,
       const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  {
    //Waiting
    wait_all(pvps.reqs.begin(),pvps.reqs.end());
  }
  
  //Reciver
  if(commDev->rank() == int(rid))
  {
    //Waiting
    wait_all(pvps.reqs.begin(),pvps.reqs.end());
    
    //Ricostruzione dati
    pVectManip<DATA,MAPITEM> manipulator;
    manipulator.mergeSimple(pvps.bufVects,rVect);
  }
}

#endif



//_______________________________________________________________________________________________________
// SEND AND RECV- RECURSIVE-COMM
//-------------------------------------------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRR(const UInt                & sid,
       const UInt                & rid,
       const Teuchos::RCP<PVECT> & sendVect,
             PVRS                & pvrs,
             sVect<request>      & reqsRR,
       const UInt                & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRR(sid,
         rid,
        *sendVect,
         pvrs,
         channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRR(const UInt           & sid,
       const UInt           & rid,
       const PVECT          & sendVect,
             PVRS           & pvrs,
             sVect<request> & reqsRR,
       const UInt           & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  if(commDev->rank() == int(sid))
  {
    //Send and sizing of the vector
    reqsRR.push_back(commDev->isend(rid, (channel-1) * mpiLayer +1, sendVect));

    //Send of the skeleton
    reqsRR.push_back(commDev->isend(rid, (channel-1) * mpiLayer +2, boost::mpi::skeleton(*sendVect.data)));
    
    //Alloc the contenent
    pvrs.bufCont = boost::mpi::get_content(*sendVect.data);
    pvrs.reqs.resize(1);
  }
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recvRR(const UInt                & sid,
       const UInt                & rid,
             Teuchos::RCP<PVECT> & rVect,
             PVRS                & pvrs,
       const UInt                & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  recvRR(sid,
         rid,
        *rVect,
         pvrs,
         channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recvRR(const UInt & sid,
       const UInt & rid,
            PVECT & rVect,
             PVRS & pvrs,
       const UInt & channel) const
{
  //Assert
  assert(commDevLoaded);

  //Receiving and sizing of the vector
  if(commDev->rank() == int(rid))
  {
    //Recv and sizing of the vector
    rVect.clear();
    
    sVect<request> reqsV(1);
    reqsV(1) = commDev->irecv(sid, (channel-1) * mpiLayer +1, rVect);
    wait_all(reqsV.begin(),reqsV.end());

    //Receiving of the skeleton
    sVect<request> reqs(1);
    reqs(1) = commDev->irecv(sid, (channel-1) * mpiLayer +2, boost::mpi::skeleton(*rVect.data));
    wait_all(reqs.begin(),reqs.end());
    
    //Alloc the contenent
    pvrs.bufCont = boost::mpi::get_content(*rVect.data);
    pvrs.reqs.resize(1);
  }
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRI(const UInt                & sid,
       const UInt                & rid,
       const Teuchos::RCP<PVECT> & sendVect,
             PVRS                & pvrs,
       const UInt                & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRI(sid,
         rid,
        *sendVect,
         pvrs,
         channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRI(const UInt  & sid,
       const UInt  & rid,
       const PVECT & sendVect,
             PVRS  & pvrs,
       const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Sender
  if(commDev->rank() == int(sid))
  { pvrs.reqs(1) = commDev->isend(rid, (channel-1) * mpiLayer +1, pvrs.bufCont); }
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recvRI(const UInt                & sid,
       const UInt                & rid,
             Teuchos::RCP<PVECT> & rVect,
             PVRS                & pvrs,
       const UInt                & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  recvRI(sid,
         rid,
        *rVect,
         pvrs,
         channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recvRI(const UInt & sid,
       const UInt & rid,
            PVECT & rVect,
             PVRS & pvrs,
       const UInt & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Reciver
  if(commDev->rank() == int(rid))
  { pvrs.reqs(1) = commDev->irecv(sid, (channel-1) * mpiLayer +1, pvrs.bufCont); }
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRO(const UInt                & sid,
       const UInt                & rid,
       const Teuchos::RCP<PVECT> & sendVect,
             PVRS                & pvrs,
       const UInt                & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRO(sid,
         rid,
        *sendVect,
         pvrs,
         channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRO(const UInt  & sid,
       const UInt  & rid,
       const PVECT & sendVect,
             PVRS  & pvrs,
       const UInt  & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Barrier
  if(commDev->rank() == int(sid))
  { wait_all(pvrs.reqs.begin(),pvrs.reqs.end()); }
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recvRO(const UInt                & sid,
       const UInt                & rid,
             Teuchos::RCP<PVECT> & rVect,
             PVRS                & pvrs,
       const UInt                & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  recvRO(sid,
         rid,
        *rVect,
         pvrs,
         channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
recvRO(const UInt & sid,
       const UInt & rid,
            PVECT & rVect,
             PVRS & pvrs,
       const UInt & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Barrier
  if(commDev->rank() == int(rid))
  { wait_all(pvrs.reqs.begin(),pvrs.reqs.end()); }
}



//_______________________________________________________________________________________________________
// MAPPED COMMUNICATIONS - NON-PENDING-COMM
//-------------------------------------------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecv(const Teuchos::RCP<SENDRECV>    & mapSend,
         const Teuchos::RCP<const PVECT> & sendVect,
               Teuchos::RCP<PVECT>       & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRecv(mapSend,*sendVect,*rVect);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecv(const SENDRECV & mapSend,
         const PVECT    & sendVect,
               PVECT    & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Data
  UInt pid    = commDev->rank();
  UInt maxPid = commDev->size();  
  
  //Segmentazione dati
  pVectManip<DATA,MAPITEM> manipulator;
  sVect<PVECT>             segmented;
  
  manipulator.segmentationMap(segmented,mapSend,maxPid,sendVect);
   
  //Sending data
  sVect<request> reqs;
  
  for(UInt rid=0; rid < maxPid; ++rid)
  {
    if(rid != pid)
    {
      isend(pid,rid,segmented[rid],reqs);
    }
  }
  
  //Receiving data
  segmented.clear();
  segmented.resize(maxPid);
  
  for(UInt sid=0; sid < maxPid; ++sid)
  {
    if(sid != pid)
    {
      irecv(sid,pid,segmented[sid]);
    }
  }
  
  //Waiting
  wait_all(reqs.begin(),reqs.end());
  
  //Merging
  rVect.clear();
  manipulator.mergeSimple(segmented,rVect);
}


template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecv(const Teuchos::RCP<const SENDRECV> & mapSend,
         const Teuchos::RCP<const SENDRECV> & mapRecv,
         const Teuchos::RCP<const PVECT>    & sendVect,
               Teuchos::RCP<PVECT>          & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRecv(*mapSend,*mapRecv,*sendVect,*rVect);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecv(const SENDRECV & mapSend,
         const SENDRECV & mapRecv,
         const PVECT    & sendVect,
               PVECT    & rVect) const
{
  //Assert
  assert(commDevLoaded);
  
  //Data
  UInt pid    = commDev->rank();
  UInt maxPid = commDev->size();
  
  //Data count
  pMapGlobalManip<pMapItemSendRecv> mapCounter(commDev);
  sVect<UInt> numSend = mapCounter.countSend(mapSend);
  sVect<UInt> numRecv = mapCounter.countRecv(mapRecv);
  
  //Segmentazione dati
  pVectManip<DATA,MAPITEM> manipulator;
  sVect<PVECT>             segmented;
  
  manipulator.segmentationMap(segmented,mapSend,maxPid,sendVect);
   
  //Sending data
  sVect<request> reqs;
  
  for(UInt rid=0; rid < maxPid; ++rid)
  {   
    assert(segmented[rid].size() == numSend[rid]);
    
    if( (rid != pid) && (numSend[rid] != 0) )
    {
      isend(pid,rid,segmented[rid],reqs);
    }
  }
  
  //Receiving data
  segmented.clear();
  segmented.resize(maxPid);
  
  for(UInt sid=0; sid < maxPid; ++sid)
  {
    if( (sid != pid) && (numRecv[sid] != 0) )
    {
      irecv(sid,pid,segmented[sid]);
    }
  }
  
  //Waiting
  wait_all(reqs.begin(),reqs.end());
  
  //Merging
  rVect.clear();
  manipulator.mergeSimple(segmented,rVect);
}


//_______________________________________________________________________________________________________
// MAPPED COMMUNICATIONS - PENDING-COMM
//-------------------------------------------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvI(const Teuchos::RCP<SENDRECV>    & mapSend,
          const Teuchos::RCP<const PVECT> & sendVect,
                Teuchos::RCP<PVECT>       & rVect,
                sVect<PVECT>              & commSegments,
                sVect<PVPS>               & sendPvps,
                sVect<PVPS>               & recvPvps,
          const UInt                      & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRecvI(*mapSend,
            *sendVect,
            *rVect,
             commSegments,
             sendPvps,
             recvPvps,
             channel);
}
 
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvO(const Teuchos::RCP<SENDRECV>    & mapSend,
          const Teuchos::RCP<const PVECT> & sendVect,
                Teuchos::RCP<PVECT>       & rVect,
                sVect<PVECT>              & commSegments,
                sVect<PVPS>               & sendPvps,
                sVect<PVPS>               & recvPvps,
          const UInt                      & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRecvO(*mapSend,
            *sendVect,
            *rVect,
             commSegments,
             sendPvps,
             recvPvps,
             channel);
}
    
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvI(const SENDRECV     & mapSend,
          const PVECT        & sendVect,
                PVECT        & rVect,
                sVect<PVECT> & commSegments,
                sVect<PVPS>  & sendPvps,
                sVect<PVPS>  & recvPvps,
          const UInt         & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Data
  UInt pid    = commDev->rank();
  UInt maxPid = commDev->size();
  
  //Segmentazione dati
  pVectManip<DATA,MAPITEM> manipulator;
  manipulator.segmentationMap(commSegments,mapSend,maxPid,sendVect);
  
  sendPvps.resize(maxPid);
  recvPvps.resize(maxPid);
   
  //Sending data  
  for(UInt rid=0; rid < maxPid; ++rid)
  {
    if(rid != pid)
    {
      sendI(pid,
            rid,
            commSegments[rid],
            sendPvps[rid],
            channel);
    }
  }
  
  //Receiving data
  commSegments.clear();
  commSegments.resize(maxPid);
  
  for(UInt sid=0; sid < maxPid; ++sid)
  {
    if(sid != pid)
    {
      recvI(sid,
            pid,
            commSegments[sid],
            recvPvps[sid],
            channel);
    }
  }
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvO(const SENDRECV     & mapSend,
          const PVECT        & sendVect,
                PVECT        & rVect,
                sVect<PVECT> & commSegments,
                sVect<PVPS>  & sendPvps,
                sVect<PVPS>  & recvPvps,
          const UInt         & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Data
  UInt pid    = commDev->rank();
  UInt maxPid = commDev->size();
  
  //Waiting receiving data  
  for(UInt sid=0; sid < maxPid; ++sid)
  {
    if(sid != pid)
    {
      recvO(sid,
            pid,
            commSegments[sid],
            recvPvps[sid],
            channel);
    }
  }
  
  //Waiting sending data  
  for(UInt rid=0; rid < maxPid; ++rid)
  {
    if(rid != pid)
    {
      sendO(pid,
            rid,
            commSegments[rid],
            sendPvps[rid],
            channel);
    }
  }
  
  //Merging
  rVect.clear();
  
  pVectManip<DATA,MAPITEM> manipulator;
  manipulator.mergeSimple(commSegments,rVect);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvI(const Teuchos::RCP<const SENDRECV> & mapSend,
          const Teuchos::RCP<const SENDRECV> & mapRecv,
          const Teuchos::RCP<const PVECT>    & sendVect,
                Teuchos::RCP<PVECT>          & rVect,
                sVect<PVECT>                 & commSegments,
                sVect<PVPS>                  & sendPvps,
                sVect<PVPS>                  & recvPvps,
          const UInt                         & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRecvI(*mapSend,
            *mapRecv,
            *sendVect,
            *rVect,
             commSegments,
             sendPvps,
             recvPvps,
             channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvO(const Teuchos::RCP<const SENDRECV> & mapSend,
          const Teuchos::RCP<const SENDRECV> & mapRecv,
          const Teuchos::RCP<const PVECT>    & sendVect,
                Teuchos::RCP<PVECT>          & rVect,
                sVect<PVECT>                 & commSegments,
                sVect<PVPS>                  & sendPvps,
                sVect<PVPS>                  & recvPvps,
          const UInt                         & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRecvO(*mapSend,
            *mapRecv,
            *sendVect,
            *rVect,
             commSegments,
             sendPvps,
             recvPvps,
             channel);
}
    
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvI(const SENDRECV       & mapSend,
          const SENDRECV       & mapRecv,
          const PVECT          & sendVect,
                PVECT          & rVect,
                sVect<PVECT>   & commSegments,
                sVect<PVPS>    & sendPvps,
                sVect<PVPS>    & recvPvps,
          const UInt           & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Data
  UInt pid    = commDev->rank();
  UInt maxPid = commDev->size();
  
  //Data count
  pMapGlobalManip<pMapItemSendRecv> mapCounter(commDev);
  sVect<UInt> numSend = mapCounter.countSend(mapSend);
  sVect<UInt> numRecv = mapCounter.countRecv(mapRecv);
  
  //Segmentazione dati
  pVectManip<DATA,MAPITEM> manipulator;
  manipulator.segmentationMap(commSegments,mapSend,maxPid,sendVect);
  
  sendPvps.resize(maxPid);
  recvPvps.resize(maxPid);
   
  //Sending data
  for(UInt rid=0; rid < maxPid; ++rid)
  {   
    assert(commSegments[rid].size() == numSend[rid]);
    
    if( (rid != pid) && (numSend[rid] != 0) )
    {
      sendI(pid,
            rid,
            commSegments[rid],
            sendPvps[rid],
            channel);
    }
  }
  
  //Receiving data
  commSegments.clear();
  commSegments.resize(maxPid);
  
  for(UInt sid=0; sid < maxPid; ++sid)
  {
    if( (sid != pid) && (numRecv[sid] != 0) )
    {
      recvI(sid,
            pid,
            commSegments[sid],
            recvPvps[sid],
            channel);
    }
  }
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvO(const SENDRECV       & mapSend,
          const SENDRECV       & mapRecv,
          const PVECT          & sendVect,
                PVECT          & rVect,
                sVect<PVECT>   & commSegments,
                sVect<PVPS>    & sendPvps,
                sVect<PVPS>    & recvPvps,
          const UInt           & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Data
  UInt pid    = commDev->rank();
  UInt maxPid = commDev->size();
  
  //Data count
  pMapGlobalManip<pMapItemSendRecv> mapCounter(commDev);
  sVect<UInt> numSend = mapCounter.countSend(mapSend);
  sVect<UInt> numRecv = mapCounter.countRecv(mapRecv);
  
  //Waiting receiving data  
  for(UInt sid=0; sid < maxPid; ++sid)
  {
    if( (sid != pid) && (numRecv[sid] != 0) )
    {
      recvO(sid,
            pid,
            commSegments[sid],
            recvPvps[sid],
            channel);
    }
  }
  
  //Waiting sending data
  for(UInt rid=0; rid < maxPid; ++rid)
  {
    if( (rid != pid) && (numSend[rid] != 0) )
    {
      sendO(pid,
            rid,
            commSegments[rid],
            sendPvps[rid],
            channel);
    }
  }
  
  //Merging
  rVect.clear();
  
  pVectManip<DATA,MAPITEM> manipulator;
  manipulator.mergeSimple(commSegments,rVect);
}


//_______________________________________________________________________________________________________
// MAPPED COMMUNICATIONS - RECURSIVE-COMM
//-------------------------------------------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvRR(const Teuchos::RCP<SENDRECV>    & mapSend,
           const Teuchos::RCP<const PVECT> & sendVect,
                 Teuchos::RCP<PVECT>       & rVect,
                 sVect<PVECT>              & sendSegments,
                 sVect<PVECT>              & recvSegments,
                 sVect<PVRS>               & sendPvrs,
                 sVect<PVRS>               & recvPvrs,
           const UInt                      & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRecvRR(*mapSend,
             *sendVect,
             *rVect,
              sendSegments,
              recvSegments,
              sendPvrs,
              recvPvrs,
              channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvRR(const SENDRECV     & mapSend,
           const PVECT        & sendVect,
                 PVECT        & rVect,
                 sVect<PVECT> & sendSegments,
                 sVect<PVECT> & recvSegments,
                 sVect<PVRS>  & sendPvrs,
                 sVect<PVRS>  & recvPvrs,
           const UInt         & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Data
  UInt pid    = commDev->rank();
  UInt maxPid = commDev->size();
  
  //Segmentazione dati
  pVectManip<DATA,MAPITEM> manipulator;
  manipulator.segmentationMap(sendSegments,mapSend,maxPid,sendVect);
  
  sVect<request> reqsRR;
  sendPvrs.resize(maxPid);
  recvPvrs.resize(maxPid);
   
  //Sending data  
  for(UInt rid=0; rid < maxPid; ++rid)
  {
    if(rid != pid)
    {
      sendRR(pid,
             rid,
             sendSegments[rid],
             sendPvrs[rid],
             reqsRR,
             channel);
    }
  }

  //Receiving data
  recvSegments.resize(maxPid);
  
  for(UInt sid=0; sid < maxPid; ++sid)
  {
    if(sid != pid)
    {
      recvRR(sid,
             pid,
             recvSegments[sid],
             recvPvrs[sid],
             channel);
    }
  }
  
  //Wait the send comm
  wait_all(reqsRR.begin(),reqsRR.end());
  
  //Merging
  rVect.clear();
  manipulator.mergeSimple(recvSegments,rVect);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvRI(const Teuchos::RCP<SENDRECV>    & mapSend,
           const Teuchos::RCP<const PVECT> & sendVect,
                 Teuchos::RCP<PVECT>       & rVect,
                 sVect<PVECT>              & sendSegments,
                 sVect<PVECT>              & recvSegments,
                 sVect<PVRS>               & sendPvrs,
                 sVect<PVRS>               & recvPvrs,
           const UInt                      & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRecvRI(*mapSend,
             *sendVect,
             *rVect,
              sendSegments,
              recvSegments,
              sendPvrs,
              recvPvrs,
              channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvRI(const SENDRECV     & mapSend,
           const PVECT        & sendVect,
                 PVECT        & rVect,
                 sVect<PVECT> & sendSegments,
                 sVect<PVECT> & recvSegments,
                 sVect<PVRS>  & sendPvrs,
                 sVect<PVRS>  & recvPvrs,
           const UInt         & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Data
  UInt pid    = commDev->rank();
  UInt maxPid = commDev->size();
  
  //Segmentazione dati
  pVectManip<DATA,MAPITEM> manipulator;
  manipulator.segmentationMapR(sendSegments,mapSend,maxPid,sendVect);
  
  //Sending data  
  for(UInt rid=0; rid < maxPid; ++rid)
  {
    if(rid != pid)
    {
      sendRI(pid,
             rid,
             sendSegments[rid],
             sendPvrs[rid],
             channel);
    }
  }
  
  //Receiving data
  for(UInt sid=0; sid < maxPid; ++sid)
  {
    if(sid != pid)
    {
      recvRI(sid,
             pid,
             recvSegments[sid],
             recvPvrs[sid],
             channel);
    }
  }
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvRO(const Teuchos::RCP<SENDRECV>    & mapSend,
           const Teuchos::RCP<const PVECT> & sendVect,
                 Teuchos::RCP<PVECT>       & rVect,
                 sVect<PVECT>              & sendSegments,
                 sVect<PVECT>              & recvSegments,
                 sVect<PVRS>               & sendPvrs,
                 sVect<PVRS>               & recvPvrs,
           const UInt                      & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRecvRO(*mapSend,
             *sendVect,
             *rVect,
              sendSegments,
              recvSegments,
              sendPvrs,
              recvPvrs,
              channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvRO(const SENDRECV     & mapSend,
           const PVECT        & sendVect,
                 PVECT        & rVect,
                 sVect<PVECT> & sendSegments,
                 sVect<PVECT> & recvSegments,
                 sVect<PVRS>  & sendPvrs,
                 sVect<PVRS>  & recvPvrs,
           const UInt         & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Data
  UInt pid    = commDev->rank();
  UInt maxPid = commDev->size();
  
  //Waiting receiving data  
  for(UInt sid=0; sid < maxPid; ++sid)
  {
    if(sid != pid)
    {
      recvRO(sid,
             pid,
             recvSegments[sid],
             recvPvrs[sid],
             channel);
    }
  }
  
  //Waiting sending data  
  for(UInt rid=0; rid < maxPid; ++rid)
  {
    if(rid != pid)
    {
      sendRO(pid,
             rid,
             sendSegments[rid],
             sendPvrs[rid],
             channel);
    }
  }
  
  //Merging
  pVectManip<DATA,MAPITEM> manipulator;
  manipulator.mergeSimpleR(recvSegments,rVect);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvRR(const Teuchos::RCP<SENDRECV>    & mapSend,
           const Teuchos::RCP<SENDRECV>    & mapRecv,
           const Teuchos::RCP<const PVECT> & sendVect,
                 Teuchos::RCP<PVECT>       & rVect,
                 sVect<PVECT>              & sendSegments,
                 sVect<PVECT>              & recvSegments,
                 sVect<PVRS>               & sendPvrs,
                 sVect<PVRS>               & recvPvrs,
           const UInt                      & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRecvRR(*mapSend,
             *mapRecv,
             *sendVect,
             *rVect,
              sendSegments,
              recvSegments,
              sendPvrs,
              recvPvrs,
              channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvRR(const SENDRECV     & mapSend,
           const SENDRECV     & mapRecv,
           const PVECT        & sendVect,
                 PVECT        & rVect,
                 sVect<PVECT> & sendSegments,
                 sVect<PVECT> & recvSegments,
                 sVect<PVRS>  & sendPvrs,
                 sVect<PVRS>  & recvPvrs,
           const UInt         & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Data
  UInt pid    = commDev->rank();
  UInt maxPid = commDev->size();
  
  //Data count
  pMapGlobalManip<pMapItemSendRecv> mapCounter(commDev);
  sVect<UInt> numSend = mapCounter.countSend(mapSend);
  sVect<UInt> numRecv = mapCounter.countRecv(mapRecv);
  
  //Segmentazione dati
  pVectManip<DATA,MAPITEM> manipulator;
  manipulator.segmentationMap(sendSegments,mapSend,maxPid,sendVect);
  
  sVect<request> reqsRR;
  sendPvrs.resize(maxPid);
  recvPvrs.resize(maxPid);
  
  //Sending data  
  for(UInt rid=0; rid < maxPid; ++rid)
  {
    if( (rid != pid) && (numSend[rid] != 0) )
    {
      sendRR(pid,
             rid,
             sendSegments[rid],
             sendPvrs[rid],
             reqsRR,
             channel);
    }
  }
  
  //Receiving data
  recvSegments.resize(maxPid);
  
  for(UInt sid=0; sid < maxPid; ++sid)
  {
    if( (sid != pid) && (numRecv[sid] != 0) )
    {
      recvRR(sid,
             pid,
             recvSegments[sid],
             recvPvrs[sid],
             channel);
    }
  }
  
  //Wait the send comm
  wait_all(reqsRR.begin(),reqsRR.end());

  //Merging
  rVect.clear();
  manipulator.mergeSimple(recvSegments,rVect);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvRI(const Teuchos::RCP<SENDRECV>    & mapSend,
           const Teuchos::RCP<SENDRECV>    & mapRecv,
           const Teuchos::RCP<const PVECT> & sendVect,
                 Teuchos::RCP<PVECT>       & rVect,
                 sVect<PVECT>              & sendSegments,
                 sVect<PVECT>              & recvSegments,
                 sVect<PVRS>               & sendPvrs,
                 sVect<PVRS>               & recvPvrs,
           const UInt                      & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRecvRI(*mapSend,
             *mapRecv,
             *sendVect,
             *rVect,
              sendSegments,
              recvSegments,
              sendPvrs,
              recvPvrs,
              channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvRI(const SENDRECV     & mapSend,
           const SENDRECV     & mapRecv,
           const PVECT        & sendVect,
                 PVECT        & rVect,
                 sVect<PVECT> & sendSegments,
                 sVect<PVECT> & recvSegments,
                 sVect<PVRS>  & sendPvrs,
                 sVect<PVRS>  & recvPvrs,
           const UInt         & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Data
  UInt pid    = commDev->rank();
  UInt maxPid = commDev->size();
  
  //Data count
  pMapGlobalManip<pMapItemSendRecv> mapCounter(commDev);
  sVect<UInt> numSend = mapCounter.countSend(mapSend);
  sVect<UInt> numRecv = mapCounter.countRecv(mapRecv);
  
  //Segmentazione dati
  pVectManip<DATA,MAPITEM> manipulator;
  manipulator.segmentationMapR(sendSegments,mapSend,maxPid,sendVect);
  
  //Sending data  
  for(UInt rid=0; rid < maxPid; ++rid)
  {
    if( (rid != pid) && (numSend[rid] != 0) )
    {
      sendRI(pid,
             rid,
             sendSegments[rid],
             sendPvrs[rid],
             channel);
    }
  }
  
  //Receiving data
  for(UInt sid=0; sid < maxPid; ++sid)
  {
    if( (sid != pid) && (numRecv[sid] != 0) )
    {
      recvRI(sid,
             pid,
             recvSegments[sid],
             recvPvrs[sid],
             channel);
    }
  }
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvRO(const Teuchos::RCP<SENDRECV>    & mapSend,
           const Teuchos::RCP<SENDRECV>    & mapRecv,
           const Teuchos::RCP<const PVECT> & sendVect,
                 Teuchos::RCP<PVECT>       & rVect,
                 sVect<PVECT>              & sendSegments,
                 sVect<PVECT>              & recvSegments,
                 sVect<PVRS>               & sendPvrs,
                 sVect<PVRS>               & recvPvrs,
           const UInt                      & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  sendRecvRO(*mapSend,
             *mapRecv,
             *sendVect,
             *rVect,
              sendSegments,
              recvSegments,
              sendPvrs,
              recvPvrs,
              channel);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
sendRecvRO(const SENDRECV     & mapSend,
           const SENDRECV     & mapRecv,
           const PVECT        & sendVect,
                 PVECT        & rVect,
                 sVect<PVECT> & sendSegments,
                 sVect<PVECT> & recvSegments,
                 sVect<PVRS>  & sendPvrs,
                 sVect<PVRS>  & recvPvrs,
           const UInt         & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Data
  UInt pid    = commDev->rank();
  UInt maxPid = commDev->size();
  
  //Data count
  pMapGlobalManip<pMapItemSendRecv> mapCounter(commDev);
  sVect<UInt> numSend = mapCounter.countSend(mapSend);
  sVect<UInt> numRecv = mapCounter.countRecv(mapRecv);
  
  //Waiting receiving data  
  for(UInt sid=0; sid < maxPid; ++sid)
  {
    if( (sid != pid) && (numRecv[sid] != 0) )
    {
      recvRO(sid,
             pid,
             recvSegments[sid],
             recvPvrs[sid],
             channel);
    }
  }
  
  //Waiting sending data  
  for(UInt rid=0; rid < maxPid; ++rid)
  {
    if( (rid != pid) && (numSend[rid] != 0) )
    {
      sendRO(pid,
             rid,
             sendSegments[rid],
             sendPvrs[rid],
             channel);
    }
  }
  
  //Merging
  pVectManip<DATA,MAPITEM> manipulator;
  manipulator.mergeSimpleR(recvSegments,rVect);
}


//_______________________________________________________________________________________________________
// PATTERN COMMUNICATIONS - NON-PENDING-COMM
//-------------------------------------------------------------------------------------------------------
template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
vectorNormal(Teuchos::RCP<PVECT> & Vect,
             const UInt & MaxGid)
{
  //Assert
  assert(commDevLoaded);

  //Comm
  vectorNormal(*Vect,MaxGid);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
vectorNormal(PVECT & Vect,
             const UInt & MaxGid)
{
  //Assert
  assert(commDevLoaded);
  
  //Alloc
  UInt proc    = commDev->rank();
  UInt numProc = commDev->size();
  
  //MaxGid
  pMapGlobalManip<MAPITEM> pMapGlobalManip(commDev);
  UInt maxGid = max(pMapGlobalManip.sizeG(Vect.getMapRcp()), MaxGid);
  
  //Segmenting
  sVect<PVECT>  sendVects;
  sVect<PVECT>  recvVects;
  PVECT         tempVect;
  
  pVectManip<DATA,MAPITEM> manipulator;
  manipulator.segmentationNormal(sendVects,numProc,maxGid,Vect);
  
  //Sending
  sVect<request> sendReqs;
  
  for(UInt k=0; k < numProc; ++k)
  {
    if(k != proc)
    {
      isend(proc,k,sendVects[k],sendReqs);
    }
  }
  
  //Reciving 
  recvVects.reserve(numProc);
  
  for(UInt k=0; k < numProc; ++k)
  {
    if(k != proc)
    {
      irecv(k,proc,tempVect);
      recvVects.push_back(tempVect);
    }
  }
  
  //Waiting
  wait_all(sendReqs.begin(),sendReqs.end());
  
  //Merging
  recvVects.push_back(sendVects[proc]);
  
  Vect.clear();
  manipulator.mergeSimple(recvVects,Vect);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
vectorPid(Teuchos::RCP<PVECT> & Vect)
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  vectorPid(*Vect);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
vectorPid(PVECT & Vect)
{
  //Assert
  assert(commDevLoaded);
  
  //Alloc
  UInt proc    = commDev->rank();
  UInt numProc = commDev->size();
  
  //Segmenting
  sVect<PVECT>  sendVects;
  sVect<PVECT>  recvVects;
  PVECT         tempVect;
  
  pVectManip<DATA,MAPITEM> manipulator;
  manipulator.segmentationPid(sendVects,numProc,Vect);
    
  //Sending
  sVect<request> sendReqs;
  
  for(UInt k=0; k < numProc; ++k)
  {
    if(k != proc)
    {
      isend(proc,k,sendVects[k],sendReqs);
    }
  }
  
  //Receiving 
  recvVects.reserve(numProc);
  
  for(UInt k=0; k < numProc; ++k)
  {
    if(k != proc)
    {
      irecv(k,proc,tempVect);
      recvVects.push_back(tempVect);
    }
  }
  
  //Waiting
  wait_all(sendReqs.begin(),sendReqs.end());
  
  //Merging
  recvVects.push_back(sendVects[proc]);
  
  Vect.clear();
  manipulator.mergeSimple(recvVects,Vect);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
vectorData(Teuchos::RCP<PVECT> & Vect)
{
  //Assert
  assert(commDevLoaded);
  
  //Comm
  vectorData(*Vect);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
vectorData(PVECT & Vect)
{
  //Assert
  assert(commDevLoaded);
  
  //Allocations
  UInt numProc = commDev->size();
  UInt proc    = commDev->rank();
  
  //Local maximum and minimum data determination
  DATA dataMax, dataMin;
  UInt dataLoaded = (Vect.size() >= 1);
  
  dataMax = dataMin;
  
  if(dataLoaded)
  {
    dataMin = Vect.getL(1);
    dataMax = Vect.getL(1);
    
    for(UInt i=1; i <= Vect.size(); ++i)
    {
      dataMin = std::min(dataMin,Vect.getL(i));
      dataMax = std::max(dataMax,Vect.getL(i));
    }
  }
  
  //Global maximum and minimum data determination
  sVect<UInt> allLoading(numProc);
  all_gather(*commDev,dataLoaded,allLoading);

  sVect<DATA> allMin(numProc);
  morgana_all_gather(*commDev,dataMin,allMin);
  
  sVect<DATA> allMax(numProc);
  morgana_all_gather(*commDev,dataMax,allMax);
  
  for(UInt i=1; i<=numProc; ++i)
  {
    if(allLoading(i))
    {
      dataMin = allMin(i);
      dataMax = allMax(i);
    }
  }
    
  for(UInt i=1; i<=numProc; ++i)
  {    
    if(allLoading(i))
    {
      dataMin = min(dataMin,allMin(i));
      dataMax = max(dataMax,allMax(i));
    }
  }
  
  //Data segmentation
  sVect<PVECT>  sendVects;
  sVect<PVECT>  recvVects;
  PVECT         tempVect;
  
  pVectManip<DATA,MAPITEM> manipulator;
  manipulator.segmentationData(sendVects,numProc,dataMin,dataMax,Vect);
   
  //Sending
  sVect<request> sendReqs;
  
  for(UInt k=0; k < numProc; ++k)
  {
    if(k != proc)
    {
      isend(proc,k,sendVects[k],sendReqs);
    }
  }
  
  //Reciving 
  recvVects.reserve(numProc);
  
  for(UInt k=0; k < numProc; ++k)
  {
    if(k != proc)
    {
      irecv(k,proc,tempVect);
      recvVects.push_back(tempVect);
    }
  }
  
  //Waiting
  wait_all(sendReqs.begin(),sendReqs.end());
  
  //Merging
  recvVects.push_back(sendVects[proc]);
  
  Vect.clear();
  manipulator.mergeSimple(recvVects,Vect);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
vectorData(Teuchos::RCP<PVECT> & Vect,
           const DATA & dataMin,
           const DATA & dataMax)
{
  //Assert
  assert(commDevLoaded);
  assert( (dataMin < dataMax) || (!(dataMin != dataMax)) );
  
  //Comm
  vectorData(*Vect,dataMin,dataMax);
}

template<typename DATA, typename MAPITEM>
void
pVectComm<DATA,MAPITEM>::
vectorData(PVECT & Vect,
           const DATA & dataMin,
           const DATA & dataMax)
{
  //Assert
  assert(commDevLoaded);
  assert( (dataMin < dataMax) || (!(dataMin != dataMax)) );

  //Checking
  for(UInt i=1; i <= Vect.size(); ++i)
  {
    assert( (Vect(i) < dataMax) || (!(Vect(i) != dataMax)) );
    assert( (dataMin < Vect(i)) || (!(Vect(i) != dataMin)) );
  }
  
  //Allocations
  UInt numProc = commDev->size();
  UInt proc    = commDev->rank();
      
  //Data segmentation
  sVect<PVECT>  sendVects;
  sVect<PVECT>  recvVects;
  PVECT         tempVect;
  
  pVectManip<DATA,MAPITEM> manipulator;
  manipulator.segmentationData(sendVects,numProc,dataMin,dataMax,Vect);
  
  //Sending
  sVect<request> sendReqs;
  
  for(UInt k=0; k < numProc; ++k)
  {
    if(k != proc)
    {
      isend(proc,k,sendVects[k],sendReqs);
    }
  }
  
  //Reciving 
  recvVects.reserve(numProc);
  
  for(UInt k=0; k < numProc; ++k)
  {
    if(k != proc)
    {
      irecv(k,proc,tempVect);
      recvVects.push_back(tempVect);
    }
  }
  
  //Waiting
  wait_all(sendReqs.begin(),sendReqs.end());
  
  //Merging
  recvVects.push_back(sendVects[proc]);
  
  Vect.clear();
  manipulator.mergeSimple(recvVects,Vect);
}

#endif
