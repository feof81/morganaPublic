/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PMAPGLOBALMANIP_H
#define PMAPGLOBALMANIP_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Version.hpp"
#include "Kokkos_DefaultNode.hpp"

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "pMapItem.h"
#include "pMapItemShare.h"
#include "pMapItemSendRecv.h"
#include "pMapManip.hpp"
#include "pMapComm.hpp"

using namespace std;
using namespace boost::mpi;


//_________________________________________________________________________________________________
// UNSPECIALIZED
//-------------------------------------------------------------------------------------------------

/*! Unspecialized - empty */
template<typename ITEM> class pMapGlobalManip
{ };


//_________________________________________________________________________________________________
// PMAPITEM
//-------------------------------------------------------------------------------------------------

/*! Perform global manipulations, retrive informations and checks the \c pMap. Specialized version for \c pMapItem */
template<> class pMapGlobalManip<pMapItem>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItem  ITEM;
    typedef pMap<ITEM> MAP;
    
    typedef Tpetra::global_size_t                    TPETRA_GLOBAL_TYPE;
    typedef Tpetra::Map<>                            TPETRA_MAP;
    typedef typename TPETRA_MAP::global_ordinal_type TPETRA_GLOBAL_ORDINAL;
    typedef Teuchos::Comm<int>                       TPETRA_COMM;
    typedef Teuchos::MpiComm<int>                    TPETRA_MPICOMM;
    //@}
  
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
    
    
    /*! @name Constructors and set functions */ //@{
  public:
    pMapGlobalManip();
    pMapGlobalManip(const Teuchos::RCP<const communicator> & CommDev);
    pMapGlobalManip(const communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    //@}
    
    
    /*! @name Data retriving functions */ //@{
  public:
    /*! Returns the maximum \c gid on all the processors */
    UInt sizeG(const Teuchos::RCP<const MAP> & Map) const;
    
    /*! Returns the maximum \c gid on all the processors */
    UInt sizeG(const MAP & Map) const;
    //@}
    
    
    /*! @name Global manipulations */ //@{
  public:
    /*! Checks that there is a consecutive numeration of the Gids and that all the gids are positive  */
    bool check(const Teuchos::RCP<const MAP> & Map) const;
    
    /*! Checks that there is a consecutive numeration of the Gids and that all the gids are positive  */
    bool check(const MAP & Map) const;
    
    /*! Modify the map eliminating all the repetitions */
    void destroyOverlap(Teuchos::RCP<MAP> & Map) const;
    
    /*! Modify the map eliminating all the repetitions */
    void destroyOverlap(MAP & Map) const;
    
    /*! Export the map to an \c Epetra_Map
    \param Map the input map
    \param epetraMap the output map 
    \param base the base  of the output map */
    void exportEpetraMap(const Teuchos::RCP<const MAP> & Map,
                         Teuchos::RCP<Epetra_Map>      & epetraMap,
                         Epetra_MpiComm                & epetraComm,
                         const UInt                    & base) const;
    
    /*! Export the map to an \c Epetra_Map
    \param Map the input map
    \param epetraMap the output map 
    \param base the base  of the output map */
    void exportEpetraMap(const MAP                & Map,
                         Teuchos::RCP<Epetra_Map> & epetraMap,
                         Epetra_MpiComm           & epetraComm,
                         const UInt               & base) const;
    
    /*! Export the map to an \c Tpetra_Map
    \param Map the input map
    \param tpetraMap the output map */
    void exportTpetraMap(const Teuchos::RCP<const MAP>  & Map,
                         Teuchos::RCP<const TPETRA_MAP> & tpetraMap,
                         const UInt                     & base = 0) const;
    
    /*! Export the map to an \c Tpetra_Map
    \param Map the input map
    \param tpetraMap the output map */
    void exportTpetraMap(const MAP                      & Map,
                         Teuchos::RCP<const TPETRA_MAP> & tpetraMap,
                         const UInt                     & base = 0) const;
    
    /*! Import the map from an \c Epetra_Map
    \param EpetraMap the input map 
    \param Map the output map*/
    void importEpetraMap(const Teuchos::RCP<const Epetra_Map> & EpetraMap,
                         Teuchos::RCP<MAP>                    & Map);
    
    /*! Import the map from an \c Epetra_Map
    \param EpetraMap the input map 
    \param Map the output map*/
    void importEpetraMap(const Epetra_Map  & EpetraMap,
                         Teuchos::RCP<MAP> & Map);
    
    /*! Import the map from an \c Tpetra_Map
    \param TpetraMap the input map 
    \param Map the output map*/
    void importTpetraMap(const Teuchos::RCP<const TPETRA_MAP> & TpetraMap,
                         Teuchos::RCP<MAP>                    & Map);
    //@}
    
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MAP>          & OldMap,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MAP>                & NewMap);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MAP          & OldMap,
                            const communicator & NewCommDev,
                                  MAP          & NewMap);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MAP>          & OldMap,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MAP>                & NewMap);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MAP          & OldMap,
                            const communicator & NewCommDev,
                                  MAP          & NewMap);
    //@}
    
    /*! @name Other functions */ //@{
  public:
    size_t memSize() const;
    //@}
};



//_________________________________________________________________________________________________
// PMAPITEMSHARE
//-------------------------------------------------------------------------------------------------

/*! Perform global manipulations, retrive informations and checks the \c pMap. Specialized version for \c pMapItemShare */
template<> class pMapGlobalManip<pMapItemShare>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItemShare ITEM;
    typedef pMap<ITEM>    MAP;
    typedef pMap<pMapItemSendRecv> SENDRECV;
    
    typedef Tpetra::global_size_t                    TPETRA_GLOBAL_TYPE;
    typedef Tpetra::Map<>                            TPETRA_MAP;
    typedef typename TPETRA_MAP::global_ordinal_type TPETRA_GLOBAL_ORDINAL;
    typedef Teuchos::Comm<int>                       TPETRA_COMM;
    typedef Teuchos::MpiComm<int>                    TPETRA_MPICOMM;
    //@}
  
    
    /*! @name Internal data and links */ //@{
  public:
    UInt maxSize;
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
    
    
    /*! @name Constructors and set functions */ //@{
  public:
    pMapGlobalManip();
    pMapGlobalManip(const Teuchos::RCP<const communicator> & CommDev);
    pMapGlobalManip(const communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    //@}
    
    
    /*! @name Data retriving functions */ //@{
  public:
    /*! Returns the maximum \c gid on all the processors */
    UInt sizeG(const Teuchos::RCP<const MAP> & Map) const;
    
    /*! Returns the maximum \c gid on all the processors */
    UInt sizeG(const MAP & Map) const;
    
    /*! Returns the number of locally shared items */
    UInt sharedL(const Teuchos::RCP<const MAP> & Map) const;
    
    /*! Returns the number of locally shared items */
    UInt sharedL(const MAP & Map) const;
    
    /*! Returns the number of globally shared items */
    UInt sharedG(const Teuchos::RCP<const MAP> & Map) const;
    
    /*! Returns the number of globally shared items */
    UInt sharedG(const MAP & Map) const;
    
    /*! Returns the number of locally owned items */
    UInt ownedL(const Teuchos::RCP<const MAP> & Map) const;
    
    /*! Returns the number of locally owned items */
    UInt ownedL(const MAP & Map) const;
    //@}
    
    
    /*! @name Global manipulations */ //@{
  public:
    /*! Checks that: all the \c gid s are present (from 1 to maxGid), any shared object is owned by only one item, not owned items are shared,
    owned-and-shared items are associated to other notOwned-and-shared and viceversa */
    bool check(const Teuchos::RCP<const MAP> & Map) const;
    
    /*! Checks that: all the \c gid s are present (from 1 to maxGid), any shared object is owned by only one item, not owned items are shared,
    owned-and-shared items are associated to other notOwned-and-shared and viceversa */
    bool check(const MAP & Map) const;
    
    /*! Modify the map eliminating all the repetitions */
    void destroyOverlap(Teuchos::RCP<MAP> & Map) const;
    
    /*! Modify the map eliminating all the repetitions */
    void destroyOverlap(MAP & Map) const;
    
    /*! Every owned and shared element updates all the other shared and not owned items.
    This function creates the communication patter for this update */
    void createSendRecvMap(const Teuchos::RCP<const MAP> & Map,
                           Teuchos::RCP<SENDRECV>        & mapSend,
                           Teuchos::RCP<SENDRECV>        & mapRecv) const;
    
    /*! Every owned and shared element updates all the other shared and not owned items.
    This function creates the communication patter for this update */
    void createSendRecvMap(const MAP & Map,
                           SENDRECV  & mapSend,
                           SENDRECV  & mapRecv) const;
    
    /*! Fix the owning-sharing, based on the gids. Assumes that the \c gids are correct and updates the 
    owning-sharing structure. Multiple gids are fixed so that only one has the ownership, the ones that are
    already ok are retained */
    void updateOwningSharing(Teuchos::RCP<MAP> & Map) const;
    
    /*! Fix the owning-sharing, based on the gids. Assumes that the \c gids are correct and updates the 
    owning-sharing structure. Multiple gids are fixed so that only one has the ownership, the ones that are
    already ok are retained */
    void updateOwningSharing(MAP & Map) const;
    
    /*! Fix the sharing, based on the \c gids. Assumes that both the \c gids and the \c owned flag are correct and 
    updates the \c shared flag only */
    void updateSharing(Teuchos::RCP<MAP> & Map) const;
    
    /*! Fix the sharing, based on the \c gids. Assumes that both the \c gids and the \c owned flag are correct and 
    updates the \c shared flag only */
    void updateSharing(MAP & Map) const;
    
    /*! Export the map to an \c Epetra_Map
    \param Map the input map
    \param epetraMap the output map 
    \param base the base  of the output map */
    void exportEpetraMap(const Teuchos::RCP<const MAP> & Map,
                         Teuchos::RCP<Epetra_Map>      & epetraMap,
                         Epetra_MpiComm                & epetraComm,
                         const UInt                    & base) const;
    
    /*! Export the map to an \c Epetra_Map
    \param Map the input map
    \param epetraMap the output map 
    \param base the base  of the output map */
    void exportEpetraMap(const MAP                & Map,
                         Teuchos::RCP<Epetra_Map> & epetraMap,
                         Epetra_MpiComm           & epetraComm,
                         const UInt               & base) const;
    
    /*! Export the map to an \c Tpetra_Map
    \param Map the input map
    \param tpetraMap the output map */
    void exportTpetraMap(const Teuchos::RCP<const MAP>  & Map,
                         Teuchos::RCP<const TPETRA_MAP> & tpetraMap,
                         const UInt                     & base = 0) const;
    
    /*! Export the map to an \c Tpetra_Map
    \param Map the input map
    \param tpetraMap the output map */
    void exportTpetraMap(const MAP                      & Map,
                         Teuchos::RCP<const TPETRA_MAP> & tpetraMap,
                         const UInt                     & base = 0) const;
    
    /*! Import the map from an \c Epetra_Map
    \param EpetraMap the input map 
    \param Map the output map*/
    void importEpetraMap(const Teuchos::RCP<const Epetra_Map> & EpetraMap,
                         Teuchos::RCP<MAP>                    & Map);
    
    /*! Import the map from an \c Epetra_Map
    \param EpetraMap the input map 
    \param Map the output map*/
    void importEpetraMap(const Epetra_Map  & EpetraMap,
                         Teuchos::RCP<MAP> & Map);
    
    /*! Import the map from an \c Tpetra_Map
    \param TpetraMap the input map 
    \param Map the output map*/
    void importTpetraMap(const Teuchos::RCP<const TPETRA_MAP> & TpetraMap,
                         Teuchos::RCP<MAP>                    & Map);
    //@}
    
    
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MAP>          & OldMap,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MAP>                & NewMap);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MAP          & OldMap,
                            const communicator & NewCommDev,
                                  MAP          & NewMap);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MAP>          & OldMap,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MAP>                & NewMap);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MAP          & OldMap,
                            const communicator & NewCommDev,
                                  MAP          & NewMap);
    //@}
    
    /*! @name Other functions */ //@{
  public:
    size_t memSize() const;
    //@}
};



//_________________________________________________________________________________________________
// PMAPITEMSENDRECV
//-------------------------------------------------------------------------------------------------
template<> class pMapGlobalManip<pMapItemSendRecv>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItemSendRecv       ITEM;
    typedef pMap<pMapItemSendRecv> MAP;
    //@}
  
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
  
    /*! @name Constructors and set functions */ //@{
  public:
    pMapGlobalManip();
    pMapGlobalManip(const Teuchos::RCP<const communicator> & CommDev);
    pMapGlobalManip(const communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    //@}
    
    /*! @name Global manipulations */ //@{
  public:
    /*! Consecutive numeration of the Gids, all the gids are positive  */
    bool check(const Teuchos::RCP<const MAP> & mapSend,
               const Teuchos::RCP<const MAP> & mapRecv) const;
    
    /*! Consecutive numeration of the Gids, all the gids are positive  */
    bool check(const MAP & mapSend,
               const MAP & mapRecv) const;
    
    /*! Returns a vector with size = maxPid containing for each element the number of items to be sent to the i-th \c pid */
    sVect<UInt> countSend(const Teuchos::RCP<const MAP> & mapSend) const;
    
    /*! Returns a vector with size = maxPid containing for each element the number of items to be sent to the i-th \c pid */
    sVect<UInt> countSend(const MAP & mapSend) const;
    
    /*! Returns a vector with size = maxPid containing for each element the number of items to be received from the i-th \c pid */
    sVect<UInt> countRecv(const Teuchos::RCP<const MAP> & mapRecv) const;
    
    /*! Returns a vector with size = maxPid containing for each element the number of items to be received from the i-th \c pid */
    sVect<UInt> countRecv(const MAP & mapRecv) const;
    //@}
    
    
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MAP>          & OldMap,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MAP>                & NewMap);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MAP          & OldMap,
                            const communicator & NewCommDev,
                                  MAP          & NewMap);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const MAP>          & OldMap,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<MAP>                & NewMap);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const MAP          & OldMap,
                            const communicator & NewCommDev,
                                  MAP          & NewMap);
    //@}
    
    /*! @name Other functions */ //@{
  public:
    size_t memSize() const;
    //@}
};


#endif

