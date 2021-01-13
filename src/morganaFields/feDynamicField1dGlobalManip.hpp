/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FEDYNAMICFIELD1DGLOBALMANIP_HPP
#define FEDYNAMICFIELD1DGLOBALMANIP_HPP

#include "pMapItem.h"
#include "pMapItemShare.h"

#include "feDynamicField1d.hpp"

using namespace std;


//_________________________________________________________________________________________________
// SERVO-SPECIALIZATION
//-------------------------------------------------------------------------------------------------

/*! Dynamic Field 1d Global Manipulator - empty general class  */
template<typename PMAPTYPE, typename FETYPE, typename DOFTYPE, dmd1d_order ORDER = dmd1d_vectMajor, dmd1d_mode MODE = dmd1d_standard>
class fdf1dGlobalManip;



/*! Dynamic Field 3d Global Manipulator - \c pMapItem specialization */
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
class fdf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItem                                    PMAPTYPE;
    typedef dofMapDynamic1d_options                     OPTIONS;
    typedef feDynamicField1d<FETYPE,DOFTYPE,ORDER,MODE> FIELD;
    
    typedef typename FIELD::MESH1D     MESH1D;
    typedef typename FIELD::CONNECT1D  CONNECT1D;
    typedef typename FIELD::DOFVECT    DOFVECT;
    typedef typename FIELD::GEOSHAPE   GEOSHAPE;
    typedef typename FETYPE::ELCARD    ELCARD;
    typedef typename FETYPE::FECARD    FECARD;
    typedef pVect<FECARD,PMAPTYPE>     FECARDS;
    //@}
  
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    fdf1dGlobalManip();
    fdf1dGlobalManip(const Teuchos::RCP<communicator> & CommDev);
    fdf1dGlobalManip(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Operator functions */ //@{
  public:
    void changeGrid(const Teuchos::RCP<MESH1D>    & newGrid,
                    const Teuchos::RCP<CONNECT1D> & newConnect,
                    const Teuchos::RCP<OPTIONS>   & options,
                    const Teuchos::RCP<FIELD>     & Field) const;
    
    void changeGrid(MESH1D    & newGrid,
                    CONNECT1D & newConnect,
                    OPTIONS   & options,
                    FIELD     & Field) const;
    
    bool matchGrid(const Teuchos::RCP<MESH1D>    & newGrid1d,
                   const Teuchos::RCP<CONNECT1D> & newConnect1d,
                   const Teuchos::RCP<OPTIONS>   & options,
                   const Teuchos::RCP<FIELD>     & field) const;
    
    bool matchGrid(MESH1D    & newGrid1d,
                   CONNECT1D & newConnect1d,
                   OPTIONS   & options,
                   FIELD     & field) const;
    //@}
   
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const FIELD>        & OldField,
                            const Teuchos::RCP<communicator>       & NewCommDev,
                            const Teuchos::RCP<MESH1D>             & NewGrid,
                            const Teuchos::RCP<CONNECT1D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                    & isActive,
                            const communicator            & OldCommDev,
                            const FIELD                   & OldField,
                                  communicator            & NewCommDev,
                            const Teuchos::RCP<MESH1D>    & NewGrid,
                            const Teuchos::RCP<CONNECT1D> & NewConnect,
                                  FIELD                   & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const FIELD>        & OldField,
                            const Teuchos::RCP<communicator>       & NewCommDev,
                            const Teuchos::RCP<MESH1D>             & NewGrid,
                            const Teuchos::RCP<CONNECT1D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                    & isActive,
                            const communicator            & OldCommDev,
                            const FIELD                   & OldField,
                                  communicator            & NewCommDev,
                            const Teuchos::RCP<MESH1D>    & NewGrid,
                            const Teuchos::RCP<CONNECT1D> & NewConnect,
                                  FIELD                   & NewField);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
fdf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
fdf1dGlobalManip()
{
  commDevLoaded = false;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
fdf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
fdf1dGlobalManip(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
fdf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
fdf1dGlobalManip(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}


//_________________________________________________________________________________________________
// OPERATOR FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(const Teuchos::RCP<MESH1D>    & newGrid,
           const Teuchos::RCP<CONNECT1D> & newConnect,
           const Teuchos::RCP<OPTIONS>   & options,
           const Teuchos::RCP<FIELD>     & Field) const
{
  assert(commDevLoaded);
  
  //Create the new field
  FIELD newField1d;
  newField1d.setCommunicator(commDev);
  newField1d.setGeometry(newGrid,newConnect);
  newField1d.setOptions(options);
  
  //Remap the feCards
  pMap<PMAPTYPE> newElMap = newGrid->getElements().getRowMap();
  FECARDS         feCards = Field->getFeCards();
  
  pVectGlobalManip<FECARD,PMAPTYPE> feCardsManip(commDev);
  feCardsManip.changeMap(feCards,newElMap);
  
  newField1d.setFeCards(feCards);
  
  //Startup
  newField1d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField1d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field->getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap));
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);
  
  //Loading the new dofVector
  newField1d.setDofVect(oldDofVect);
  *Field = newField1d;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(MESH1D    & newGrid,
           CONNECT1D & newConnect,
           OPTIONS   & options,
           FIELD     & Field) const
{
  assert(commDevLoaded);
  
  //Create the new field
  FIELD newField1d;
  newField1d.setCommunicator(commDev);
  newField1d.setGeometry(newGrid,newConnect);
  newField1d.setOptions(options);
  
  //Remap the feCards
  pMap<PMAPTYPE> newElMap = newGrid.getElements().getRowMap();
  FECARDS         feCards = Field.getFeCards();
  
  pVectGlobalManip<FECARD,PMAPTYPE> feCardsManip(commDev);
  
  assert(feCardsManip.check(feCards));
  feCardsManip.changeMap(feCards,newElMap); 
  
  newField1d.setFeCards(feCards);  
  
  //Startup
  newField1d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField1d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field.getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap)); 
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);  

  //Loading the new dofVector
  newField1d.setDofVect(oldDofVect);
  Field = newField1d;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
bool
fdf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(const Teuchos::RCP<MESH1D>    & newGrid1d,
          const Teuchos::RCP<CONNECT1D> & newConnect1d,
          const Teuchos::RCP<OPTIONS>   & options,
          const Teuchos::RCP<FIELD>     & field) const
{
  return( changeGrid(newGrid1d,
                     newConnect1d,
                     options,
                     field) );
}
    
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
bool
fdf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(MESH1D    & newGrid1d,
          CONNECT1D & newConnect1d,
          OPTIONS   & options,
          FIELD     & field) const
{
  changeGrid(newGrid1d,
             newConnect1d,
             options,
             field);
  
  return(true);
}


//_________________________________________________________________________________________________
// COMMUNICATOR MANIPULATIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH1D>             & NewGrid,
                   const Teuchos::RCP<CONNECT1D>          & NewConnect,
                         Teuchos::RCP<FIELD>              & NewField)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldField,
                    *NewCommDev,
                     NewGrid,
                     NewConnect,
                    *NewField);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH1D>    & NewGrid,
                   const Teuchos::RCP<CONNECT1D> & NewConnect,
                         FIELD                   & NewField)
{
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::DOFVECT DOFVECT;
  typedef typename FIELD::FECARDS FECARDS;
  typedef typename FIELD::FECARD  FECARD;
  
  typedef pVectGlobalManip<DOFTYPE,DOFTYPE>  MANIP_DOFVECT;
  typedef pVectGlobalManip<FECARD,DOFTYPE>   MANIP_FECARDS;
  
  //Manipulators-----------------------------------------------------
  DOFVECT newDofVect;
  FECARDS newFeCards;

  MANIP_DOFVECT manipDofVect;
  MANIP_FECARDS manipFeCards;
  
  OPTIONS options = OldField.getOptions();
  
  //Copy to new communicator-----------------------------------------
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() < OldCommDev.size());
    
    //Reduce comm
    manipDofVect.reduceCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getDofVect(),
                                    NewCommDev,
                                    newDofVect);
    
    manipFeCards.reduceCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getFeCards(),
                                    NewCommDev,
                                    newFeCards);
    
    //New Field
    NewField.setCommunicator(NewCommDev);
    NewField.setGeometry(NewGrid,NewConnect);
    NewField.setOptions(options);
    NewField.setFeCards(newFeCards);
    NewField.startup();
    NewField.setDofVect(newDofVect);
  }
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH1D>             & NewGrid,
                   const Teuchos::RCP<CONNECT1D>          & NewConnect,
                         Teuchos::RCP<FIELD>              & NewField)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldField,
                    *NewCommDev,
                     NewGrid,
                     NewConnect,
                    *NewField);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH1D>    & NewGrid,
                   const Teuchos::RCP<CONNECT1D> & NewConnect,
                         FIELD                   & NewField)
{
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::DOFVECT DOFVECT;
  typedef typename FIELD::FECARDS FECARDS;
  typedef typename FIELD::FECARD  FECARD;
  
  typedef pVectGlobalManip<DOFTYPE,DOFTYPE>  MANIP_DOFVECT;
  typedef pVectGlobalManip<FECARD,DOFTYPE>   MANIP_FECARDS;
  
  //Manipulators-----------------------------------------------------
  DOFVECT newDofVect;
  FECARDS newFeCards;

  MANIP_DOFVECT manipDofVect;
  MANIP_FECARDS manipFeCards;
  
  OPTIONS options = OldField.getOptions();
  
  //Copy to new communicator-----------------------------------------
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() > OldCommDev.size());
    
    //Reduce comm
    manipDofVect.expandCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getDofVect(),
                                    NewCommDev,
                                    newDofVect);
    
    manipFeCards.expandCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getFeCards(),
                                    NewCommDev,
                                    newFeCards);
    
    //New Field
    NewField.setCommunicator(NewCommDev);
    NewField.setGeometry(NewGrid,NewConnect);
    NewField.setOptions(options);
    NewField.setFeCards(newFeCards);
    NewField.startup();
    NewField.setDofVect(newDofVect);
  }
}



//_________________________________________________________________________________________________
// PMAPITEMSHARE SPECIALIZATION
//-------------------------------------------------------------------------------------------------

/*! Dynamic Field 3d Global Manipulator - \c pMapItemShare specialization */
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
class fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItemShare                               PMAPTYPE;
    typedef feDynamicField1d<FETYPE,DOFTYPE,ORDER,MODE> FIELD;
    typedef pMap<pMapItemSendRecv>                      SENDRECV;
    typedef dofMapDynamic1d_options                     OPTIONS;
    
    typedef typename FIELD::MESH1D     MESH1D;
    typedef typename FIELD::CONNECT1D  CONNECT1D;
    typedef typename FIELD::DOFVECT    DOFVECT;
    typedef typename FIELD::GEOSHAPE   GEOSHAPE;
    typedef typename FETYPE::ELCARD    ELCARD;
    typedef typename FETYPE::FECARD    FECARD;
    
    typedef pVect<FECARD,PMAPTYPE>                      FECARDS;
    typedef pVectGlobalManip<DOFTYPE,pMapItemShare>     PVGLOBMANIP;
    typedef typename PVGLOBMANIP::PVPS                  PVPS;
    typedef typename PVGLOBMANIP::PVUR                  PVUR;
    //@}
  
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    fdf1dGlobalManip();
    fdf1dGlobalManip(const Teuchos::RCP<communicator> & CommDev);
    fdf1dGlobalManip(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Operator functions */ //@{
  public:
    void changeGrid(const Teuchos::RCP<MESH1D>    & newGrid,
                    const Teuchos::RCP<CONNECT1D> & newConnect,
                    const Teuchos::RCP<OPTIONS>   & options,
                    const Teuchos::RCP<FIELD>     & Field) const;
    
    void changeGrid(MESH1D    & newGrid,
                    CONNECT1D & newConnect,
                    OPTIONS   & options,
                    FIELD     & Field) const;
    
    bool matchGrid(const Teuchos::RCP<MESH1D>    & newGrid1d,
                   const Teuchos::RCP<CONNECT1D> & newConnect1d,
                   const Teuchos::RCP<OPTIONS>   & options,
                   const Teuchos::RCP<FIELD>     & field) const;
    
    bool matchGrid(MESH1D    & newGrid1d,
                   CONNECT1D & newConnect1d,
                   OPTIONS   & options,
                   FIELD     & field) const;
    
    void createSendRecvMap(const Teuchos::RCP<FIELD>    & Field,
                                 Teuchos::RCP<SENDRECV> & mapSend,
                                 Teuchos::RCP<SENDRECV> & mapRecv) const;
    
    void createSendRecvMap(const FIELD    & Field,
                                 SENDRECV & mapSend,
                                 SENDRECV & mapRecv) const;
    
    void updateData(const Teuchos::RCP<FIELD>    & Field,
                    const Teuchos::RCP<SENDRECV> & mapSend) const;
    
    void updateData(FIELD    & Field,
                    SENDRECV & mapSend) const;
    
    void updateData(const Teuchos::RCP<FIELD>    & Field,
                    const Teuchos::RCP<SENDRECV> & mapSend,
                    const Teuchos::RCP<SENDRECV> & mapRecv) const;
    
    void updateData(FIELD    & Field,
                    SENDRECV & mapSend,
                    SENDRECV & mapRecv) const;
    //@}
     
    /*! @name Non-pending update functions */ //@{
  public:
    void updateDataI(const Teuchos::RCP<FIELD>    & Field,
                     const Teuchos::RCP<SENDRECV> & mapSend,
                           sVect<DOFVECT>         & bufSegments,
                           sVect<PVPS>            & sendPvps,
                           sVect<PVPS>            & recvPvps,
                     const UInt                   & channel = 2) const;
    
    void updateDataO(const Teuchos::RCP<FIELD>    & Field,
                     const Teuchos::RCP<SENDRECV> & mapSend,
                           sVect<DOFVECT>         & bufSegments,
                           sVect<PVPS>            & sendPvps,
                           sVect<PVPS>            & recvPvps,
                     const UInt                   & channel = 2) const;
     
    void updateDataI(FIELD          & Field,
                     SENDRECV       & mapSend,
                     sVect<DOFVECT> & bufSegments,
                     sVect<PVPS>    & sendPvps,
                     sVect<PVPS>    & recvPvps,
                     const UInt     & channel = 2) const;
   
    void updateDataO(FIELD          & Field,
                     SENDRECV       & mapSend,
                     sVect<DOFVECT> & bufSegments,
                     sVect<PVPS>    & sendPvps,
                     sVect<PVPS>    & recvPvps,
                     const UInt     & channel = 2) const;
     
    void updateDataI(const Teuchos::RCP<FIELD>    & Field,
                     const Teuchos::RCP<SENDRECV> & mapSend,
                     const Teuchos::RCP<SENDRECV> & mapRecv,
                           sVect<DOFVECT>         & bufSegments,
                           sVect<PVPS>            & sendPvps,
                           sVect<PVPS>            & recvPvps,
                     const UInt                   & channel = 2) const;
    
    void updateDataO(const Teuchos::RCP<FIELD>    & Field,
                     const Teuchos::RCP<SENDRECV> & mapSend,
                     const Teuchos::RCP<SENDRECV> & mapRecv,
                           sVect<DOFVECT>         & bufSegments,
                           sVect<PVPS>            & sendPvps,
                           sVect<PVPS>            & recvPvps,
                     const UInt                   & channel = 2) const;
     
    void updateDataI(FIELD          & Field,
                     SENDRECV       & mapSend,
                     SENDRECV       & mapRecv,
                     sVect<DOFVECT> & bufSegments,
                     sVect<PVPS>    & sendPvps,
                     sVect<PVPS>    & recvPvps,
                     const UInt     & channel = 2) const;
   
    void updateDataO(FIELD          & Field,
                     SENDRECV       & mapSend,
                     SENDRECV       & mapRecv,
                     sVect<DOFVECT> & bufSegments,
                     sVect<PVPS>    & sendPvps,
                     sVect<PVPS>    & recvPvps,
                     const UInt     & channel = 2) const;
    //@}
     
    /*! @name Recursive update functions */ //@{
  public:
    void updateDataRR(const Teuchos::RCP<FIELD>          & Field,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
      
    void updateDataRR(      FIELD    & Field,
                      const SENDRECV & mapSend,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
      
    void updateDataRI(const Teuchos::RCP<FIELD>          & Field,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
     
    void updateDataRI(      FIELD    & Field,
                      const SENDRECV & mapSend,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
      
    void updateDataRO(const Teuchos::RCP<FIELD>          & Field,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
     
    void updateDataRO(      FIELD    & Field,
                      const SENDRECV & mapSend,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
     
    void updateDataRR(const Teuchos::RCP<FIELD>          & Field,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                      const Teuchos::RCP<const SENDRECV> & mapRecv,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
      
    void updateDataRR(      FIELD    & Field,
                      const SENDRECV & mapSend,
                      const SENDRECV & mapRecv,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
      
    void updateDataRI(const Teuchos::RCP<FIELD>          & Field,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                      const Teuchos::RCP<const SENDRECV> & mapRecv,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
     
    void updateDataRI(      FIELD    & Field,
                      const SENDRECV & mapSend,
                      const SENDRECV & mapRecv,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
      
    void updateDataRO(const Teuchos::RCP<FIELD>          & Field,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                      const Teuchos::RCP<const SENDRECV> & mapRecv,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
     
    void updateDataRO(      FIELD    & Field,
                      const SENDRECV & mapSend,
                      const SENDRECV & mapRecv,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
    //@}
      
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const FIELD>        & OldField,
                            const Teuchos::RCP<communicator>       & NewCommDev,
                            const Teuchos::RCP<MESH1D>             & NewGrid,
                            const Teuchos::RCP<CONNECT1D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                    & isActive,
                            const communicator            & OldCommDev,
                            const FIELD                   & OldField,
                                  communicator            & NewCommDev,
                            const Teuchos::RCP<MESH1D>    & NewGrid,
                            const Teuchos::RCP<CONNECT1D> & NewConnect,
                                  FIELD                   & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const FIELD>        & OldField,
                            const Teuchos::RCP<communicator>       & NewCommDev,
                            const Teuchos::RCP<MESH1D>             & NewGrid,
                            const Teuchos::RCP<CONNECT1D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                     & isActive,
                            const communicator             & OldCommDev,
                            const FIELD                    & OldField,
                                  communicator             & NewCommDev,
                            const Teuchos::RCP<MESH1D>     & NewGrid,
                            const Teuchos::RCP<CONNECT1D>  & NewConnect,
                                  FIELD                    & NewField);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
fdf1dGlobalManip()
{
  commDevLoaded = false;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
fdf1dGlobalManip(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
fdf1dGlobalManip(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}


//_________________________________________________________________________________________________
// OPERATOR FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(const Teuchos::RCP<MESH1D>    & newGrid,
           const Teuchos::RCP<CONNECT1D> & newConnect,
           const Teuchos::RCP<OPTIONS>   & options,
           const Teuchos::RCP<FIELD>     & Field) const
{
  assert(commDevLoaded);
  assert(newGrid->getElements().colIsLocal());
  
  //Create the new field
  FIELD newField1d;
  newField1d.setCommunicator(commDev);
  newField1d.setGeometry(newGrid,newConnect);
  newField1d.setOptions(options);
  
  //Remap the feCards
  pMap<PMAPTYPE> newElMap = newGrid->getElements().getRowMap();
  FECARDS         feCards = Field->getFeCards();
  
  pVectGlobalManip<FECARD,PMAPTYPE> feCardsManip(commDev);
  feCardsManip.changeMap(feCards,newElMap);
  
  newField1d.setFeCards(feCards);
  
  //Startup
  newField1d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField1d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field->getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap));
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);
  
  //Loading the new dofVector
  newField1d.setDofVect(oldDofVect);
  *Field = newField1d;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(MESH1D    & newGrid,
           CONNECT1D & newConnect,
           OPTIONS   & options,
           FIELD     & Field) const
{
  assert(commDevLoaded);
  assert(newGrid.getElements().colIsLocal());
  
  //Create the new field
  FIELD newField1d;
  newField1d.setCommunicator(commDev);
  newField1d.setGeometry(newGrid,newConnect);
  newField1d.setOptions(options);
  
  //Remap the feCards
  pMap<PMAPTYPE> newElMap = newGrid.getElements().getRowMap();
  FECARDS         feCards = Field.getFeCards();
  
  pVectGlobalManip<FECARD,PMAPTYPE> feCardsManip(commDev);
  
  assert(feCardsManip.check(feCards));
  feCardsManip.changeMap(feCards,newElMap); 
  
  newField1d.setFeCards(feCards);  
  
  //Startup
  newField1d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField1d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field.getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap)); 
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);  

  //Loading the new dofVector
  newField1d.setDofVect(oldDofVect);
  Field = newField1d;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
bool
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(const Teuchos::RCP<MESH1D>    & newGrid1d,
          const Teuchos::RCP<CONNECT1D> & newConnect1d,
          const Teuchos::RCP<OPTIONS>   & options,
          const Teuchos::RCP<FIELD>     & field) const
{
  return( changeGrid(newGrid1d,
                     newConnect1d,
                     options,
                     field) );
}
    
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
bool
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(MESH1D    & newGrid1d,
          CONNECT1D & newConnect1d,
          OPTIONS   & options,
          FIELD     & field) const
{
  return( changeGrid(newGrid1d,
                     newConnect1d,
                     options,
                     field) );
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
createSendRecvMap(const Teuchos::RCP<FIELD>    & Field,
                        Teuchos::RCP<SENDRECV> & mapSend,
                        Teuchos::RCP<SENDRECV> & mapRecv) const
{
  assert(commDevLoaded);
  
  pMapGlobalManip<PMAPTYPE> mapManipulator(commDev);
  mapManipulator.createSendRecvMap(Field.getDofVect().getMapRcp(), mapSend, mapRecv);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
createSendRecvMap(const FIELD    & Field,
                        SENDRECV & mapSend,
                        SENDRECV & mapRecv) const
{
  assert(commDevLoaded);
  
  pMapGlobalManip<PMAPTYPE> mapManipulator(commDev);
  mapManipulator.createSendRecvMap(Field.getDofVect().getMapRef(), mapSend, mapRecv);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateData(const Teuchos::RCP<FIELD>    & Field,
           const Teuchos::RCP<SENDRECV> & mapSend) const
{
  //Checks
  assert(commDevLoaded);
    
  //Update the data
  pVect<DOFTYPE,pMapItemShare> dofsVect = Field->getDofVect();
  
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateData(dofsVect,*mapSend);
  
  //Reload the dofs
  Field.setDofVect(dofsVect);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateData(FIELD    & Field,
           SENDRECV & mapSend) const
{
  //Checks
  assert(commDevLoaded);
  
  //Update the data
  pVect<DOFTYPE,pMapItemShare> dofsVect = Field.getDofVect();
  
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateData(dofsVect,mapSend);
  
  //Reload the dofs
  Field.setDofVect(dofsVect);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateData(const Teuchos::RCP<FIELD>    & Field,
           const Teuchos::RCP<SENDRECV> & mapSend,
           const Teuchos::RCP<SENDRECV> & mapRecv) const
{
  //Checks
  assert(commDevLoaded);
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  cheker.check(mapSend,mapRecv);
  
  //Update the data
  pVect<DOFTYPE,pMapItemShare> dofsVect = Field->getDofVect();
  
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateData(dofsVect,*mapSend,*mapRecv);
  
  //Reload the dofs
  Field.setDofVect(dofsVect);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateData(FIELD    & Field,
           SENDRECV & mapSend,
           SENDRECV & mapRecv) const
{
  //Checks
  assert(commDevLoaded);
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  cheker.check(mapSend,mapRecv);
  
  //Update the data
  pVect<DOFTYPE,pMapItemShare> dofsVect = Field.getDofVect();
  
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateData(dofsVect,mapSend,mapRecv);
  
  //Reload the dofs
  Field.setDofVect(dofsVect);
}



//_________________________________________________________________________________________________
// NON-PENDING UPDATE FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataI(const Teuchos::RCP<FIELD>    & Field,
            const Teuchos::RCP<SENDRECV> & mapSend,
                  sVect<DOFVECT>         & bufSegments,
                  sVect<PVPS>            & sendPvps,
                  sVect<PVPS>            & recvPvps,
            const UInt                   & channel) const
{
  assert(commDevLoaded);
  updateDataI(*Field,
              *mapSend,
               bufSegments,
               sendPvps,
               recvPvps,
               channel);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataO(const Teuchos::RCP<FIELD>    & Field,
            const Teuchos::RCP<SENDRECV> & mapSend,
                  sVect<DOFVECT>         & bufSegments,
                  sVect<PVPS>            & sendPvps,
                  sVect<PVPS>            & recvPvps,
            const UInt                   & channel) const
{
  assert(commDevLoaded);
  updateDataO(*Field,
              *mapSend,
               bufSegments,
               sendPvps,
               recvPvps,
               channel);
}
   
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>     
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataI(FIELD          & Field,
            SENDRECV       & mapSend,
            sVect<DOFVECT> & bufSegments,
            sVect<PVPS>    & sendPvps,
            sVect<PVPS>    & recvPvps,
	    const UInt     & channel) const
{
  //Checks
  assert(commDevLoaded);
  
  //Dofs Vect
  pVect<DOFTYPE,pMapItemShare> & dofsVect = Field.dofVect;
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataI(dofsVect,
                          mapSend,
                          bufSegments,
                          sendPvps,
                          recvPvps,
                          channel);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataO(FIELD          & Field,
            SENDRECV       & mapSend,
            sVect<DOFVECT> & bufSegments,
            sVect<PVPS>    & sendPvps,
            sVect<PVPS>    & recvPvps,
            const UInt     & channel) const
{
  //Checks
  assert(commDevLoaded);
  
  //Dofs Vect
  pVect<DOFTYPE,pMapItemShare> & dofsVect = Field.dofVect;
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataO(dofsVect,
                          mapSend,
                          bufSegments,
                          sendPvps,
                          recvPvps,
                          channel);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataI(const Teuchos::RCP<FIELD>    & Field,
            const Teuchos::RCP<SENDRECV> & mapSend,
            const Teuchos::RCP<SENDRECV> & mapRecv,
                  sVect<DOFVECT>         & bufSegments,
                  sVect<PVPS>            & sendPvps,
                  sVect<PVPS>            & recvPvps,
            const UInt                   & channel) const
{
  assert(commDevLoaded);
  updateDataI(*Field,
              *mapSend,
              *mapRecv,
               bufSegments,
               sendPvps,
               recvPvps,
               channel);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataO(const Teuchos::RCP<FIELD>    & Field,
            const Teuchos::RCP<SENDRECV> & mapSend,
            const Teuchos::RCP<SENDRECV> & mapRecv,
                  sVect<DOFVECT>         & bufSegments,
                  sVect<PVPS>            & sendPvps,
                  sVect<PVPS>            & recvPvps,
            const UInt                   & channel) const
{
  assert(commDevLoaded);
  updateDataO(*Field,
              *mapSend,
              *mapRecv,
               bufSegments,
               sendPvps,
               recvPvps,
               channel);
}
   
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>     
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataI(FIELD          & Field,
            SENDRECV       & mapSend,
            SENDRECV       & mapRecv,
            sVect<DOFVECT> & bufSegments,
            sVect<PVPS>    & sendPvps,
            sVect<PVPS>    & recvPvps,
            const UInt     & channel) const
{
  //Checks
  assert(commDevLoaded);
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  cheker.check(mapSend,mapRecv);
  
  //Dofs Vect
  pVect<DOFTYPE,pMapItemShare> & dofsVect = Field.dofVect;
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataI(dofsVect,
                          mapSend,
                          mapRecv,
                          bufSegments,
                          sendPvps,
                          recvPvps, 
                          channel);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataO(FIELD          & Field,
            SENDRECV       & mapSend,
            SENDRECV       & mapRecv,
            sVect<DOFVECT> & bufSegments,
            sVect<PVPS>    & sendPvps,
            sVect<PVPS>    & recvPvps,
            const UInt     & channel) const
{
  //Checks
  assert(commDevLoaded);
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  cheker.check(mapSend,mapRecv);
  
  //Dofs Vect
  pVect<DOFTYPE,pMapItemShare> & dofsVect = Field.dofVect;
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataO(dofsVect,
                          mapSend,
                          mapRecv,
                          bufSegments,
                          sendPvps,
                          recvPvps,
                          channel);
}


//_________________________________________________________________________________________________
// RECURSIVE UPDATE FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRR(const Teuchos::RCP<FIELD>          & Field,
             const Teuchos::RCP<const SENDRECV> & mapSend,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRR(*Field,
               *mapSend,
                commBuffer,
                channel);
}
      
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRR(      FIELD    & Field,
             const SENDRECV & mapSend,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataRR(Field.dofVect,
                           mapSend,
                           commBuffer,
                           channel);
}
      
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRI(const Teuchos::RCP<FIELD>          & Field,
             const Teuchos::RCP<const SENDRECV> & mapSend,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRI(*Field,
               *mapSend,
                commBuffer,
                channel);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRI(      FIELD    & Field,
             const SENDRECV & mapSend,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataRI(Field.dofVect,
                           mapSend,
                           commBuffer,
                           channel);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRO(const Teuchos::RCP<FIELD>          & Field,
             const Teuchos::RCP<const SENDRECV> & mapSend,
                    PVUR & commBuffer,
              const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRO(*Field,
               *mapSend,
                commBuffer,
                channel);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRO(      FIELD    & Field,
             const SENDRECV & mapSend,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataRO(Field.dofVect,
                           mapSend,
                           commBuffer,
                           channel);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRR(const Teuchos::RCP<FIELD>          & Field,
             const Teuchos::RCP<const SENDRECV> & mapSend,
             const Teuchos::RCP<const SENDRECV> & mapRecv,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRR(*Field,
               *mapSend,
               *mapRecv,
                commBuffer,
                channel);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRR(      FIELD    & Field,
             const SENDRECV & mapSend,
             const SENDRECV & mapRecv,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  
  assert(commDevLoaded);
  assert(cheker.check(mapSend,mapRecv));
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataRR(Field.dofVect,
                           mapSend,
                           mapRecv,
                           commBuffer,
                           channel);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRI(const Teuchos::RCP<FIELD>          & Field,
             const Teuchos::RCP<const SENDRECV> & mapSend,
             const Teuchos::RCP<const SENDRECV> & mapRecv,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRI(*Field,
               *mapSend,
               *mapRecv,
                commBuffer,
                channel);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRI(      FIELD    & Field,
             const SENDRECV & mapSend,
             const SENDRECV & mapRecv,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  
  assert(commDevLoaded);
  assert(cheker.check(mapSend,mapRecv));
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataRI(Field.dofVect,
                           mapSend,
                           mapRecv,
                           commBuffer,
                           channel);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRO(const Teuchos::RCP<FIELD>          & Field,
             const Teuchos::RCP<const SENDRECV> & mapSend,
             const Teuchos::RCP<const SENDRECV> & mapRecv,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRO(*Field,
               *mapSend,
               *mapRecv,
                commBuffer,
                channel);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>  
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRO(      FIELD    & Field,
             const SENDRECV & mapSend,
             const SENDRECV & mapRecv,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  
  assert(commDevLoaded);
  assert(cheker.check(mapSend,mapRecv));
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataRO(Field.dofVect,
                           mapSend,
                           mapRecv,
                           commBuffer,
                           channel);
}


//_________________________________________________________________________________________________
// COMMUNICATOR MANIPULATIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>  
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH1D>             & NewGrid,
                   const Teuchos::RCP<CONNECT1D>          & NewConnect,
                         Teuchos::RCP<FIELD>              & NewField)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldField,
                    *NewCommDev,
                     NewGrid,
                     NewConnect,
                    *NewField);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>  
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH1D>    & NewGrid,
                   const Teuchos::RCP<CONNECT1D> & NewConnect,
                         FIELD                   & NewField)
{
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::DOFVECT DOFVECT;
  typedef typename FIELD::FECARDS FECARDS;
  typedef typename FIELD::FECARD  FECARD;
  
  typedef pVectGlobalManip<DOFTYPE,PMAPTYPE>  MANIP_DOFVECT;
  typedef pVectGlobalManip<FECARD,PMAPTYPE>   MANIP_FECARDS;
  
  //Manipulators-----------------------------------------------------
  DOFVECT newDofVect;
  FECARDS newFeCards;

  MANIP_DOFVECT manipDofVect;
  MANIP_FECARDS manipFeCards;
  
  OPTIONS options = OldField.getOptions();
  
  //Copy to new communicator-----------------------------------------
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() < OldCommDev.size());
    
    //Reduce comm
    manipDofVect.reduceCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getDofVect(),
                                    NewCommDev,
                                    newDofVect);
    
    manipFeCards.reduceCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getFeCards(),
                                    NewCommDev,
                                    newFeCards);
    
    //New Field
    NewField.setCommunicator(NewCommDev);
    NewField.setGeometry(NewGrid,NewConnect);
    NewField.setOptions(options);
    NewField.setFeCards(newFeCards);
    NewField.startup();
    NewField.setDofVect(newDofVect);
  }
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>  
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH1D>             & NewGrid,
                   const Teuchos::RCP<CONNECT1D>          & NewConnect,
                         Teuchos::RCP<FIELD>              & NewField)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldField,
                    *NewCommDev,
                     NewGrid,
                     NewConnect,
                    *NewField);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>  
void
fdf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH1D>    & NewGrid,
                   const Teuchos::RCP<CONNECT1D> & NewConnect,
                         FIELD                   & NewField)
{
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::DOFVECT DOFVECT;
  typedef typename FIELD::FECARDS FECARDS;
  typedef typename FIELD::FECARD  FECARD;
  
  typedef pVectGlobalManip<DOFTYPE,DOFTYPE>  MANIP_DOFVECT;
  typedef pVectGlobalManip<FECARD,DOFTYPE>   MANIP_FECARDS;
  
  //Manipulators-----------------------------------------------------
  DOFVECT newDofVect;
  FECARDS newFeCards;

  MANIP_DOFVECT manipDofVect;
  MANIP_FECARDS manipFeCards;
  
  OPTIONS options = OldField.getOptions();
  
  //Copy to new communicator-----------------------------------------
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() > OldCommDev.size());
    
    //Reduce comm
    manipDofVect.expandCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getDofVect(),
                                    NewCommDev,
                                    newDofVect);
    
    manipFeCards.expandCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getFeCards(),
                                    NewCommDev,
                                    newFeCards);
    
    //New Field
    NewField.setCommunicator(NewCommDev);
    NewField.setGeometry(NewGrid,NewConnect);
    NewField.setOptions(options);
    NewField.setFeCards(newFeCards);
    NewField.startup();
    NewField.setDofVect(newDofVect);
  }
}



//_________________________________________________________________________________________________
// SUMMARY CLASS
//-------------------------------------------------------------------------------------------------

/*! Interface class for the manipulation of Dynamic Field 1d */
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER = dmd1d_vectMajor, dmd1d_mode MODE = dmd1d_standard>
class feDynamicField1dGlobalManip : public fdf1dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>
{
     /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE PMAPTYPE;
    typedef feDynamicField1d<FETYPE,DOFTYPE,ORDER,MODE> FIELD;
    typedef typename FETYPE::GEOSHAPE  GEOSHAPE;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    feDynamicField1dGlobalManip();
    feDynamicField1dGlobalManip(const Teuchos::RCP<communicator> & CommDev);
    feDynamicField1dGlobalManip(communicator & CommDev);
    //@}
};


template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
feDynamicField1dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>::
feDynamicField1dGlobalManip() : fdf1dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>()
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
feDynamicField1dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>::
feDynamicField1dGlobalManip(const Teuchos::RCP<communicator> & CommDev) : fdf1dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>(CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
feDynamicField1dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>::
feDynamicField1dGlobalManip(communicator & CommDev) : fdf1dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>(CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
}

#endif
