/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FEDYNAMICFIELD2DGLOBALMANIP_HPP
#define FEDYNAMICFIELD2DGLOBALMANIP_HPP

#include "pMapItem.h"
#include "pMapItemShare.h"

#include "ggmGridSearch.hpp"
#include "feDynamicField2d.hpp"

using namespace std;


//_________________________________________________________________________________________________
// SERVO-SPECIALIZATION
//-------------------------------------------------------------------------------------------------

/*! Dynamic Field 2d Global Manipulator - empty general class  */
template<typename PMAPTYPE, typename FETYPE, typename DOFTYPE, dmd2d_order ORDER = dmd2d_vectMajor, dmd2d_mode MODE = dmd2d_standard>
class fdf2dGlobalManip;



//_________________________________________________________________________________________________
// PMAPITEM
//-------------------------------------------------------------------------------------------------

/*! Dynamic Field 2d Global Manipulator - \c pMapItem specialization */
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
class fdf2dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItem                                    PMAPTYPE;
    typedef dofMapDynamic2d_options                     OPTIONS;
    typedef feDynamicField2d<FETYPE,DOFTYPE,ORDER,MODE> FIELD;
    
    typedef typename FIELD::MESH2D     MESH2D;
    typedef typename FIELD::CONNECT2D  CONNECT2D;
    typedef typename FIELD::DOFVECT    DOFVECT;
    typedef typename FIELD::GEOSHAPE   GEOSHAPE;
    typedef typename FETYPE::ELCARD    ELCARD;
    typedef typename FETYPE::FECARD    FECARD;
    typedef pVect<FECARD,PMAPTYPE>     FECARDS;
    
    typedef typename GEOSHAPE::GEOBSHAPE            GEOSHAPE1D;
    typedef mesh1d<GEOSHAPE1D,PMAPTYPE,PMAPTYPE>    MESH1D;
    typedef connect1d<GEOSHAPE1D,PMAPTYPE,PMAPTYPE> CONNECT1D;
    //@}
  
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    fdf2dGlobalManip();
    fdf2dGlobalManip(const Teuchos::RCP<communicator> & CommDev);
    fdf2dGlobalManip(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Operator functions */ //@{
  public:
    void changeGrid(const Teuchos::RCP<MESH2D>    & newGrid,
                    const Teuchos::RCP<CONNECT2D> & newConnect,
                    const Teuchos::RCP<OPTIONS>   & options,
                    const Teuchos::RCP<FIELD>     & Field) const;
    
    void changeGrid(MESH2D    & newGrid,
                    CONNECT2D & newConnect,
                    OPTIONS   & options,
                    FIELD     & Field) const;
    
    bool matchGrid(const Teuchos::RCP<MESH2D>    & newGrid2d,
                   const Teuchos::RCP<CONNECT2D> & newConnect2d,
                   const Teuchos::RCP<OPTIONS>   & options,
                   const Teuchos::RCP<FIELD>     & field) const;
    
    bool matchGrid(MESH2D    & newGrid2d,
                   CONNECT2D & newConnect2d,
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
                            const Teuchos::RCP<MESH2D>             & NewGrid,
                            const Teuchos::RCP<CONNECT2D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                    & isActive,
                            const communicator            & OldCommDev,
                            const FIELD                   & OldField,
                                  communicator            & NewCommDev,
                            const Teuchos::RCP<MESH2D>    & NewGrid,
                            const Teuchos::RCP<CONNECT2D> & NewConnect,
                                  FIELD                   & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const FIELD>        & OldField,
                            const Teuchos::RCP<communicator>       & NewCommDev,
                            const Teuchos::RCP<MESH2D>             & NewGrid,
                            const Teuchos::RCP<CONNECT2D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                    & isActive,
                            const communicator            & OldCommDev,
                            const FIELD                   & OldField,
                                  communicator            & NewCommDev,
                            const Teuchos::RCP<MESH2D>    & NewGrid,
                            const Teuchos::RCP<CONNECT2D> & NewConnect,
                                  FIELD                   & NewField);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
fdf2dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
fdf2dGlobalManip()
{
  commDevLoaded = false;
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
fdf2dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
fdf2dGlobalManip(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
fdf2dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
fdf2dGlobalManip(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}


//_________________________________________________________________________________________________
// OPERATOR FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(const Teuchos::RCP<MESH2D>    & newGrid,
           const Teuchos::RCP<CONNECT2D> & newConnect,
           const Teuchos::RCP<OPTIONS>   & options,
           const Teuchos::RCP<FIELD>     & Field) const
{
  assert(commDevLoaded);
  assert(newGrid->getElements().colIsLocal());
  
  //Create the new field
  FIELD newField2d;
  newField2d.setCommunicator(commDev);
  newField2d.setGeometry(newGrid,newConnect);
  newField2d.setOptions(options);
  
  //Remap the feCards
  pMap<PMAPTYPE> newElMap = newGrid->getElements().getRowMap();
  FECARDS         feCards = Field->getFeCards();
  
  pVectGlobalManip<FECARD,PMAPTYPE> feCardsManip(commDev);
  feCardsManip.changeMap(feCards,newElMap);
  
  newField2d.setFeCards(feCards);
  
  //Startup
  newField2d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField2d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field->getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap));
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);
  
  //Loading the new dofVector
  newField2d.setDofVect(oldDofVect);
  *Field = newField2d;
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(MESH2D    & newGrid,
           CONNECT2D & newConnect,
           OPTIONS   & options,
           FIELD     & Field) const
{
  assert(commDevLoaded);
  assert(newGrid->getElements().colIsLocal());
  
  //Create the new field
  FIELD newField2d;
  newField2d.setCommunicator(commDev);
  newField2d.setGeometry(newGrid,newConnect);
  newField2d.setOptions(options);
  
  //Remap the feCards
  pMap<PMAPTYPE> newElMap = newGrid.getElements().getRowMap();
  FECARDS         feCards = Field.getFeCards();
  
  pVectGlobalManip<FECARD,PMAPTYPE> feCardsManip(commDev);
  
  assert(feCardsManip.check(feCards));
  feCardsManip.changeMap(feCards,newElMap); 
  
  newField2d.setFeCards(feCards);  
  
  //Startup
  newField2d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField2d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field.getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap)); 
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);  

  //Loading the new dofVector
  newField2d.setDofVect(oldDofVect);
  Field = newField2d;
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
bool
fdf2dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(const Teuchos::RCP<MESH2D>    & newGrid2d,
          const Teuchos::RCP<CONNECT2D> & newConnect2d,
          const Teuchos::RCP<OPTIONS>   & options,
          const Teuchos::RCP<FIELD>     & field) const
{
  assert(commDevLoaded);
  
  matchGrid(*newGrid2d,
            *newConnect2d,
            *options,
            *field);
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
bool
fdf2dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(MESH2D    & newGrid2d,
          CONNECT2D & newConnect2d,
          OPTIONS   & options,
          FIELD     & field) const
{
  //Typdefs----------------------------------------------------------------------------------------
  typedef ggmGridSearch<MESH1D,MESH1D>  EDGESEARCH;
  typedef typename EDGESEARCH::OUTVECT  EDGEVECT;
  
  typedef pMap<PMAPTYPE> EDGEMAP;
  
  typedef typename FIELD::OPTIONS FIELD_OPTIONS;
  
  //Assert-----------------------------------------------------------------------------------------
  assert(commDevLoaded);
  
  //Download data----------------------------------------------------------------------------------
  Teuchos::RCP<MESH2D>    oldGrid2d(new MESH2D(*field.getMesh2d()));
  Teuchos::RCP<CONNECT2D> oldConnect2d(new CONNECT2D(*field.getMeshConnect()));
  
  //Alloc------------------------------------------------------------------------------------------
  bool success = true;
  sVect<point3d> inLocCoords;
  
  //Map the edges----------------------------------------------------------------------------------
  Teuchos::RCP<MESH1D> oldEdgeGrid1d(new MESH1D(oldGrid2d->getNodes(),
                                                oldGrid2d->getIsVertex(),
                                                oldGrid2d->getEdges()));
  
  Teuchos::RCP<MESH1D> newEdgeGrid1d(new MESH1D(newGrid2d.getNodes(),
                                                newGrid2d.getIsVertex(),
                                                newGrid2d.getEdges()));
  
  Teuchos::RCP<CONNECT1D> oldEdgeConnect1d(new CONNECT1D(commDev));
  oldEdgeConnect1d->setMesh1d(oldEdgeGrid1d);
  
  Teuchos::RCP<CONNECT1D> newEdgeConnect1d(new CONNECT1D(commDev));
  newEdgeConnect1d->setMesh1d(newEdgeGrid1d);
  
  EDGESEARCH edgeMatcher(commDev);
  edgeMatcher.setMesh(oldEdgeGrid1d, //tgt
                      oldEdgeConnect1d,
                      newEdgeGrid1d, //src
                      newEdgeConnect1d);
  
  EDGEVECT edgeVect = edgeMatcher.search(inLocCoords);
  EDGEMAP  edgeMap(edgeVect.size());
  
  for(UInt i=1; i <= edgeVect.size(); ++i)
  {
    success   = success && edgeVect(i).getIsNested();
    edgeMap(i) = edgeVect(i).getElMap();
  }
  
  if(!success) { return(success); }
  
  assert(oldGrid2d->getNumEdges() == edgeMap.size());
  
  //Change the gids--------------------------------------------------------------------------------  
  for(UInt d=1; d <= oldGrid2d->getNumEdges(); ++d)
  { oldGrid2d->getEdges().getRowMapL(d).setGid(edgeMap(d).getGid()); }
  
  oldConnect2d->getVertexToEdge().setColMap(oldGrid2d->getEdges().getRowMap());
  oldConnect2d->getElementToEdge().setColMap(oldGrid2d->getEdges().getRowMap());
  
  //Build old field--------------------------------------------------------------------------------
  Teuchos::RCP<FIELD_OPTIONS> oldOptions(new FIELD_OPTIONS(field.getOptions())); 
  
  Teuchos::RCP<FIELD> oldField(new FIELD);
  oldField->setCommunicator(commDev);
  oldField->setGeometry(oldGrid2d,oldConnect2d);
  oldField->setFeCards(field.getFeCards());
  oldField->setOptions(oldOptions);
  oldField->startup();
  
  for(UInt lid=1; lid <= oldField->getDofVect().size(); ++lid)
  { oldField->setDofL(lid, field.getDofL(lid)); }
  
  //Change the maps--------------------------------------------------------------------------------
  changeGrid(newGrid2d,
             newConnect2d,
             options,
            *oldField);
  
  field = *oldField;
  
  return(success);
}


//_________________________________________________________________________________________________
// COMMUNICATOR MANIPULATIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH2D>             & NewGrid,
                   const Teuchos::RCP<CONNECT2D>          & NewConnect,
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH2D>    & NewGrid,
                   const Teuchos::RCP<CONNECT2D> & NewConnect,
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH2D>             & NewGrid,
                   const Teuchos::RCP<CONNECT2D>          & NewConnect,
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH2D>    & NewGrid,
                   const Teuchos::RCP<CONNECT2D> & NewConnect,
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
// PMAPITEMSHARE
//-------------------------------------------------------------------------------------------------

/*! Dynamic Field 2d Global Manipulator - \c pMapItemShare specialization */
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
class fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItemShare                               PMAPTYPE;
    typedef feDynamicField2d<FETYPE,DOFTYPE,ORDER,MODE> FIELD;
    typedef pMap<pMapItemSendRecv>                      SENDRECV;
    typedef dofMapDynamic2d_options                     OPTIONS;
    
    typedef typename FIELD::MESH2D     MESH2D;
    typedef typename FIELD::CONNECT2D  CONNECT2D;
    typedef typename FIELD::DOFVECT    DOFVECT;
    typedef typename FIELD::GEOSHAPE   GEOSHAPE;
    typedef typename FETYPE::ELCARD    ELCARD;
    typedef typename FETYPE::FECARD    FECARD;
    
    typedef pVect<FECARD,PMAPTYPE>                  FECARDS;
    typedef pVectGlobalManip<DOFTYPE,pMapItemShare> PVGLOBMANIP;
    typedef typename PVGLOBMANIP::PVPS              PVPS;
    typedef typename PVGLOBMANIP::PVUR              PVUR;
    
    typedef typename GEOSHAPE::GEOBSHAPE            GEOSHAPE1D;
    typedef mesh1d<GEOSHAPE1D,PMAPTYPE,PMAPTYPE>    MESH1D;
    typedef connect1d<GEOSHAPE1D,PMAPTYPE,PMAPTYPE> CONNECT1D;
    //@}
  
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    fdf2dGlobalManip();
    fdf2dGlobalManip(const Teuchos::RCP<communicator> & CommDev);
    fdf2dGlobalManip(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Operator functions */ //@{
  public:
    void changeGrid(const Teuchos::RCP<MESH2D>    & newGrid,
                    const Teuchos::RCP<CONNECT2D> & newConnect,
                    const Teuchos::RCP<OPTIONS>   & options,
                    const Teuchos::RCP<FIELD>     & Field) const;
    
    void changeGrid(MESH2D    & newGrid,
                    CONNECT2D & newConnect,
                    OPTIONS   & options,
                    FIELD     & Field) const;
   
    bool matchGrid(const Teuchos::RCP<MESH2D>    & newGrid2d,
                   const Teuchos::RCP<CONNECT2D> & newConnect2d,
                   const Teuchos::RCP<OPTIONS>   & options,
                   const Teuchos::RCP<FIELD>     & field) const;
    
    bool matchGrid(MESH2D    & newGrid2d,
                   CONNECT2D & newConnect2d,
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
                            const Teuchos::RCP<MESH2D>             & NewGrid,
                            const Teuchos::RCP<CONNECT2D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                    & isActive,
                            const communicator            & OldCommDev,
                            const FIELD                   & OldField,
                                  communicator            & NewCommDev,
                            const Teuchos::RCP<MESH2D>    & NewGrid,
                            const Teuchos::RCP<CONNECT2D> & NewConnect,
                                  FIELD                   & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const FIELD>        & OldField,
                            const Teuchos::RCP<communicator>       & NewCommDev,
                            const Teuchos::RCP<MESH2D>             & NewGrid,
                            const Teuchos::RCP<CONNECT2D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                     & isActive,
                            const communicator             & OldCommDev,
                            const FIELD                    & OldField,
                                  communicator             & NewCommDev,
                            const Teuchos::RCP<MESH2D>     & NewGrid,
                            const Teuchos::RCP<CONNECT2D>  & NewConnect,
                                  FIELD                    & NewField);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
fdf2dGlobalManip()
{
  commDevLoaded = false;
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
fdf2dGlobalManip(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
fdf2dGlobalManip(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}


//_________________________________________________________________________________________________
// OPERATOR FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(const Teuchos::RCP<MESH2D>    & newGrid,
           const Teuchos::RCP<CONNECT2D> & newConnect,
           const Teuchos::RCP<OPTIONS>   & options,
           const Teuchos::RCP<FIELD>     & Field) const
{
  assert(commDevLoaded);
  
  //Create the new field
  FIELD newField2d;
  newField2d.setCommunicator(commDev);
  newField2d.setGeometry(newGrid,newConnect);
  newField2d.setOptions(options);
  
  //Remap the feCards
  pMap<PMAPTYPE> newElMap = newGrid->getElements().getRowMap();
  FECARDS         feCards = Field->getFeCards();
  
  pVectGlobalManip<FECARD,PMAPTYPE> feCardsManip(commDev);
  feCardsManip.changeMap(feCards,newElMap);
  
  newField2d.setFeCards(feCards);
  
  //Startup
  newField2d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField2d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field->getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap));
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);
  
  //Loading the new dofVector
  newField2d.setDofVect(oldDofVect);
  *Field = newField2d;
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(MESH2D    & newGrid,
           CONNECT2D & newConnect,
           OPTIONS   & options,
           FIELD     & Field) const
{
  assert(commDevLoaded);
  
  //Create the new field
  FIELD newField2d;
  newField2d.setCommunicator(commDev);
  newField2d.setGeometry(newGrid,newConnect);
  newField2d.setOptions(options);
  
  //Remap the feCards
  pMap<PMAPTYPE> newElMap = newGrid.getElements().getRowMap();
  FECARDS         feCards = Field.getFeCards();
  
  pVectGlobalManip<FECARD,PMAPTYPE> feCardsManip(commDev);
  
  assert(feCardsManip.check(feCards));
  feCardsManip.changeMap(feCards,newElMap); 
  
  newField2d.setFeCards(feCards);  
  
  //Startup
  newField2d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField2d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field.getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap)); 
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);  

  //Loading the new dofVector
  newField2d.setDofVect(oldDofVect);
  Field = newField2d;
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
bool
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(const Teuchos::RCP<MESH2D>    & newGrid2d,
          const Teuchos::RCP<CONNECT2D> & newConnect2d,
          const Teuchos::RCP<OPTIONS>   & options,
          const Teuchos::RCP<FIELD>     & field) const
{
  assert(commDevLoaded);
  
  matchGrid(*newGrid2d,
            *newConnect2d,
            *options,
            *field);
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
bool
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(MESH2D    & newGrid2d,
          CONNECT2D & newConnect2d,
          OPTIONS   & options,
          FIELD     & field) const
{
  //Typdefs----------------------------------------------------------------------------------------
  typedef ggmGridSearch<MESH1D,MESH1D>  EDGESEARCH;
  typedef typename EDGESEARCH::OUTVECT  EDGEVECT;
  
  typedef pMap<PMAPTYPE> EDGEMAP;
  
  typedef typename FIELD::OPTIONS FIELD_OPTIONS;
  
  //Assert-----------------------------------------------------------------------------------------
  assert(commDevLoaded);
  
  //Download data----------------------------------------------------------------------------------
  Teuchos::RCP<MESH2D>    oldGrid2d(new MESH2D(*field.getMesh2d()));
  Teuchos::RCP<CONNECT2D> oldConnect2d(new CONNECT2D(*field.getMeshConnect()));
  
  //Alloc------------------------------------------------------------------------------------------
  bool success = true;
  sVect<point3d> inLocCoords;
  
  //Map the edges----------------------------------------------------------------------------------
  Teuchos::RCP<MESH1D> oldEdgeGrid1d(new MESH1D(oldGrid2d->getNodes(),
                                                oldGrid2d->getIsVertex(),
                                                oldGrid2d->getEdges()));
  
  Teuchos::RCP<MESH1D> newEdgeGrid1d(new MESH1D(newGrid2d.getNodes(),
                                                newGrid2d.getIsVertex(),
                                                newGrid2d.getEdges()));
  
  Teuchos::RCP<CONNECT1D> oldEdgeConnect1d(new CONNECT1D(commDev));
  oldEdgeConnect1d->setMesh1d(oldEdgeGrid1d);
  
  Teuchos::RCP<CONNECT1D> newEdgeConnect1d(new CONNECT1D(commDev));
  newEdgeConnect1d->setMesh1d(newEdgeGrid1d);
  
  EDGESEARCH edgeMatcher(commDev);
  edgeMatcher.setMesh(oldEdgeGrid1d, //tgt
                      oldEdgeConnect1d,
                      newEdgeGrid1d, //src
                      newEdgeConnect1d);
  
  EDGEVECT edgeVect = edgeMatcher.search(inLocCoords);
  EDGEMAP  edgeMap(edgeVect.size());
  
  for(UInt i=1; i <= edgeVect.size(); ++i)
  {
    success   = success && edgeVect(i).getIsNested();
    edgeMap(i) = edgeVect(i).getElMap();
  }
  
  if(!success) { return(success); }
  
  assert(oldGrid2d->getNumEdges() == edgeMap.size());
  
  //Change the gids--------------------------------------------------------------------------------  
  for(UInt d=1; d <= oldGrid2d->getNumEdges(); ++d)
  { oldGrid2d->getEdges().getRowMapL(d).setGid(edgeMap(d).getGid()); }
  
  oldConnect2d->getVertexToEdge().setColMap(oldGrid2d->getEdges().getRowMap());
  oldConnect2d->getElementToEdge().setColMap(oldGrid2d->getEdges().getRowMap());
  
  //Build old field--------------------------------------------------------------------------------
  Teuchos::RCP<FIELD_OPTIONS> oldOptions(new FIELD_OPTIONS(field.getOptions())); 
  
  Teuchos::RCP<FIELD> oldField(new FIELD);
  oldField->setCommunicator(commDev);
  oldField->setGeometry(oldGrid2d,oldConnect2d);
  oldField->setFeCards(field.getFeCards());
  oldField->setOptions(oldOptions);
  oldField->startup();
  
  for(UInt lid=1; lid <= oldField->getDofVect().size(); ++lid)
  { oldField->setDofL(lid, field.getDofL(lid)); }
  
  //Change the maps--------------------------------------------------------------------------------
  changeGrid(newGrid2d,
             newConnect2d,
             options,
            *oldField);
  
  field = *oldField;
  
  return(success);
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
createSendRecvMap(const Teuchos::RCP<FIELD>    & Field,
                        Teuchos::RCP<SENDRECV> & mapSend,
                        Teuchos::RCP<SENDRECV> & mapRecv) const
{
  assert(commDevLoaded);
  
  pMapGlobalManip<PMAPTYPE> mapManipulator(commDev);
  mapManipulator.createSendRecvMap(Field.getDofVect().getMapRcp(), mapSend, mapRecv);
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
createSendRecvMap(const FIELD    & Field,
                        SENDRECV & mapSend,
                        SENDRECV & mapRecv) const
{
  assert(commDevLoaded);
  
  pMapGlobalManip<PMAPTYPE> mapManipulator(commDev);
  mapManipulator.createSendRecvMap(Field.getDofVect().getMapRef(), mapSend, mapRecv);
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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
   
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>     
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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
   
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>     
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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
      
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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
      
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>  
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>  
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH2D>             & NewGrid,
                   const Teuchos::RCP<CONNECT2D>          & NewConnect,
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>  
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH2D>    & NewGrid,
                   const Teuchos::RCP<CONNECT2D> & NewConnect,
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>  
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH2D>             & NewGrid,
                   const Teuchos::RCP<CONNECT2D>          & NewConnect,
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

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>  
void
fdf2dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH2D>    & NewGrid,
                   const Teuchos::RCP<CONNECT2D> & NewConnect,
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

/*! Interface class for the manipulation of Dynamic Field 2d */
template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER = dmd2d_vectMajor, dmd2d_mode MODE = dmd2d_standard>
class feDynamicField2dGlobalManip : public fdf2dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>
{
     /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE PMAPTYPE;
    typedef feDynamicField2d<FETYPE,DOFTYPE,ORDER,MODE> FIELD;
    typedef typename FETYPE::GEOSHAPE  GEOSHAPE;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    feDynamicField2dGlobalManip();
    feDynamicField2dGlobalManip(const Teuchos::RCP<communicator> & CommDev);
    feDynamicField2dGlobalManip(communicator & CommDev);
    //@}
};


template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
feDynamicField2dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>::
feDynamicField2dGlobalManip() : fdf2dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>()
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
feDynamicField2dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>::
feDynamicField2dGlobalManip(const Teuchos::RCP<communicator> & CommDev) : fdf2dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>(CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
}

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
feDynamicField2dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>::
feDynamicField2dGlobalManip(communicator & CommDev) : fdf2dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>(CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
}

#endif
