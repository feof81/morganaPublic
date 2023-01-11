/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTERPOLATORMASS_HPP
#define INTERPOLATORMASS_HPP

#include "morganaFields.hpp"
#include "traitsInterpolator.hpp"

#include "fnEL3dB.hpp"
#include "localVectorEL3d.hpp"
#include "vectorBuilder.hpp"
#include "epetraVector_to_pVect.h"
#include "traitsColMapFixer.hpp"


/*! Mass preserving interpolator. The \TARGETFIELD must be Lagrangian, must statisfy the partition of the unity. */
template<typename SOURCEFIELD, typename TARGETFIELD>
class interpolatorMass : public interpolatorTrait<TARGETFIELD>::SEARCH
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename SOURCEFIELD::GEOSHAPE      SOURCE_GEOSHAPE;
    typedef typename SOURCEFIELD::OUTTYPE       SOURCE_OUTTYPE;
    typedef typename SOURCEFIELD::FIELD_FETYPE  SOURCE_FETYPE;
    typedef typename SOURCEFIELD::PMAPTYPE      SOURCE_PMAPTYPE;
    typedef typename SOURCEFIELD::FIELD_DOFTYPE SOURCE_DOFTYPE;
    typedef typename SOURCEFIELD::FEINTERFACE   SOURCE_FEINTERFACE;
    typedef typename interpolatorTrait<SOURCEFIELD>::MESH    SOURCE_MESH;
    typedef typename interpolatorTrait<SOURCEFIELD>::CONNECT SOURCE_CONNECT;
    typedef typename interpolatorTrait<SOURCEFIELD>::SUPPORT SOURCE_SUPPORT;
    typedef typename interpolatorTrait<SOURCEFIELD>::MANIP   SOURCE_MANIP;
    
    typedef typename TARGETFIELD::GEOSHAPE      TARGET_GEOSHAPE;
    typedef typename TARGETFIELD::OUTTYPE       TARGET_OUTTYPE;
    typedef typename TARGETFIELD::FIELD_FETYPE  TARGET_FETYPE;
    typedef typename TARGETFIELD::PMAPTYPE      TARGET_PMAPTYPE;
    typedef typename TARGETFIELD::FIELD_DOFTYPE TARGET_DOFTYPE;
    typedef typename TARGETFIELD::FEINTERFACE   TARGET_FEINTERFACE;
    typedef typename interpolatorTrait<TARGETFIELD>::MESH    TARGET_MESH;
    typedef typename interpolatorTrait<TARGETFIELD>::CONNECT TARGET_CONNECT;
    typedef typename interpolatorTrait<TARGETFIELD>::SUPPORT TARGET_SUPPORT;
    typedef typename interpolatorTrait<TARGETFIELD>::MANIP   TARGET_MANIP;
    
    typedef typename interpolatorTrait<TARGETFIELD>::SEARCH SEARCH;
    typedef searchData<SOURCE_PMAPTYPE>                     SEARCHDATA;
    typedef pVect<SEARCHDATA,SOURCE_PMAPTYPE>               SEARCHVECT;
    typedef pVect<point3d,SOURCE_PMAPTYPE>                  SEARCHPOINTS;
    
    typedef pVect<SOURCE_OUTTYPE,SOURCE_PMAPTYPE> PIECE_VECT;
    typedef SOURCE_OUTTYPE                        PIECE_DOF;
    
    typedef pMap<pMapItemSendRecv> SENDRECV;
    //@}
   
    /*! @name Internal links and flags */ //@{
  public:
    bool commLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    sVect<SEARCHVECT> dataVect;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    interpolatorMass();
    interpolatorMass(const Teuchos::RCP<communicator> & CommDev);
    interpolatorMass(communicator & CommDev);
    void setCommDev(const Teuchos::RCP<communicator> & CommDev);
    void setCommDev(communicator & CommDev);
    //@}
    
    /*! @name Computing functions */ //@{
  public:
    void findDofs(const sVect<point3d> & Y,
                  const Teuchos::RCP<SOURCEFIELD> & SourceField);
                    
    void findDofs(const sVect<point3d> & Y,
                  const SOURCEFIELD    & SourceField);
    
    void exchangeData(const sVect<point3d> & Y,
                      const Teuchos::RCP<SOURCEFIELD> & SourceField,
                      const Teuchos::RCP<TARGETFIELD> & TargetField);
                      
    void exchangeData(const sVect<point3d> & Y,
                      const SOURCEFIELD & SourceField,
                            TARGETFIELD & TargetField);
    
    void exchangeDataToll(const sVect<point3d> & Y,
                          const Teuchos::RCP<SOURCEFIELD> & SourceField,
                          const Teuchos::RCP<TARGETFIELD> & TargetField,
                          const Real & distToll = 1.0e-10);
                      
    void exchangeDataToll(const sVect<point3d> & Y,
                          const SOURCEFIELD & SourceField,
                                TARGETFIELD & TargetField,
                          const Real & distToll = 1.0e-10);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename SOURCEFIELD, typename TARGETFIELD>
interpolatorMass<SOURCEFIELD,TARGETFIELD>::
interpolatorMass() : interpolatorTrait<TARGETFIELD>::SEARCH()
{
  assert(staticAssert< TARGET_FETYPE::isNodal                >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feBaseType == scalar   >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feClass    == feStatic >::returnValue);
  
  assert(staticAssert< traitsBasic<SOURCE_OUTTYPE>::myType  == traitsBasic<TARGET_OUTTYPE>::myType >::returnValue);
  assert(staticAssert< SOURCE_PMAPTYPE::parallelType        == TARGET_PMAPTYPE::parallelType>::returnValue);
  
  commLoaded = false;
}

template<typename SOURCEFIELD, typename TARGETFIELD>
interpolatorMass<SOURCEFIELD,TARGETFIELD>::
interpolatorMass(const Teuchos::RCP<communicator> & CommDev) : interpolatorTrait<TARGETFIELD>::SEARCH(CommDev)
{
  assert(staticAssert< TARGET_FETYPE::isNodal                >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feBaseType == scalar   >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feClass    == feStatic >::returnValue);
  
  assert(staticAssert< traitsBasic<SOURCE_OUTTYPE>::myType  ==  traitsBasic<TARGET_OUTTYPE>::myType >::returnValue);
  assert(staticAssert< SOURCE_PMAPTYPE::parallelType        == TARGET_PMAPTYPE::parallelType>::returnValue);
  
  commLoaded = true;
  commDev    = CommDev;
}

template<typename SOURCEFIELD, typename TARGETFIELD>
interpolatorMass<SOURCEFIELD,TARGETFIELD>::
interpolatorMass(communicator & CommDev) : interpolatorTrait<TARGETFIELD>::SEARCH(CommDev)
{
  assert(staticAssert< TARGET_FETYPE::isNodal                >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feBaseType == scalar   >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feClass    == feStatic >::returnValue);
  
  assert(staticAssert< traitsBasic<SOURCE_OUTTYPE>::myType == traitsBasic<TARGET_OUTTYPE>::myType >::returnValue);
  assert(staticAssert< SOURCE_PMAPTYPE::parallelType       == TARGET_PMAPTYPE::parallelType>::returnValue);
  
  commLoaded = true;
  commDev    = Teuchos::rcpFromRef(CommDev);
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorMass<SOURCEFIELD,TARGETFIELD>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commLoaded = true;
  commDev    = CommDev;
  SEARCH::setCommunicator(CommDev);
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorMass<SOURCEFIELD,TARGETFIELD>::
setCommDev(communicator & CommDev)
{
  commLoaded = true;
  commDev    = Teuchos::rcpFromRef(CommDev);
  SEARCH::setCommunicator(CommDev);
}


//_________________________________________________________________________________________________
// COMPUTING FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorMass<SOURCEFIELD,TARGETFIELD>::
findDofs(const sVect<point3d> & Y,
         const Teuchos::RCP<SOURCEFIELD> & SourceField)
{
  assert(commLoaded);
  findDofs(Y,*SourceField);
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorMass<SOURCEFIELD,TARGETFIELD>::
findDofs(const sVect<point3d> & Y,
         const SOURCEFIELD    & SourceField)
{
  //Assert-----------------------------------------------------------
  assert(commLoaded);
  
  //Alloc------------------------------------------------------------
  dataVect.resize(Y.size());
  sVect<SEARCHPOINTS> nodesVect(Y.size());
  sVect<point3d> elementNodes;
  point3d P;
  
  Teuchos::RCP<SOURCE_MESH> sourceGrid = SourceField.getMesh();
  
  //Build nodes vect-------------------------------------------------
  for(UInt k=1; k <= Y.size(); ++k)
  {
    nodesVect(k).setMap(sourceGrid->getElements().getRowMap());
    nodesVect(k).setData(sVect<point3d>(sourceGrid->getNumElements()));
    nodesVect(k).updateFinder();

    for(UInt el=1; el <= sourceGrid->getNumElements(); ++el)
    {
      elementNodes = sourceGrid->getElementNodesL(el);
      P = sourceGrid->getPosition(elementNodes,Y(k));
      nodesVect(k)(el) = P;
    }
    
    SEARCH::findGlobal(nodesVect(k),dataVect(k));
    dataVect(k).bufferLids();
  }
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorMass<SOURCEFIELD,TARGETFIELD>::
exchangeData(const sVect<point3d> & Y,
             const Teuchos::RCP<SOURCEFIELD> & SourceField,
             const Teuchos::RCP<TARGETFIELD> & TargetField)
{
  exchangeData(Y,
              *SourceField,
              *TargetField);
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorMass<SOURCEFIELD,TARGETFIELD>::
exchangeData(const sVect<point3d> & Y,
             const SOURCEFIELD & SourceField,
                   TARGETFIELD & TargetField)
{
  //Assert-----------------------------------------------------------------------------------------
  assert(commLoaded);
  
  //Typedef----------------------------------------------------------------------------------------
  typedef typename TARGETFIELD::DOFCARD     TARGET_DOFCARD;
  typedef fnEL3dB<TARGETFIELD>              FUNCTIONAL;
  typedef typename FUNCTIONAL::TEST_OPTIONS TEST_OPTIONS;
  
  typedef localVectorEL3d<TARGETFIELD>      LOCALVECTOR;
  typedef vectorBuilder<LOCALVECTOR>        FEBUILDER;
  typedef typename FEBUILDER::EPETRA_VECTOR EPETRA_VECTOR;
  
  
  //Alloc------------------------------------------------------------------------------------------
  Teuchos::RCP<SOURCE_MESH> sourceGrid = SourceField.getMesh();
  Teuchos::RCP<TARGET_MESH> targetGrid = TargetField.getMesh();
  
  Teuchos::RCP<TARGET_CONNECT> targetConnect = TargetField.getMeshConnect();
  
  SOURCE_SUPPORT sourceSupport;
  TARGET_SUPPORT targetSupport;
  
  UInt el3d, numBasis, dofLid;
  Real Vs, Vt;
  point3d Yt;
  sVect<Real> basisVal;
  PIECE_DOF dofS, dofT;
  TARGET_DOFCARD targetDofCard;
  
  //Reset the target field-------------------------------------------------------------------------
  PIECE_DOF val = traitsBasic<PIECE_DOF>::getZero();
  
  for(UInt i=1; i <= TargetField.size(); ++i)
  { TargetField.setDofL(i,val); }
  
  //Build elements map-----------------------------------------------------------------------------
  pMap<pMapItemShare> elMapOwning;
  
  upgradeMap_share(sourceGrid->getElements().getRowMap(),
                   elMapOwning,
                  *commDev);
   
  //Build data pieces------------------------------------------------------------------------------
  sVect<PIECE_VECT> pieceVect(Y.size());
  sVect<PIECE_VECT> recvPieceVect(Y.size());
  sVect<SEARCHVECT> recvSearchVect(Y.size());
  sVect<point3d>    elementNodes;
  
  for(UInt k=1; k <= Y.size(); ++k)
  {
    //Piece vect startup
    pieceVect(k).setMap(sourceGrid->getElements().getRowMap());
    pieceVect(k).setData(sVect<PIECE_DOF>(sourceGrid->getElements().getRowMap().size()));
    pieceVect(k).updateFinder();
    
    //Main loop
    for(UInt el=1; el <= sourceGrid->getNumElements(); ++el)
    {
      elementNodes = sourceGrid->getElementNodesL(el);
      
      sourceSupport.setPoints(elementNodes);
      Vs = sourceSupport.template volume<STANDARD,2>();
      
      SourceField.evalL(el,Y(k),dofS);
      pieceVect(k)(el) = dofS * (Vs / Real(Y.size())) * Real(elMapOwning(el).getOwned());
    }
    
    assert(pieceVect(k).size() == dataVect(k).size());
  }
  
  //Send data pieces-------------------------------------------------------------------------------
  pMap<pMapItemSendRecv>                    sendingMap;
  pVectComm<SOURCE_OUTTYPE,SOURCE_PMAPTYPE> vectorPieceComm(commDev);
  pVectComm<SEARCHDATA,SOURCE_PMAPTYPE>     searchDataComm(commDev);
  
  for(UInt k=1; k <= Y.size(); ++k)
  {
    for(UInt el=1; el <= sourceGrid->getNumElements(); ++el)
    {
      pieceVect(k).getMapL(el).setPid(dataVect(k)(el).getElMap().getPid());
      dataVect(k).getMapL(el).setPid(dataVect(k)(el).getElMap().getPid());
    }
    
    recvPieceVect(k)  = pieceVect(k);
    recvSearchVect(k) = dataVect(k);
    
    vectorPieceComm.vectorPid(recvPieceVect(k));
    searchDataComm.vectorPid(recvSearchVect(k));
  }
  
  //Re-assemble data pieces------------------------------------------------------------------------
  for(UInt k=1; k <= Y.size(); ++k)
  {
    for(UInt i=1; i <= recvSearchVect(k).size(); ++i)
    {
      //Download searchData
      el3d = recvSearchVect(k)(i).getElMap().getLid();
      Yt   = recvSearchVect(k)(i).getLocCoord();
      
      //Area
      elementNodes = targetGrid->getElementNodesL(el3d);
      
      targetSupport.setPoints(elementNodes);
      Vt = targetSupport.template volume<STANDARD,2>();
      
      //Estimate basis data
      targetDofCard.setLocalElId(el3d);
      
      TargetField.feInterface.setCards(TargetField.dofMapper->getFeCards().getL(el3d),
                                       TargetField.feeder.getCardLocal(el3d) );
      
      numBasis = TargetField.dofMapper->isActive(targetDofCard)
               * TargetField.feInterface.getNumBasis();
               
      //Evaluate the basis
      basisVal.resize(numBasis);
      TargetField.feInterface.globalEval(Yt,basisVal);
      
      //Add the pieces
      for(UInt j=1; j <= numBasis; ++j)
      {
        targetDofCard = TargetField.feInterface.getDofCard(j);
        targetDofCard.setLocalElId(el3d);
        
        dofLid = TargetField.dofMapper->mapDofL(targetDofCard);
        val    = TargetField.getDofL(dofLid) + recvPieceVect(k)(i) * basisVal(j);

        TargetField.setDofL(dofLid,val);
      }
    }
  }
  
  //Reduce-----------------------------------------------------------------------------------------
  typedef std::plus<Real> OP;
  typedef typename TARGETFIELD::DOFVECT                    TARGETFIELD_DOFVECT;
  typedef pVectGlobalManip<TARGET_DOFTYPE,TARGET_PMAPTYPE> DOFVECT_MANIP;
  
  OP op;
  TARGETFIELD_DOFVECT dofVect = TargetField.getDofVect();
  
  DOFVECT_MANIP dofVectManip(commDev);
  dofVectManip.allReduce(dofVect,op);
  
  TargetField.setDofVect(dofVect);
  
  //Normalize the dofs-----------------------------------------------------------------------------
  
  //Operators
  TEST_OPTIONS testOptions = TargetField.getOptions();
  
  Teuchos::RCP<FUNCTIONAL> funzionale(new FUNCTIONAL(*commDev,
                                                     *targetGrid,
                                                     *targetConnect));
  funzionale->setOptions_test(testOptions);
  funzionale->startup();
  
  //Vector
  FEBUILDER feVector;
  feVector.setFunctional(funzionale);
  feVector.setCommDev(commDev);
  
  Teuchos::RCP<EPETRA_VECTOR> vector;
  feVector.buildEpetraVector(vector);
  
  epetraVector_to_pVect<SOURCE_PMAPTYPE> converter(commDev);
  pMap<SOURCE_PMAPTYPE>    listMapVol = TargetField.getListMap();
  pVect<Real,TARGET_PMAPTYPE> volVect = converter.convertAndChange(*vector,listMapVol);
  
  //Normalize
  assert(volVect.size() == TargetField.size());
  
  for(UInt i=1; i <= volVect.size(); ++i)
  {
    val = TargetField.getDofL(i) / volVect(i);
    TargetField.setDofL(i,val);
  }
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorMass<SOURCEFIELD,TARGETFIELD>::
exchangeDataToll(const sVect<point3d> & Y,
                 const Teuchos::RCP<SOURCEFIELD> & SourceField,
                 const Teuchos::RCP<TARGETFIELD> & TargetField,
                 const Real & distToll)
{
  exchangeDataToll(Y,
                  *SourceField,
                  *TargetField,
                   distToll);
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorMass<SOURCEFIELD,TARGETFIELD>::
exchangeDataToll(const sVect<point3d> & Y,
                 const SOURCEFIELD & SourceField,
                       TARGETFIELD & TargetField,
                 const Real & distToll)
{
  //Assert-----------------------------------------------------------------------------------------
  assert(commLoaded);
  
  //Typedef----------------------------------------------------------------------------------------
  typedef typename TARGETFIELD::DOFCARD     TARGET_DOFCARD;
  typedef fnEL3dB<TARGETFIELD>              FUNCTIONAL;
  typedef typename FUNCTIONAL::TEST_OPTIONS TEST_OPTIONS;
  
  typedef localVectorEL3d<TARGETFIELD>      LOCALVECTOR;
  typedef vectorBuilder<LOCALVECTOR>        FEBUILDER;
  typedef typename FEBUILDER::EPETRA_VECTOR EPETRA_VECTOR;
  
  
  //Alloc------------------------------------------------------------------------------------------
  Teuchos::RCP<SOURCE_MESH> sourceGrid = SourceField.getMesh();
  Teuchos::RCP<TARGET_MESH> targetGrid = TargetField.getMesh();
  
  Teuchos::RCP<TARGET_CONNECT> targetConnect = TargetField.getMeshConnect();
  
  SOURCE_SUPPORT sourceSupport;
  TARGET_SUPPORT targetSupport;
  
  UInt el3d, numBasis, dofLid;
  Real Vs, Vt, d;
  point3d Yt;
  sVect<Real> basisVal;
  PIECE_DOF dofS, dofT;
  TARGET_DOFCARD targetDofCard;
  
  //Reset the target field-------------------------------------------------------------------------
  PIECE_DOF val = traitsBasic<PIECE_DOF>::getZero();
  
  for(UInt i=1; i <= TargetField.size(); ++i)
  { TargetField.setDofL(i,val); }
  
  //Build elements map-----------------------------------------------------------------------------
  pMap<pMapItemShare> elMapOwning;
  
  upgradeMap_share(sourceGrid->getElements().getRowMap(),
                   elMapOwning,
                  *commDev);
   
  //Build data pieces------------------------------------------------------------------------------
  sVect<PIECE_VECT> pieceVect(Y.size());
  sVect<PIECE_VECT> recvPieceVect(Y.size());
  sVect<SEARCHVECT> recvSearchVect(Y.size());
  sVect<point3d>    elementNodes;
  
  for(UInt k=1; k <= Y.size(); ++k)
  {
    //Piece vect startup
    pieceVect(k).setMap(sourceGrid->getElements().getRowMap());
    pieceVect(k).setData(sVect<PIECE_DOF>(sourceGrid->getElements().getRowMap().size()));
    pieceVect(k).updateFinder();
    
    //Main loop
    for(UInt el=1; el <= sourceGrid->getNumElements(); ++el)
    {
      elementNodes = sourceGrid->getElementNodesL(el);
      
      sourceSupport.setPoints(elementNodes);
      Vs = sourceSupport.template volume<STANDARD,2>();
      
      SourceField.evalL(el,Y(k),dofS);
      
      d = dataVect(k).getDataL(el).getDistance();
      pieceVect(k)(el) = dofS * (Vs / Real(Y.size())) * Real(elMapOwning(el).getOwned()) * Real(d <= distToll);
    }
    
    assert(pieceVect(k).size() == dataVect(k).size());
  }
  
  //Send data pieces-------------------------------------------------------------------------------
  pMap<pMapItemSendRecv>                    sendingMap;
  pVectComm<SOURCE_OUTTYPE,SOURCE_PMAPTYPE> vectorPieceComm(commDev);
  pVectComm<SEARCHDATA,SOURCE_PMAPTYPE>     searchDataComm(commDev);
  
  for(UInt k=1; k <= Y.size(); ++k)
  {
    for(UInt el=1; el <= sourceGrid->getNumElements(); ++el)
    {
      pieceVect(k).getMapL(el).setPid(dataVect(k)(el).getElMap().getPid());
      dataVect(k).getMapL(el).setPid(dataVect(k)(el).getElMap().getPid());
    }
    
    recvPieceVect(k)  = pieceVect(k);
    recvSearchVect(k) = dataVect(k);
    
    vectorPieceComm.vectorPid(recvPieceVect(k));
    searchDataComm.vectorPid(recvSearchVect(k));
  }
  
  //Re-assemble data pieces------------------------------------------------------------------------
  for(UInt k=1; k <= Y.size(); ++k)
  {
    for(UInt i=1; i <= recvSearchVect(k).size(); ++i)
    {
      //Download searchData
      el3d = recvSearchVect(k)(i).getElMap().getLid();
      Yt   = recvSearchVect(k)(i).getLocCoord();
      
      //Area
      elementNodes = targetGrid->getElementNodesL(el3d);
      
      targetSupport.setPoints(elementNodes);
      Vt = targetSupport.template volume<STANDARD,2>();
      
      //Estimate basis data
      targetDofCard.setLocalElId(el3d);
      
      TargetField.feInterface.setCards(TargetField.dofMapper->getFeCards().getL(el3d),
                                       TargetField.feeder.getCardLocal(el3d) );
      
      numBasis = TargetField.dofMapper->isActive(targetDofCard)
               * TargetField.feInterface.getNumBasis();
               
      //Evaluate the basis
      basisVal.resize(numBasis);
      TargetField.feInterface.globalEval(Yt,basisVal);
      
      //Add the pieces
      for(UInt j=1; j <= numBasis; ++j)
      {
        targetDofCard = TargetField.feInterface.getDofCard(j);
        targetDofCard.setLocalElId(el3d);
        
        dofLid = TargetField.dofMapper->mapDofL(targetDofCard);
        val    = TargetField.getDofL(dofLid) + recvPieceVect(k)(i) * basisVal(j);

        TargetField.setDofL(dofLid,val);
      }
    }
  }
  
  //Reduce-----------------------------------------------------------------------------------------
  typedef std::plus<Real> OP;
  typedef typename TARGETFIELD::DOFVECT                    TARGETFIELD_DOFVECT;
  typedef pVectGlobalManip<TARGET_DOFTYPE,TARGET_PMAPTYPE> DOFVECT_MANIP;
  
  OP op;
  TARGETFIELD_DOFVECT dofVect = TargetField.getDofVect();
  
  DOFVECT_MANIP dofVectManip(commDev);
  dofVectManip.allReduce(dofVect,op);
  
  TargetField.setDofVect(dofVect);
  
  //Normalize the dofs-----------------------------------------------------------------------------
  
  //Operators
  TEST_OPTIONS testOptions = TargetField.getOptions();
  
  Teuchos::RCP<FUNCTIONAL> funzionale(new FUNCTIONAL(*commDev,
                                                     *targetGrid,
                                                     *targetConnect));
  funzionale->setOptions_test(testOptions);
  funzionale->startup();
  
  //Vector
  FEBUILDER feVector;
  feVector.setFunctional(funzionale);
  feVector.setCommDev(commDev);
  
  Teuchos::RCP<EPETRA_VECTOR> vector;
  feVector.buildEpetraVector(vector);
  
  epetraVector_to_pVect<SOURCE_PMAPTYPE> converter(commDev);
  pMap<SOURCE_PMAPTYPE>    listMapVol = TargetField.getListMap();
  pVect<Real,TARGET_PMAPTYPE> volVect = converter.convertAndChange(*vector,listMapVol);
  
  //Normalize
  assert(volVect.size() == TargetField.size());
  
  for(UInt i=1; i <= volVect.size(); ++i)
  {
    val = TargetField.getDofL(i) / volVect(i);
    TargetField.setDofL(i,val);
  }
}

#endif
