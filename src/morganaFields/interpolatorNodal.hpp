/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTERPOLATORNODAL_HPP
#define INTERPOLATORNODAL_HPP

#include "morganaFields.hpp"
#include "traitsInterpolator.hpp"



/*! Interpolator for the nodal static fields. The \c sourceField should be \c static and \c nodal since the finite element should
 implement the reference nodes. The \c targetField should only implement the \c eval method, therefore all the fields work.*/
template<typename SOURCEFIELD, typename TARGETFIELD>
class interpolatorNodal : public interpolatorTrait<SOURCEFIELD>::SEARCH
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename SOURCEFIELD::GEOSHAPE     SOURCE_GEOSHAPE;
    typedef typename SOURCEFIELD::OUTTYPE      SOURCE_OUTTYPE;
    typedef typename SOURCEFIELD::FIELD_FETYPE SOURCE_FETYPE;
    typedef typename SOURCEFIELD::PMAPTYPE     SOURCE_PMAPTYPE;
    typedef typename SOURCEFIELD::FEINTERFACE  SOURCE_FEINTERFACE;
    typedef typename interpolatorTrait<SOURCEFIELD>::MESH    SOURCE_MESH;
    typedef typename interpolatorTrait<SOURCEFIELD>::CONNECT SOURCE_CONNECT;
    
    typedef typename TARGETFIELD::GEOSHAPE     TARGET_GEOSHAPE;
    typedef typename TARGETFIELD::OUTTYPE      TARGET_OUTTYPE;
    typedef typename TARGETFIELD::FIELD_FETYPE TARGET_FETYPE;
    typedef typename TARGETFIELD::PMAPTYPE     TARGET_PMAPTYPE;
    typedef typename TARGETFIELD::FEINTERFACE  TARGET_FEINTERFACE;
    typedef typename interpolatorTrait<TARGETFIELD>::MESH    TARGET_MESH;
    typedef typename interpolatorTrait<TARGETFIELD>::CONNECT TARGET_CONNECT;
    
    typedef typename interpolatorTrait<SOURCEFIELD>::SEARCH SEARCH;
    typedef searchData<SOURCE_PMAPTYPE>                     SEARCHDATA;
    //@}
   
    /*! @name Internal links and flags */ //@{
  public:
    bool commLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    pVect<SEARCHDATA,SOURCE_PMAPTYPE> dataVect;
    //@}  
    
    /*! @name Constructors and set functions */ //@{
  public:
    interpolatorNodal();
    interpolatorNodal(const Teuchos::RCP<communicator> & CommDev);
    interpolatorNodal(communicator & CommDev);
    void setCommDev(const Teuchos::RCP<communicator> & CommDev);
    void setCommDev(communicator & CommDev);
    //@}
    
    /*! @name Computing functions */ //@{
  public:
    void findDofs(const Teuchos::RCP<TARGETFIELD> & TargetField);
    void findDofs(const TARGETFIELD & TargetField);
    void exchangeData(const Teuchos::RCP<SOURCEFIELD> & SourceField, const Teuchos::RCP<TARGETFIELD> & TargetField);
    void exchangeData(const SOURCEFIELD & SourceField, const TARGETFIELD & TargetField);
    //@}
};


template<typename SOURCEFIELD, typename TARGETFIELD>
interpolatorNodal<SOURCEFIELD,TARGETFIELD>::
interpolatorNodal() : interpolatorTrait<SOURCEFIELD>::SEARCH()
{
  assert(staticAssert< TARGET_FETYPE::isNodal                >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feBaseType == scalar   >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feClass    == feStatic >::returnValue);
  
  assert(staticAssert< traitsBasic<SOURCE_OUTTYPE>::myType  ==  traitsBasic<TARGET_OUTTYPE>::myType >::returnValue);
  assert(staticAssert< SOURCE_PMAPTYPE::parallelType        == TARGET_PMAPTYPE::parallelType>::returnValue);
  
  commLoaded = false;
}

template<typename SOURCEFIELD, typename TARGETFIELD>
interpolatorNodal<SOURCEFIELD,TARGETFIELD>::
interpolatorNodal(const Teuchos::RCP<communicator> & CommDev) : interpolatorTrait<SOURCEFIELD>::SEARCH(CommDev)
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
interpolatorNodal<SOURCEFIELD,TARGETFIELD>::
interpolatorNodal(communicator & CommDev) : interpolatorTrait<SOURCEFIELD>::SEARCH(CommDev)
{
  assert(staticAssert< TARGET_FETYPE::isNodal                >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feBaseType == scalar   >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feClass    == feStatic >::returnValue);
  
  assert(staticAssert< traitsBasic<SOURCE_OUTTYPE>::myType  ==  traitsBasic<TARGET_OUTTYPE>::myType >::returnValue);
  assert(staticAssert< SOURCE_PMAPTYPE::parallelType        == TARGET_PMAPTYPE::parallelType>::returnValue);
  
  commLoaded = true;
  commDev    = Teuchos::rcpFromRef(CommDev);
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorNodal<SOURCEFIELD,TARGETFIELD>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commLoaded = true;
  commDev    = CommDev;
  SEARCH::setCommDev(CommDev);
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorNodal<SOURCEFIELD,TARGETFIELD>::
setCommDev(communicator & CommDev)
{
  commLoaded = true;
  commDev    = Teuchos::rcpFromRef(CommDev);
  SEARCH::setCommDev(CommDev);
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorNodal<SOURCEFIELD,TARGETFIELD>::
findDofs(const Teuchos::RCP<TARGETFIELD> & TargetField)
{
  assert(commLoaded);
  findDofs(*TargetField);
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorNodal<SOURCEFIELD,TARGETFIELD>::
findDofs(const TARGETFIELD & TargetField)
{
  //Typedefs---------------------------------------------------------
  typedef SOURCE_PMAPTYPE                     PMAPTYPE;
  typedef geoMapInterface<TARGET_GEOSHAPE>    GEOMAPINTERFACE;
  typedef pMap<pMapItemSendRecv>              SENDRECV;
  typedef typename TARGETFIELD::DOFCARD       TARGET_DOFCARD;
  typedef typename TARGETFIELD::FIELD_DOFTYPE TARGET_DOFTYPE;
  
  //Assert-----------------------------------------------------------
  assert(commLoaded);
  assert(TargetField.getStartupOk());
  
  //Alloc structures-------------------------------------------------
  Teuchos::RCP<TARGET_MESH> targetGrid = TargetField.getMesh();
  
  //Queries generation-----------------------------------------------
  TARGET_DOFCARD          targetDofCard;
  GEOMAPINTERFACE         geoInterface;
  pVect<point3d,PMAPTYPE> nodesVect;
  
  sVect<point3d> points;
  point3d        Y, P;
  UInt lid, gid, sid, rid;
  
  nodesVect.setMap(TargetField.getDofVect().getMapRef());
  nodesVect.setData(sVect<point3d>(TargetField.getDofVect().size()));
  nodesVect.updateFinder();
  
  //Loop on the elements
  for(UInt i=1; i <= targetGrid->getNumElements(); ++i)
  {
    points = targetGrid->getElementNodesL(i);
    
   //Loop on the local basis
   for(UInt j=1; j <= TARGET_FETYPE::numBasis; ++j)
   {
      Y = TARGET_FETYPE::getRefNode(j);
      P = geoInterface.getPosition(points,Y);
      
      targetDofCard = TARGET_FETYPE::getDofCard(j);
      targetDofCard.setLocalElId(i);
      
      if(TargetField.getDofMapper().isActive(targetDofCard))
      {
        lid = TargetField.getDofMapper().mapDofL(targetDofCard);

        assert(lid <= nodesVect.size());
        nodesVect(lid) = P;
      }
    }
  }

  //Searching--------------------------------------------------------
  pVect<SEARCHDATA,PMAPTYPE> tempDataVect;
  
  SEARCH::findGlobal(nodesVect,tempDataVect);
  tempDataVect.bufferLids();

  //Sending to the correct process-----------------------------------
  SENDRECV sendMap;
  dataVect.clear();
  
  for(UInt i=1; i <= nodesVect.size(); ++i)
  {
    rid = tempDataVect(i).getElMap().getPid();
    sid = commDev->rank();
    lid = tempDataVect.getMapL(i).getLid();
    gid = tempDataVect.getMapL(i).getGid();
    
    if(int(rid) != commDev->rank())
    { sendMap.push_back(pMapItemSendRecv(lid,gid,sid,rid)); }
    else
    { dataVect.push_back(tempDataVect.getMapL(i),tempDataVect.getDataL(i)); }
  }
  
  pVect<SEARCHDATA,PMAPTYPE> inVect;                    // Transmit the data to the correct pid
  
  pVectComm<SEARCHDATA,PMAPTYPE> dataComm(commDev);
  dataComm.sendRecv(sendMap,tempDataVect,inVect);
  inVect.restoreLids();
  
  for(UInt i=1 ; i <= inVect.size(); ++i)               // Merging with the original data
  { dataVect.push_back(inVect.getMapL(i), inVect.getDataL(i)); }
  
  dataVect.bufferLids();
  
  for(UInt i=1; i <= dataVect.size(); ++i)
  { assert(dataVect(i).getElMap().getLid() != 0); }
}


template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorNodal<SOURCEFIELD,TARGETFIELD>::
exchangeData(const Teuchos::RCP<SOURCEFIELD> & SourceField, const Teuchos::RCP<TARGETFIELD> & TargetField)
{
  assert(commLoaded);
  exchangeData(*SourceField,*TargetField);
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorNodal<SOURCEFIELD,TARGETFIELD>::
exchangeData(const SOURCEFIELD & SourceField, const TARGETFIELD & TargetField)
{
  //Typedefs---------------------------------------------------------
  typedef SOURCE_PMAPTYPE                     PMAPTYPE;
  typedef geoMapInterface<TARGET_GEOSHAPE>    GEOMAPINTERFACE;
  typedef pMap<pMapItemSendRecv>              SENDRECV;
  typedef typename TARGETFIELD::DOFCARD       TARGET_DOFCARD;
  typedef typename TARGETFIELD::FIELD_DOFTYPE TARGET_DOFTYPE;
  
  //Evaluation of the field------------------------------------------
  pVect<TARGET_DOFTYPE,PMAPTYPE> newDofVector(dataVect.size());
  newDofVector.setMap(dataVect.getMapRef());
  newDofVector.updateFinder();
  
  UInt elL;
  point3d Y;
  TARGET_DOFTYPE dofVal;
  
  for(UInt i=1; i <= dataVect.size(); ++i)
  {
    assert(commDev->rank() == dataVect(i).getElMap().getPid());
    
    elL = dataVect(i).getElMap().getLid();
    Y   = dataVect(i).getLocCoord();
    
    assert(elL != 0); 
    SourceField.evalL(elL,Y,dofVal);
    
    newDofVector(i) = dofVal;
  }
  
  //Transmit back to correct pid-------------------------------------
  pVectComm<TARGET_DOFTYPE,PMAPTYPE> dofComm(commDev);
  dofComm.vectorPid(newDofVector);
  
  newDofVector.restoreLids();  
  
  pVectManip<TARGET_DOFTYPE,PMAPTYPE> dofManipulator;
  dofManipulator.setIndexing(newDofVector);
  
  //Loading in the field---------------------------------------------
  TargetField.setDofVect(newDofVector);
}

#endif
