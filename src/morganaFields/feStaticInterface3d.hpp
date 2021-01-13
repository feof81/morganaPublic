/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FESTATICINTERFACE3D_HPP
#define FESTATICINTERFACE3D_HPP

#include "traitsMultiply.h"
#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"


/*! Interface class for the three dimensional finite static element fields */
template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER = dms3d_vectMajor, dms3d_mode MODE = dms3d_allMode>
class feStaticInterface3d : public FETYPE, public dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename  FETYPE::FETYPE_PMAPTYPE PMAPTYPE;
    typedef typename  FETYPE::GEOSHAPE        GEOSHAPE;
    typedef typename  FETYPE::FECARD          FECARD;
    typedef typename  FETYPE::ELCARD          ELCARD;
    typedef typename  FETYPE::BASETYPE        BASETYPE;
    typedef typename  FETYPE::DOFCARD         DOFCARD;
    
    typedef mesh3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH3D;
    typedef connect3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT3D;
    
    typedef typename traitsMultiply<BASETYPE,DOFTYPE>::DIRECTTYPE  OUTTYPE;
    typedef dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>              DOFMAPPER;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    elCardFeeder3d<GEOSHAPE,PMAPTYPE> feeder;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    feStaticInterface3d();
    feStaticInterface3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    feStaticInterface3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
    feStaticInterface3d(const feStaticInterface3d & Interface);
    feStaticInterface3d operator=(const feStaticInterface3d & Interface);
    void clone(const DOFMAPPER & DofMapper);
    void setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    void setGeometry(MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
    //@}
    
    /*! @name Info */ //@{
  public:
    bool feTest() const;
    UInt getNumBasis() const;
    sVect<DOFCARD> getDofCards() const;
    //@}
    
    /*! @name Evaluation functions */ //@{
  public:
    void eval(const UInt & elL, const UInt & I, const UInt & J, const FECARD & FECard, const ELCARD & ELCard, const point3d & Y, sVect<OUTTYPE> & val);
    void grad(const UInt & elL, const UInt & I, const UInt & J, const FECARD & FECard, const ELCARD & ELCard, const point3d & Y, sVect<OUTTYPE> & gradX, sVect<OUTTYPE> & gradY, sVect<OUTTYPE> & gradZ);
    //@}
    
    /*! @name Fe interface */ //@{
  public:
    const FECARD & getFeCardL(const UInt & elL);
    const FECARD & getFeCardG(const UInt & elG);
    UInt  feNumBasisL(const UInt & elL);
    UInt  feNumBasisG(const UInt & elG);
    void  feIndicesLL(const UInt & elL, const UInt & I, const UInt & J, sVect<UInt> & indices);
    void  feIndicesGL(const UInt & elG, const UInt & I, const UInt & J, sVect<UInt> & indices);
    void  feIndicesLG(const UInt & elL, const UInt & I, const UInt & J, sVect<UInt> & indices);
    void  feIndicesGG(const UInt & elG, const UInt & I, const UInt & J, sVect<UInt> & indices);
    void  feEvalL(const UInt & elL, const UInt & I, const UInt & J, const point3d & Y, sVect<OUTTYPE> & val);
    void  feEvalG(const UInt & elG, const UInt & I, const UInt & J, const point3d & Y, sVect<OUTTYPE> & val);
    void  feGradL(const UInt & elL, const UInt & I, const UInt & J, const point3d & Y, sVect<OUTTYPE> & gradX, sVect<OUTTYPE> & gradY, sVect<OUTTYPE> & gradZ);
    void  feGradG(const UInt & elG, const UInt & I, const UInt & J, const point3d & Y, sVect<OUTTYPE> & gradX, sVect<OUTTYPE> & gradY, sVect<OUTTYPE> & gradZ);
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticInterface3d() :
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>()
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticInterface3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d) :
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  feeder.setGeometry(Grid3d,ConnedGrid3d);
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticInterface3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d) :
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  feeder.setGeometry(Grid3d,ConnedGrid3d);
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticInterface3d(const feStaticInterface3d & Interface) : dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>(Interface)
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  feeder.setGeometry(Interface.grid3d, Interface.connedGrid3d);
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
operator=(const feStaticInterface3d & Interface)
{
  feeder.setGeometry(Interface.grid3d, Interface.connectGrid3d);
  DOFMAPPER::operator=(Interface);
  
  return(*this);
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
void
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
clone(const DOFMAPPER & DofMapper)
{
  feeder.setGeometry(DofMapper.grid3d, DofMapper.connectGrid3d);
  DOFMAPPER::operator=(DofMapper);
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
void 
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d)
{
  feeder.setGeometry(Grid3d,ConnectGrid3d);
  DOFMAPPER::setGeometry(Grid3d,ConnectGrid3d);
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
void
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
setGeometry(MESH3D & Grid3d, CONNECT3D & ConnectGrid3d)
{
  feeder.setGeometry(Grid3d,ConnectGrid3d);
  DOFMAPPER::setGeometry(Grid3d,ConnectGrid3d);
}



//_________________________________________________________________________________________________
// INFO
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
bool
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feTest() const
{  
  UInt dofVertex = 0;
  UInt dofEdge   = 0;
  UInt dofFace   = 0;
  UInt dofVolume = 0;
  
  for(UInt i=1; i <= FETYPE::numBasis; ++i)
  {
    switch(FETYPE::getDofCard(i).getGeoType())
    {
      case VERTEX :
	++dofVertex;
	break;
	
      case EDGE :
	++dofEdge;
	break;
	
      case FACE :
	++dofFace;
        break;
	
      case VOLUME :
	++dofVolume;
	break;
    }
  }
  
  return(
    ( (geoMapInterface<GEOSHAPE>::getNumVertices() ) * FETYPE::dofPerVertex == dofVertex) &&
    ( (geoMapInterface<GEOSHAPE>::getNumEdges()    ) * FETYPE::dofPerEdge   == dofEdge)   &&
    ( (geoMapInterface<GEOSHAPE>::getNumFaces()    ) * FETYPE::dofPerFace   == dofFace)   &&
                                                       FETYPE::dofPerVolume == dofVolume  );
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
UInt
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>:: 
getNumBasis() const
{
  return(FETYPE::numBasis);
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
sVect<typename feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::DOFCARD>
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>:: 
getDofCards() const
{
  sVect<DOFCARD> outVect(FETYPE::numBasis);
  
  for(UInt i=1; i <= FETYPE::numBasis; ++i)
  {
    outVect(i) = FETYPE::getDofCard(i);
  }
  
  return(outVect);
}



//_________________________________________________________________________________________________
// EVALUATIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
void
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
eval(const UInt & elL, const UInt & I, const UInt & J, const FECARD & FECard, const ELCARD & ELCard, const point3d & Y, sVect<OUTTYPE> & val)
{
  sVect<BASETYPE> basis(FETYPE::numBasis);
  FETYPE::setCards(FECard,ELCard);
  FETYPE::globalEval(Y,basis);
  
  DOFTYPE dof = traitsBasic<DOFTYPE>::getUnityIJ(I,J);
  
  for(UInt i=1; i <= feNumBasisL(elL); ++i)
  {
    val(i) = traitsMultiply<BASETYPE,DOFTYPE>::multiply(basis(i),dof);
  }
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
void
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>:: 
grad(const UInt & elL, const UInt & I, const UInt & J, const FECARD & FECard, const ELCARD & ELCard, const point3d & Y, sVect<OUTTYPE> & gradX, sVect<OUTTYPE> & gradY, sVect<OUTTYPE> & gradZ)
{
  sVect<BASETYPE> basisX(FETYPE::numBasis);
  sVect<BASETYPE> basisY(FETYPE::numBasis);
  sVect<BASETYPE> basisZ(FETYPE::numBasis);
  
  FETYPE::setCards(FECard,ELCard);
  FETYPE::globalGrad(Y,basisX,basisY,basisZ);

  DOFTYPE dof = traitsBasic<DOFTYPE>::getUnityIJ(I,J);
  
  for(UInt i=1; i <= feNumBasisL(elL); ++i)
  {
    gradX(i) = traitsMultiply<BASETYPE,DOFTYPE>::multiply(basisX(i),dof);
    gradY(i) = traitsMultiply<BASETYPE,DOFTYPE>::multiply(basisY(i),dof);
    gradZ(i) = traitsMultiply<BASETYPE,DOFTYPE>::multiply(basisZ(i),dof);
  }
}



//_________________________________________________________________________________________________
// FE INTERFACE
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
const typename feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::FECARD &
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
getFeCardL(const UInt & elL)
{
  typedef dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(elL >= 1);
  assert(elL <= DOFMAPPER::feCards.size());
  
  return(DOFMAPPER::feCards.getDataL(elL));
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
const typename feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::FECARD &
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
getFeCardG(const UInt & elG)
{
  typedef dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(DOFMAPPER::feCards.isG(elG));
  
  return(DOFMAPPER::feCards.getDataG(elG));
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
UInt
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feNumBasisL(const UInt & elL)
{
  typedef dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  DOFCARD dofCard;
  dofCard.setLocalElId(elL);
  
  return(FETYPE::numBasis * DOFMAPPER::isActive(dofCard));
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
UInt
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feNumBasisG(const UInt & elG)
{
  typedef dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  assert(DOFMAPPER::feCards.isG(elG));
  
  UInt elL = DOFMAPPER::grid3d->getElements().getRowMapG(elG).getLid();
  DOFCARD dofCard;
  dofCard.setLocalElId(elL);
  
  return(FETYPE::numBasis * DOFMAPPER::isActive(dofCard));
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
void
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feIndicesLL(const UInt & elL, const UInt & I, const UInt & J, sVect<UInt> & indices)
{
  typedef dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(indices.size() == feNumBasisL(elL));
  assert(elL >= 1);
  assert(elL <= DOFMAPPER::feCards.size());
  
  //Alloc
  DOFCARD dofCard;
  
  //Filling the index vector
  for(UInt i=1; i <= feNumBasisL(elL); ++i)
  {
    dofCard = FETYPE::getDofCard(i);
    dofCard.setLocalElId(elL);
    indices(i) = DOFMAPPER::mapListL(I,J,dofCard);
  }
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
void
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feIndicesGL(const UInt & elG, const UInt & I, const UInt & J, sVect<UInt> & indices)
{
  typedef dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(indices.size() == feNumBasisG(elG));
  assert(DOFMAPPER::feCards.isG(elG));
  
  //Alloc
  DOFCARD dofCard;
  UInt elL = DOFMAPPER::grid3d->getElements().getRowMapG(elG).getLid();
  
  //Filling the index vector
  for(UInt i=1; i <= feNumBasisL(elL); ++i)
  {
    dofCard = FETYPE::getDofCard(i);
    dofCard.setLocalElId(elL);
    indices(i) = DOFMAPPER::mapListL(I,J,dofCard);
  }
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
void
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feIndicesLG(const UInt & elL, const UInt & I, const UInt & J, sVect<UInt> & indices)
{
  typedef dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(indices.size() == feNumBasisL(elL));
  assert(elL >= 1);
  assert(elL <= DOFMAPPER::feCards.size());
  
  //Alloc
  DOFCARD dofCard;
  
  //Filling the index vector
  for(UInt i=1; i <= feNumBasisL(elL); ++i)
  {
    dofCard = FETYPE::getDofCard(i);
    dofCard.setLocalElId(elL);
    indices(i) = DOFMAPPER::mapListG(I,J,dofCard);
  }
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
void
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feIndicesGG(const UInt & elG, const UInt & I, const UInt & J, sVect<UInt> & indices)
{
  typedef dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(indices.size() == feNumBasisG(elG));
  assert(DOFMAPPER::feCards.isG(elG));
  
  //Alloc
  DOFCARD dofCard;
  UInt elL = DOFMAPPER::grid3d->getElements().getRowMapG(elG).getLid();
  
  //Filling the index vector
  for(UInt i=1; i <= feNumBasisL(elL); ++i)
  {
    dofCard = FETYPE::getDofCard(i);
    dofCard.setLocalElId(elL);
    indices(i) = DOFMAPPER::mapListG(I,J,dofCard);
  }
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
void
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feEvalL(const UInt & elL, const UInt & I, const UInt & J, const point3d & Y, sVect<OUTTYPE> & val)
{
  typedef dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(elL >= 1);
  assert(elL <= DOFMAPPER::feCards.size());
  assert(val.size() == feNumBasisL(elL));
  
  //Evaluation
  eval(elL, I,J, DOFMAPPER::feCards.getDataL(elL), feeder.getCardLocal(elL), Y, val);
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
void
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feEvalG(const UInt & elG, const UInt & I, const UInt & J, const point3d & Y, sVect<OUTTYPE> & val)
{
  typedef dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(val.size() == feNumBasisG(elG));
  assert(DOFMAPPER::feCards.isG(elG));
  
  //Alloc
  UInt elL = DOFMAPPER::grid3d->getElements().getRowMapG(elG).getLid();
  
  //Evaluation
  eval(elL,I,J, DOFMAPPER::feCards.getDataL(elL), feeder.getCardLocal(elL), Y, val);
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
void
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feGradL(const UInt & elL, const UInt & I, const UInt & J, const point3d & Y, sVect<OUTTYPE> & gradX, sVect<OUTTYPE> & gradY, sVect<OUTTYPE> & gradZ)
{
  typedef dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(elL >= 1);
  assert(elL <= DOFMAPPER::feCards.size());
  assert(gradX.size() == feNumBasisL(elL));
  assert(gradY.size() == feNumBasisL(elL));
  assert(gradZ.size() == feNumBasisL(elL));
  
  //Evaluation
  sVect<BASETYPE> basisX(FETYPE::numBasis);
  sVect<BASETYPE> basisY(FETYPE::numBasis);
  sVect<BASETYPE> basisZ(FETYPE::numBasis);
  
  FETYPE::setCards(DOFMAPPER::feCards.getDataL(elL), feeder.getCardLocal(elL));
  FETYPE::globalGrad(Y,basisX,basisY,basisZ);

  DOFTYPE dof = traitsBasic<DOFTYPE>::getUnityIJ(I,J);
  
  for(UInt i=1; i <= feNumBasisL(elL); ++i)
  {
    gradX(i) = traitsMultiply<BASETYPE,DOFTYPE>::multiply(basisX(i),dof);
    gradY(i) = traitsMultiply<BASETYPE,DOFTYPE>::multiply(basisY(i),dof);
    gradZ(i) = traitsMultiply<BASETYPE,DOFTYPE>::multiply(basisZ(i),dof);
  }
}

template<typename FETYPE, typename DOFTYPE,  dms3d_order ORDER, dms3d_mode MODE>
void
feStaticInterface3d<FETYPE,DOFTYPE,ORDER,MODE>::
feGradG(const UInt & elG, const UInt & I, const UInt & J, const point3d & Y, sVect<OUTTYPE> & gradX, sVect<OUTTYPE> & gradY, sVect<OUTTYPE> & gradZ)
{
  typedef dofMapStatic3d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(DOFMAPPER::feCards.isG(elG));
  assert(gradX.size() == feNumBasisG(elG));
  assert(gradY.size() == feNumBasisG(elG));
  assert(gradZ.size() == feNumBasisG(elG));
  
  //Mapping
  UInt elL = DOFMAPPER::grid3d->getElements().getRowMapG(elG).getLid();
  
  //Evaluation  
  sVect<BASETYPE> basisX(FETYPE::numBasis);
  sVect<BASETYPE> basisY(FETYPE::numBasis);
  sVect<BASETYPE> basisZ(FETYPE::numBasis);
  
  FETYPE::setCards(DOFMAPPER::feCards.getDataL(elL), feeder.getCardLocal(elL));
  FETYPE::globalGrad(Y,basisX,basisY,basisZ);

  DOFTYPE dof = traitsBasic<DOFTYPE>::getUnityIJ(I,J);
  
  for(UInt i=1; i <= feNumBasisL(elL); ++i)
  {
    gradX(i) = traitsMultiply<BASETYPE,DOFTYPE>::multiply(basisX(i),dof);
    gradY(i) = traitsMultiply<BASETYPE,DOFTYPE>::multiply(basisY(i),dof);
    gradZ(i) = traitsMultiply<BASETYPE,DOFTYPE>::multiply(basisZ(i),dof);
  }
}

#endif
