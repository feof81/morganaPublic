/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FESTATICINTERFACE2D_HPP
#define FESTATICINTERFACE2D_HPP

#include "traitsMultiply.h"
#include "dofMapStatic2d.hpp"
#include "elCardFeeder2d.hpp"


/*! Interface class for the three dimensional finite static element fields */
template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER = dms2d_vectMajor, dms2d_mode MODE = dms2d_allMode>
class feStaticInterface2d : public FETYPE, public dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename  FETYPE::FETYPE_PMAPTYPE PMAPTYPE;
    typedef typename  FETYPE::GEOSHAPE        GEOSHAPE;
    typedef typename  FETYPE::FECARD          FECARD;
    typedef typename  FETYPE::ELCARD          ELCARD;
    typedef typename  FETYPE::BASETYPE        BASETYPE;
    typedef typename  FETYPE::DOFCARD         DOFCARD;
    
    typedef mesh2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH2D;
    typedef connect2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT2D;
    
    typedef typename traitsMultiply<BASETYPE,DOFTYPE>::DIRECTTYPE  OUTTYPE;
    typedef dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>              DOFMAPPER;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    elCardFeeder2d<GEOSHAPE,PMAPTYPE> feeder;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    feStaticInterface2d();
    feStaticInterface2d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d);
    feStaticInterface2d(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d);
    feStaticInterface2d(const feStaticInterface2d & Interface);
    feStaticInterface2d operator=(const feStaticInterface2d & Interface);
    void clone(const DOFMAPPER & DofMapper);
    void setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d);
    void setGeometry(MESH2D & Grid2d, CONNECT2D & ConnedGrid2d);    
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
template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticInterface2d() :
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>()
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticInterface2d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d) :
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>(CommDev,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticInterface2d(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d) :
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>(CommDev,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticInterface2d(const feStaticInterface2d & Interface) : dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>(Interface)
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  feeder.setGeometry(Interface.grid2d, Interface.connedGrid2d);
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
operator=(const feStaticInterface2d & Interface)
{
  feeder.setGeometry(Interface.grid2d, Interface.connectGrid2d);
  DOFMAPPER::operator=(Interface);
  
  return(*this);
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
void
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
clone(const DOFMAPPER & DofMapper)
{
  feeder.setGeometry(DofMapper.grid2d, DofMapper.connectGrid2d);
  DOFMAPPER::operator=(DofMapper);
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
void 
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d)
{
  feeder.setGeometry(Grid2d,ConnedGrid2d);
  DOFMAPPER::setGeometry(Grid2d,ConnedGrid2d);
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
void
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
setGeometry(MESH2D & Grid2d, CONNECT2D & ConnedGrid2d)
{
  feeder.setGeometry(Grid2d,ConnedGrid2d);
  DOFMAPPER::setGeometry(Grid2d,ConnedGrid2d);
}


//_________________________________________________________________________________________________
// INFO
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
bool
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
UInt
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>:: 
getNumBasis() const
{
  return(FETYPE::numBasis);
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
sVect<typename feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::DOFCARD>
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>:: 
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
template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
void
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
void
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>:: 
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
template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
const typename feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::FECARD &
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
getFeCardL(const UInt & elL)
{
  typedef dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(elL >= 1);
  assert(elL <= DOFMAPPER::feCards.size());
  
  return(DOFMAPPER::feCards.getDataL(elL));
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
const typename feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::FECARD &
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
getFeCardG(const UInt & elG)
{
  typedef dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(DOFMAPPER::feCards.isG(elG));
  
  return(DOFMAPPER::feCards.getDataG(elG));
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
UInt
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
feNumBasisL(const UInt & elL)
{
  typedef dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  DOFCARD dofCard;
  dofCard.setLocalElId(elL);
  
  return(FETYPE::numBasis * DOFMAPPER::isActive(dofCard));
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
UInt
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
feNumBasisG(const UInt & elG)
{
  typedef dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  assert(DOFMAPPER::feCards.isG(elG));
  
  UInt elL = DOFMAPPER::grid2d->getElements().getRowMapG(elG).getLid();
  DOFCARD dofCard;
  dofCard.setLocalElId(elL);
  
  return(FETYPE::numBasis * DOFMAPPER::isActive(dofCard));
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
void
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
feIndicesLL(const UInt & elL, const UInt & I, const UInt & J, sVect<UInt> & indices)
{
  typedef dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
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

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
void
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
feIndicesGL(const UInt & elG, const UInt & I, const UInt & J, sVect<UInt> & indices)
{
  typedef dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(indices.size() == feNumBasisG(elG));
  assert(DOFMAPPER::feCards.isG(elG));
  
  //Alloc
  DOFCARD dofCard;
  UInt elL = DOFMAPPER::grid2d->getElements().getRowMapG(elG).getLid();
  
  //Filling the index vector
  for(UInt i=1; i <= feNumBasisL(elL); ++i)
  {
    dofCard = FETYPE::getDofCard(i);
    dofCard.setLocalElId(elL);
    indices(i) = DOFMAPPER::mapListL(I,J,dofCard);
  }
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
void
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
feIndicesLG(const UInt & elL, const UInt & I, const UInt & J, sVect<UInt> & indices)
{
  typedef dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
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

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
void
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
feIndicesGG(const UInt & elG, const UInt & I, const UInt & J, sVect<UInt> & indices)
{
  typedef dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(indices.size() == feNumBasisG(elG));
  assert(DOFMAPPER::feCards.isG(elG));
  
  //Alloc
  DOFCARD dofCard;
  UInt elL = DOFMAPPER::grid2d->getElements().getRowMapG(elG).getLid();
  
  //Filling the index vector
  for(UInt i=1; i <= feNumBasisL(elL); ++i)
  {
    dofCard = FETYPE::getDofCard(i);
    dofCard.setLocalElId(elL);
    indices(i) = DOFMAPPER::mapListG(I,J,dofCard);
  }
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
void
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
feEvalL(const UInt & elL, const UInt & I, const UInt & J, const point3d & Y, sVect<OUTTYPE> & val)
{
  typedef dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(elL >= 1);
  assert(elL <= DOFMAPPER::feCards.size());  
  assert(val.size() == feNumBasisL(elL));
  
  //Evaluation
  sVect<BASETYPE> basis(FETYPE::numBasis);
  FETYPE::setCards(DOFMAPPER::feCards.getDataL(elL), feeder.getCardLocal(elL));
  FETYPE::globalEval(Y,basis);
  
  DOFTYPE dof = traitsBasic<DOFTYPE>::getUnityIJ(I,J);
  
  for(UInt i=1; i <= feNumBasisL(elL); ++i)
  {
    val(i) = traitsMultiply<BASETYPE,DOFTYPE>::multiply(basis(i),dof);
  }
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
void
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
feEvalG(const UInt & elG, const UInt & I, const UInt & J, const point3d & Y, sVect<OUTTYPE> & val)
{
  typedef dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(val.size() == feNumBasisG(elG));
  assert(DOFMAPPER::feCards.isG(elG));
  
  //Alloc
  UInt elL = DOFMAPPER::grid2d->getElements().getRowMapG(elG).getLid();
  
  //Evaluation
  sVect<BASETYPE> basis(FETYPE::numBasis);
  FETYPE::setCards(DOFMAPPER::feCards.getDataL(elL), feeder.getCardLocal(elL));
  FETYPE::globalEval(Y,basis);
  
  DOFTYPE dof = traitsBasic<DOFTYPE>::getUnityIJ(I,J);
  
  for(UInt i=1; i <= feNumBasisL(elL); ++i)
  {
    val(i) = traitsMultiply<BASETYPE,DOFTYPE>::multiply(basis(i),dof);
  }
}

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
void
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
feGradL(const UInt & elL, const UInt & I, const UInt & J, const point3d & Y, sVect<OUTTYPE> & gradX, sVect<OUTTYPE> & gradY, sVect<OUTTYPE> & gradZ)
{
  typedef dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
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

template<typename FETYPE, typename DOFTYPE,  dms2d_order ORDER, dms2d_mode MODE>
void
feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>::
feGradG(const UInt & elG, const UInt & I, const UInt & J, const point3d & Y, sVect<OUTTYPE> & gradX, sVect<OUTTYPE> & gradY, sVect<OUTTYPE> & gradZ)
{
  typedef dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>  DOFMAPPER;
  
  //Assert
  assert(DOFMAPPER::feCards.isG(elG));
  assert(gradX.size() == feNumBasisG(elG));
  assert(gradY.size() == feNumBasisG(elG));
  assert(gradZ.size() == feNumBasisG(elG));
  
  //Alloc
  UInt elL = DOFMAPPER::grid2d->getElements().getRowMapG(elG).getLid();
  
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
