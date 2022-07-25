/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTEGRATOROEL3D_STC_HPP
#define INTEGRATOROEL3D_STC_HPP

#include "intPolicySTC.h"
#include "elCard3d.hpp"
#include "intStaticCards.h"
#include "geoMapInterface.hpp"

#include "morganaIntegrator.hpp"
#include "morganaFiniteElements.hpp"


/*! Conditional 3d integrator, some elements could be de-activated (operator) */
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
class integratorOEL3d_STC
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename OPERATOR::FIELD_FECARD      FIELD_FECARD;
    typedef typename OPERATOR::FIELD_GEOSHAPE    FIELD_GEOSHAPE;
    typedef typename OPERATOR::FIELD_PMAPTYPE    FIELD_PMAPTYPE;
    
    typedef typename OPERATOR::TEST_FECARD       TEST_FECARD;
    typedef typename OPERATOR::TEST_GEOSHAPE     TEST_GEOSHAPE;
    typedef typename OPERATOR::TEST_PMAPTYPE     TEST_PMAPTYPE;
    
    typedef TEST_GEOSHAPE               GEOSHAPE;
    typedef TEST_PMAPTYPE               PMAPTYPE;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef intPolicySTC                INTCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const intClass intFlag = OEL3d_STC;
    static const OPType   opFlag  = opEL3d;
    
    INTCARD intCard;
    
    geoMapInterface<GEOSHAPE> geoInterface;
    //@}
    
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    integratorOEL3d_STC();
    
    /*! Constructor */
    integratorOEL3d_STC(const INTCARD & IntCard);
    
    /*! Set integration card */
    void setIntCard(const INTCARD & IntCard);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(OPERATOR & op);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(const Teuchos::RCP<OPERATOR> & op);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
integratorOEL3d_STC<OPERATOR,TYPE,PRECISION>::
integratorOEL3d_STC()
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 3>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName      == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
integratorOEL3d_STC<OPERATOR,TYPE,PRECISION>::
integratorOEL3d_STC(const INTCARD & IntCard)
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 3>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName      == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  intCard = IntCard;
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
void
integratorOEL3d_STC<OPERATOR,TYPE,PRECISION>::
setIntCard(const INTCARD & IntCard)
{
  intCard = IntCard;
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorOEL3d_STC<OPERATOR,TYPE,PRECISION>::
integration(OPERATOR & op)
{
  typedef intStaticCard<GEOSHAPE,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = op.numIndex_field() * op.numIndex_test();
  sVect<Real> bufMat(matSize);
  Real jac;
  
  sVect<Real> mat(matSize);
  
  //Integration loop
  for(UInt i=1; i <= (INTSTATIC::N * UInt(intCard.getIsActive())); ++i)
  {
    //Evaluate field
    op.eval(INTSTATIC::getYn(i), bufMat);
    
    //Evaluate jacobian
    jac = geoInterface.getGradientDet(op.getElCard().getNodes(), INTSTATIC::getYn(i));

    //Sum
    for(UInt j=1; j <= matSize; ++j)
    {
      mat(j) += bufMat(j) * jac * INTSTATIC::getWn(i);
    }
  }
  
  return(mat);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorOEL3d_STC<OPERATOR,TYPE,PRECISION>::
integration(const Teuchos::RCP<OPERATOR> & op)
{
  typedef intStaticCard<GEOSHAPE,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = op->numIndex_field() * op->numIndex_test();
  sVect<Real> bufMat(matSize);
  Real jac;
  
  sVect<Real> mat(matSize);
  
  //Integration loop
  for(UInt i=1; i <= (INTSTATIC::N * UInt(intCard.getIsActive())); ++i)
  {
    //Evaluate field
    op->eval(INTSTATIC::getYn(i), bufMat);
    
    //Evaluate jacobian
    jac = geoInterface.getGradientDet(op->getElCard().getNodes(), INTSTATIC::getYn(i));    
    
    //Sum
    for(UInt j=1; j <= matSize; ++j)
    {
      mat(j) += bufMat(j) * jac * INTSTATIC::getWn(i);
    }
  }
  
  return(mat);
}



#endif
