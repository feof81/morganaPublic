/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTEGRATOROHY3D_STC_HPP
#define INTEGRATOROHY3D_STC_HPP

#include "intPolicySTC.h"
#include "elCard3d.hpp"
#include "intStaticCards.h"
#include "geoMapInterface.hpp"

#include "morganaIntegrator.hpp"
#include "morganaFiniteElements.hpp"



/*! Hybrid 3d integrator, integration on the internal faces (operator) */
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
class integratorOHY3d_STC
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename OPERATOR::FIELD_FECARD      FIELD_FECARD;
    typedef typename OPERATOR::FIELD_GEOSHAPE    FIELD_GEOSHAPE;
    typedef typename OPERATOR::FIELD_PMAPTYPE    FIELD_PMAPTYPE;
    
    typedef typename OPERATOR::TEST_FECARD       TEST_FECARD;
    typedef typename OPERATOR::TEST_GEOSHAPE     TEST_GEOSHAPE;
    typedef typename OPERATOR::TEST_PMAPTYPE     TEST_PMAPTYPE;
    
    typedef typename OPERATOR::INTGRID    INTGRID;
    typedef typename INTGRID::GRID_ELMAP  PMAPTYPE;
    
    typedef typename INTGRID::GEOSHAPE3D   GEOSHAPE3D;
    typedef typename GEOSHAPE3D::GEOBSHAPE GEOSHAPE2D;
    
    typedef elCard3d<GEOSHAPE3D,PMAPTYPE> ELCARD;
    typedef intPolicySTC                  INTCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const intClass intFlag = OHY3d_STC;
    static const OPType   opFlag  = opHY3d;
    
    INTCARD intCard;
    
    geoMapInterface<GEOSHAPE2D> geoInterface;
    //@}
    
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    integratorOHY3d_STC();
    
    /*! Constructor */
    integratorOHY3d_STC(const INTCARD & IntCard);
    
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
integratorOHY3d_STC<OPERATOR,TYPE,PRECISION>::
integratorOHY3d_STC()
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 3>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName      == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  assert(staticAssert<TYPE == HYBRID>::returnValue);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
integratorOHY3d_STC<OPERATOR,TYPE,PRECISION>::
integratorOHY3d_STC(const INTCARD & IntCard)
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 3>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName      == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  assert(staticAssert<TYPE == HYBRID>::returnValue);
  
  intCard = IntCard;
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
void
integratorOHY3d_STC<OPERATOR,TYPE,PRECISION>::
setIntCard(const INTCARD & IntCard)
{
  intCard = IntCard;
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorOHY3d_STC<OPERATOR,TYPE,PRECISION>::
integration(OPERATOR & op)
{
  typedef intStaticCard<GEOSHAPE2D,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = op.numIndex_field() * op.numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  point3d Vx, Vy, N;
  
  //Integration loop
  for(UInt j=1; j <= GEOSHAPE3D::numFaces; ++j)
  {
    op.setLocFace(j);
    
    for(UInt i=1; i <= (INTSTATIC::N * UInt(intCard.getIsActive())); ++i)
    {
      //Evaluate field
      op.eval(INTSTATIC::getYn(i), bufMat);
      
      //Evaluate volume element
      Vx = geoInterface.getDerX(op.getGlobFacePoints(), INTSTATIC::getYn(i));
      Vy = geoInterface.getDerY(op.getGlobFacePoints(), INTSTATIC::getYn(i));
      N  = Vx ^ Vy; 
  
      //Sum
      for(UInt k=1; k <= matSize; ++k)
      {
        mat(k) += bufMat(k) * N.norm2() * INTSTATIC::getWn(i);
      }
    }
  }
  
  return(mat);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorOHY3d_STC<OPERATOR,TYPE,PRECISION>::
integration(const Teuchos::RCP<OPERATOR> & op)
{
  typedef intStaticCard<GEOSHAPE2D,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = op->numIndex_field() * op->numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  point3d Vx, Vy, N;
  
  //Integration loop
  for(UInt j=1; j <= GEOSHAPE3D::numFaces; ++j)
  {
    op->setLocFace(j);
    
    for(UInt i=1; i <= (INTSTATIC::N * UInt(intCard.getIsActive())); ++i)
    {
      //Evaluate field
      op->eval(INTSTATIC::getYn(i), bufMat);
      
      //Evaluate volume element
      Vx = geoInterface.getDerX(op->getGlobFacePoints(), INTSTATIC::getYn(i));
      Vy = geoInterface.getDerY(op->getGlobFacePoints(), INTSTATIC::getYn(i));
      N  = Vx ^ Vy; 
  
      //Sum
      for(UInt k=1; k <= matSize; ++k)
      {
        mat(k) += bufMat(k) * N.norm2() * INTSTATIC::getWn(i);
      }
    }
  }  
  
  return(mat);
}



#endif
