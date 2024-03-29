/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTEGRATOROEL2D_SKT_HPP
#define INTEGRATOROEL2D_SKT_HPP

#include "glRule.h"

#include "intPolicySTD.hpp"
#include "elCard2d.hpp"
#include "intStaticCards.h"
#include "geoMapInterface.hpp"

#include "morganaIntegrator.hpp"
#include "morganaFiniteElements.hpp"


/*! Spectral 2d integrator */
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
class integratorOEL2d_SKT
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
    typedef elCard2d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef intPolicySTD                INTCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const intClass intFlag = OEL2d_SKT;
    static const OPType   opFlag  = opEL2d;
    
    geoMapInterface<GEOSHAPE> geoInterface;
    //@}
    
    /*! @name Functions */ //@{
  public:
    integratorOEL2d_SKT();
    
    /*! Constructor */
    integratorOEL2d_SKT(const INTCARD & IntCard);
    
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
integratorOEL2d_SKT<OPERATOR,TYPE,PRECISION>::
integratorOEL2d_SKT()
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 2>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 2>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName      == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
integratorOEL2d_SKT<OPERATOR,TYPE,PRECISION>::
integratorOEL2d_SKT(const INTCARD & IntCard)
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 2>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 2>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName      == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  assert( (&IntCard) == (&IntCard) );
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
void
integratorOEL2d_SKT<OPERATOR,TYPE,PRECISION>::
setIntCard(const INTCARD & IntCard)
{
  assert( (&IntCard) == (&IntCard) );
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorOEL2d_SKT<OPERATOR,TYPE,PRECISION>::
integration(OPERATOR & op)
{
  //Degree determination
  UInt rx = PRECISION * op.getFeCardL_field().getRx() * op.getFeCardL_test().getRx() * UInt(op.getFeCardL_field().getIsActive()) * UInt(op.getFeCardL_test().getIsActive());
  UInt ry = PRECISION * op.getFeCardL_field().getRy() * op.getFeCardL_test().getRy() * UInt(op.getFeCardL_field().getIsActive()) * UInt(op.getFeCardL_test().getIsActive());
  
  //Set toll
  Real toll = 1e-10;
  
  //Nodes
  sVect<Real>  xgl,  ygl;
  sVect<Real> wxgl, wygl;
  
  glRule ruleX( UInt((rx + 1) / 2.0) + 1 ); 
  glRule ruleY( UInt((ry + 1) / 2.0) + 1 );
  
  ruleX.compute(xgl,wxgl,toll);
  ruleY.compute(ygl,wygl,toll);
  
  
  //Alloc
  UInt matSize = op.numIndex_field() * op.numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  point3d Vx, Vy, N, Y;
  
  
  //Integration loop
  for(UInt ix = 1; ix <= xgl.size(); ++ix)
  {
    Y.setX(xgl(ix));
    
    for(UInt iy = 1; iy <= ygl.size(); ++iy)
    {
      Y.setY(ygl(iy));
      
      //Evaluate field
      op.eval(Y, bufMat);
      
      //Evaluate volume element
      Vx = geoInterface.getDerX(op.getElCard().getNodes(), Y);
      Vy = geoInterface.getDerY(op.getElCard().getNodes(), Y);
      N  = Vx ^ Vy;
      
      //Sum
      for(UInt j=1; j <= matSize; ++j)
      {
        mat(j) += bufMat(j) * N.norm2() * wxgl(ix) * wygl(iy);
      }
    }
  }
  
  return(mat);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorOEL2d_SKT<OPERATOR,TYPE,PRECISION>::
integration(const Teuchos::RCP<OPERATOR> & op)
{
  //Degree determination
  UInt rx = PRECISION * op->getFeCardL_field().getRx() * op->getFeCardL_test().getRx() * UInt(op->getFeCardL_field().getIsActive()) * UInt(op->getFeCardL_test().getIsActive());
  UInt ry = PRECISION * op->getFeCardL_field().getRy() * op->getFeCardL_test().getRy() * UInt(op->getFeCardL_field().getIsActive()) * UInt(op->getFeCardL_test().getIsActive());
  
  //Set toll
  Real toll = 1e-10;
  
  //Nodes
  sVect<Real>  xgl,  ygl;
  sVect<Real> wxgl, wygl;
  
  glRule ruleX( UInt((rx + 1) / 2.0) + 1 ); 
  glRule ruleY( UInt((ry + 1) / 2.0) + 1 );
  
  ruleX.compute(xgl,wxgl,toll);
  ruleY.compute(ygl,wygl,toll);
  
  
  //Alloc
  UInt matSize = op->numIndex_field() * op->numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  point3d Vx, Vy, N, Y;
  
  
  //Integration loop
  for(UInt ix = 1; ix <= xgl.size(); ++ix)
  {
    Y.setX(xgl(ix));
    
    for(UInt iy = 1; iy <= ygl.size(); ++iy)
    {
      Y.setY(ygl(iy));
      
      //Evaluate field
      op->eval(Y, bufMat);
      
      //Evaluate volume element
      Vx = geoInterface.getDerX(op->getElCard().getNodes(), Y);
      Vy = geoInterface.getDerY(op->getElCard().getNodes(), Y);
      N  = Vx ^ Vy;
      
      //Sum
      for(UInt j=1; j <= matSize; ++j)
      {
        mat(j) += bufMat(j) * N.norm2() * wxgl(ix) * wygl(iy);
      }
    }
  }
  
  return(mat);
}

#endif
