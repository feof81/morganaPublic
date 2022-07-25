/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTEGRATOROHY2D_STD_HPP
#define INTEGRATOROHY2D_STD_HPP

#include "intPolicySTD.hpp"
#include "elCard2d.hpp"
#include "intStaticCards.h"
#include "geoMapInterface.hpp"

#include "morganaIntegrator.hpp"
#include "morganaFiniteElements.hpp"



/*! Hybrid 3d integrator, integration on the internal faces (operator) */
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
class integratorOHY2d_STD
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
    
    typedef typename INTGRID::GEOSHAPE2D   GEOSHAPE2D;
    typedef typename GEOSHAPE2D::GEOBSHAPE GEOSHAPE1D;
    
    typedef elCard2d<GEOSHAPE2D,PMAPTYPE> ELCARD;
    typedef intPolicySTD                  INTCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const intClass intFlag = OHY2d_STD;
    static const OPType   opFlag  = opHY2d;
    
    geoMapInterface<GEOSHAPE1D> geoInterface;
    //@}
    
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    integratorOHY2d_STD();
    
    /*! Constructor */
    integratorOHY2d_STD(const INTCARD & IntCard);
    
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
integratorOHY2d_STD<OPERATOR,TYPE,PRECISION>::
integratorOHY2d_STD()
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 2>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 2>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName      == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  assert(staticAssert<TYPE == HYBRID>::returnValue);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
integratorOHY2d_STD<OPERATOR,TYPE,PRECISION>::
integratorOHY2d_STD(const INTCARD & IntCard)
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 2>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 2>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName      == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  assert(staticAssert<TYPE == HYBRID>::returnValue);
  
  assert( (&IntCard) == (&IntCard) );
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
void
integratorOHY2d_STD<OPERATOR,TYPE,PRECISION>::
setIntCard(const INTCARD & IntCard)
{
  assert( (&IntCard) == (&IntCard) );
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorOHY2d_STD<OPERATOR,TYPE,PRECISION>::
integration(OPERATOR & op)
{
  typedef intStaticCard<GEOSHAPE2D,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = op.numIndex_field() * op.numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  point3d Vx;
  
  //Integration loop
  for(UInt j=1; j <= GEOSHAPE2D::numEdges; ++j)
  {
    op.setLocEdge(j);
    
    for(UInt i=1; i <= INTSTATIC::N; ++i)
    {
      //Evaluate field
      op.eval(INTSTATIC::getYn(i), bufMat);
      
      //Evaluate volume element
      Vx = geoInterface.getDerX(op.getGlobEdgePoints(), INTSTATIC::getYn(i)); 
  
      //Sum
      for(UInt k=1; k <= matSize; ++k)
      {
        mat(k) += bufMat(k) * Vx.norm2() * INTSTATIC::getWn(i);
      }
    }
  }
  
  return(mat);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorOHY2d_STD<OPERATOR,TYPE,PRECISION>::
integration(const Teuchos::RCP<OPERATOR> & op)
{
  typedef intStaticCard<GEOSHAPE2D,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = op->numIndex_field() * op->numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  point3d Vx;
  
  //Integration loop
  for(UInt j=1; j <= GEOSHAPE2D::numEdges; ++j)
  {
    op->setLocEdge(j);
    
    for(UInt i=1; i <= INTSTATIC::N; ++i)
    {
      //Evaluate field
      op->eval(INTSTATIC::getYn(i), bufMat);
      
      //Evaluate volume element
      Vx = geoInterface.getDerX(op->getGlobEdgePoints(), INTSTATIC::getYn(i));
  
      //Sum
      for(UInt k=1; k <= matSize; ++k)
      {
        mat(k) += bufMat(k) * Vx.norm2() * INTSTATIC::getWn(i);
      }
    }
  }
  
  return(mat);
}

#endif
