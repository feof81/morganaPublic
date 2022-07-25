/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTEGRATOROEL1D_STD_HPP
#define INTEGRATOROEL1D_STD_HPP

#include "intPolicySTD.hpp"
#include "elCard1d.hpp"
#include "intStaticCards.h"
#include "geoMapInterface.hpp"

#include "morganaIntegrator.hpp"
#include "morganaFiniteElements.hpp"


/*! Standard 3d integrator (operator) */
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
class integratorOEL1d_STD
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename OPERATOR::FIELD_PMAPTYPE    FIELD_PMAPTYPE;
    typedef typename OPERATOR::TEST_PMAPTYPE     TEST_PMAPTYPE;
    
    typedef typename OPERATOR::INTGRID    INTGRID;
    typedef typename INTGRID::GRID_ELMAP  PMAPTYPE;
    typedef typename INTGRID::GEOSHAPE1D  GEOSHAPE;
    
    typedef elCard1d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef intPolicySTD                INTCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const intClass intFlag = OEL1d_STD;
    static const OPType   opFlag  = opEL1d;
    
    geoMapInterface<GEOSHAPE> geoInterface;
    //@}
    
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    integratorOEL1d_STD();
    
    /*! Constructor */
    integratorOEL1d_STD(const INTCARD & IntCard);
    
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
integratorOEL1d_STD<OPERATOR,TYPE,PRECISION>::
integratorOEL1d_STD()
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
integratorOEL1d_STD<OPERATOR,TYPE,PRECISION>::
integratorOEL1d_STD(const INTCARD & IntCard)
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
void
integratorOEL1d_STD<OPERATOR,TYPE,PRECISION>::
setIntCard(const INTCARD & IntCard)
{
  assert( (&IntCard) == (&IntCard) );
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorOEL1d_STD<OPERATOR,TYPE,PRECISION>::
integration(OPERATOR & op)
{
  typedef intStaticCard<GEOSHAPE,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = op.numIndex_field() * op.numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  point3d Vx;
  
  //Integration loop
  for(UInt i=1; i <= INTSTATIC::N; ++i)
  {
    //Evaluate field
    op.eval(INTSTATIC::getYn(i), bufMat);
    
    //Evaluate volume element
    Vx = geoInterface.getDerX(op.getElCard().getNodes(), INTSTATIC::getYn(i));
  
    //Sum
    for(UInt j=1; j <= matSize; ++j)
    {
      mat(j) += bufMat(j) * Vx.norm2() * INTSTATIC::getWn(i);
    }    
  }
  
  return(mat);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorOEL1d_STD<OPERATOR,TYPE,PRECISION>::
integration(const Teuchos::RCP<OPERATOR> & op)
{
  typedef intStaticCard<GEOSHAPE,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = op->numIndex_field() * op->numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  point3d Vx;
  
  //Integration loop
  for(UInt i=1; i <= INTSTATIC::N; ++i)
  {
    //Evaluate field
    op->eval(INTSTATIC::getYn(i), bufMat);
    
    //Evaluate volume element
    Vx = geoInterface.getDerX(op->getElCard().getNodes(), INTSTATIC::getYn(i));
  
    //Sum
    for(UInt j=1; j <= matSize; ++j)
    {
      mat(j) += bufMat(j) * Vx.norm2() * INTSTATIC::getWn(i);
    }    
  }
  
  return(mat);
}


#endif
