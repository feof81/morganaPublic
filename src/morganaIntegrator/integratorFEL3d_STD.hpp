/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTEGRATORFEL3D_STD_HPP
#define INTEGRATORFEL3D_STD_HPP

#include "intPolicySTD.hpp"
#include "elCard3d.hpp"
#include "intStaticCards.h"
#include "geoMapInterface.hpp"

#include "morganaIntegrator.hpp"
#include "morganaFiniteElements.hpp"


/*! Standard 3d integrator (functional) */
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
class integratorFEL3d_STD
{
    /*! @name Typedefs */ //@{
  public:    
    typedef typename FUNCTIONAL::TEST_FECARD       TEST_FECARD;
    typedef typename FUNCTIONAL::TEST_GEOSHAPE     TEST_GEOSHAPE;
    typedef typename FUNCTIONAL::TEST_PMAPTYPE     TEST_PMAPTYPE;
    
    typedef TEST_GEOSHAPE               GEOSHAPE;
    typedef TEST_PMAPTYPE               PMAPTYPE;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef intPolicySTD                INTCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const intClass intFlag = FEL3d_STD;
    static const OPType   opFlag  = fnEL3d;
    
    geoMapInterface<GEOSHAPE> geoInterface;
    //@}
    
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    integratorFEL3d_STD();
    
    /*! Constructor */
    integratorFEL3d_STD(const INTCARD & IntCard);
    
    /*! Set integration card */
    void setIntCard(const INTCARD & IntCard);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(FUNCTIONAL & op);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(const Teuchos::RCP<FUNCTIONAL> & op);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
integratorFEL3d_STD<FUNCTIONAL,TYPE,PRECISION>::
integratorFEL3d_STD()
{
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
integratorFEL3d_STD<FUNCTIONAL,TYPE,PRECISION>::
integratorFEL3d_STD(const INTCARD & IntCard)
{
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
void
integratorFEL3d_STD<FUNCTIONAL,TYPE,PRECISION>::
setIntCard(const INTCARD & IntCard)
{
  assert( (&IntCard) == (&IntCard) );
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorFEL3d_STD<FUNCTIONAL,TYPE,PRECISION>::
integration(FUNCTIONAL & op)
{
  typedef intStaticCard<GEOSHAPE,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = op.numIndex_test();
  sVect<Real> bufMat(matSize);
  Real jac;
  
  sVect<Real> mat(matSize);
  
  //Integration loop
  for(UInt i=1; i <= INTSTATIC::N; ++i)
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

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorFEL3d_STD<FUNCTIONAL,TYPE,PRECISION>::
integration(const Teuchos::RCP<FUNCTIONAL> & op)
{
  typedef intStaticCard<GEOSHAPE,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = op->numIndex_test();
  sVect<Real> bufMat(matSize);
  Real jac;
  
  sVect<Real> mat(matSize);
  
  //Integration loop
  for(UInt i=1; i <= INTSTATIC::N; ++i)
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
