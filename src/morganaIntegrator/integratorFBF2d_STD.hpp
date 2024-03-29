/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTEGRATORFBF2D_STD_HPP
#define INTEGRATORFBF2D_STD_HPP

#include "intPolicySTD.hpp"
#include "elCard2d.hpp"
#include "intStaticCards.h"
#include "geoMapInterface.hpp"
#include "morganaIntegrator.hpp"
#include "morganaFiniteElements.hpp"


/*! Integrates on all the boundary edges of a 2d mesh */
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
class integratorFBF2d_STD
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename FUNCTIONAL::TEST_FECARD    TEST_FECARD;
    typedef typename FUNCTIONAL::TEST_GEOSHAPE  TEST_GEOSHAPE;
    typedef typename FUNCTIONAL::TEST_PMAPTYPE  TEST_PMAPTYPE;
    
    typedef typename FUNCTIONAL::INTGRID  INTGRID;
    typedef typename INTGRID::GRID_ELMAP  PMAPTYPE;
    
    typedef typename INTGRID::GEOSHAPE2D   GEOSHAPE2D;
    typedef typename GEOSHAPE2D::GEOBSHAPE GEOSHAPE1D;
    
    typedef elCard2d<GEOSHAPE2D,PMAPTYPE> ELCARD;
    typedef intPolicySTD                  INTCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const intClass intFlag = FBF2d_STD;
    static const OPType   opFlag  = fnBF2d;
    
    geoMapInterface<GEOSHAPE1D> geoInterface;
    //@}
    
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    integratorFBF2d_STD();
    
    /*! Constructor */
    integratorFBF2d_STD(const INTCARD & IntCard);
    
    /*! Set integration card */
    void setIntCard(const INTCARD & IntCard);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(FUNCTIONAL & fn);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(const Teuchos::RCP<FUNCTIONAL> & fn);
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
integratorFBF2d_STD<FUNCTIONAL,TYPE,PRECISION>::
integratorFBF2d_STD()
{
  assert(staticAssert<GEOSHAPE2D::nDim == 2>::returnValue);
}
    
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
integratorFBF2d_STD<FUNCTIONAL,TYPE,PRECISION>::
integratorFBF2d_STD(const INTCARD & IntCard)
{
  assert(staticAssert<GEOSHAPE2D::nDim == 2>::returnValue);
  assert( (&IntCard) == (&IntCard) );
}    

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
void
integratorFBF2d_STD<FUNCTIONAL,TYPE,PRECISION>::
setIntCard(const INTCARD & IntCard)
{
  assert( (&IntCard) == (&IntCard) );
}    

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorFBF2d_STD<FUNCTIONAL,TYPE,PRECISION>::
integration(FUNCTIONAL & fn)
{
  typedef intStaticCard<GEOSHAPE1D,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = fn.numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  point3d Vx;
  
  //Integration loop
  for(UInt j=1; j <= GEOSHAPE2D::numEdges; ++j)
  {
    fn.setLocEdge(j);
    
    if(fn.isBoundary())
    {
      for(UInt i=1; i <= INTSTATIC::N; ++i)
      {
        //Evaluate field
        fn.eval(INTSTATIC::getYn(i), bufMat);
      
        //Evaluate volume element
        Vx = geoInterface.getDerX(fn.getGlobEdgePoints(), INTSTATIC::getYn(i));
  
        //Sum
        for(UInt k=1; k <= matSize; ++k)
        { mat(k) += bufMat(k) * Vx.norm2() * INTSTATIC::getWn(i);}
      }
    }
  }
  
  return(mat);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorFBF2d_STD<FUNCTIONAL,TYPE,PRECISION>::
integration(const Teuchos::RCP<FUNCTIONAL> & fn)
{
  typedef intStaticCard<GEOSHAPE1D,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = fn->numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  point3d Vx;
  
  //Integration loop
  for(UInt j=1; j <= GEOSHAPE2D::numEdges; ++j)
  {
    fn->setLocEdge(j);
    
    if(fn->isBoundary())
    {
      for(UInt i=1; i <= INTSTATIC::N; ++i)
      {
        //Evaluate field
        fn->eval(INTSTATIC::getYn(i), bufMat);
      
        //Evaluate volume element
        Vx = geoInterface.getDerX(fn->getGlobEdgePoints(), INTSTATIC::getYn(i));
  
        //Sum
        for(UInt k=1; k <= matSize; ++k)
        { mat(k) += bufMat(k) * Vx.norm2() * INTSTATIC::getWn(i); }
      }
    }
  }
  
  return(mat);
}


#endif
