/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTEGRATORFLA3D_STC_HPP
#define INTEGRATORFLA3D_STC_HPP

#include "intPolicySTC.h"
#include "elCard3d.hpp"
#include "intStaticCards.h"
#include "geoMapInterface.hpp"
#include "morganaIntegrator.hpp"
#include "morganaFiniteElements.hpp"


/*! Integrates on a single boundary face of a 3d mesh (face is pre-set) */
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
class integratorFLA3d_STC
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename FUNCTIONAL::TEST_FECARD    TEST_FECARD;
    typedef typename FUNCTIONAL::TEST_GEOSHAPE  TEST_GEOSHAPE;
    typedef typename FUNCTIONAL::TEST_PMAPTYPE  TEST_PMAPTYPE;
    
    typedef typename FUNCTIONAL::INTGRID  INTGRID;
    typedef typename INTGRID::GRID_ELMAP  PMAPTYPE;
    
    typedef typename INTGRID::GEOSHAPE3D   GEOSHAPE3D;
    typedef typename GEOSHAPE3D::GEOBSHAPE GEOSHAPE2D;
    
    typedef elCard3d<GEOSHAPE3D,PMAPTYPE> ELCARD;
    typedef intPolicySTC                  INTCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const intClass intFlag = FLA3d_STC;
    static const OPType   opFlag  = fnLA3d;
    
    INTCARD intCard;
    
    geoMapInterface<GEOSHAPE2D> geoInterface;
    //@}
    
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    integratorFLA3d_STC();
    
    /*! Constructor */
    integratorFLA3d_STC(const INTCARD & IntCard);
    
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
integratorFLA3d_STC<FUNCTIONAL,TYPE,PRECISION>::
integratorFLA3d_STC()
{
  assert(staticAssert<GEOSHAPE3D::nDim == 3>::returnValue);
}
    
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
integratorFLA3d_STC<FUNCTIONAL,TYPE,PRECISION>::
integratorFLA3d_STC(const INTCARD & IntCard)
{
  assert(staticAssert<GEOSHAPE3D::nDim == 3>::returnValue);
  intCard = IntCard;
}    

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
void
integratorFLA3d_STC<FUNCTIONAL,TYPE,PRECISION>::
setIntCard(const INTCARD & IntCard)
{
  intCard = IntCard;
}    

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorFLA3d_STC<FUNCTIONAL,TYPE,PRECISION>::
integration(FUNCTIONAL & fn)
{
  typedef intStaticCard<GEOSHAPE2D,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = fn.numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  point3d Vx, Vy, N;
  
  //Integration loop
  for(UInt i=1; i <= (INTSTATIC::N * UInt(intCard.getIsActive())); ++i)
  {
    //Evaluate field
    fn.eval(INTSTATIC::getYn(i), bufMat);
      
    //Evaluate volume element
    Vx = geoInterface.getDerX(fn.getGlobFacePoints(), INTSTATIC::getYn(i));
    Vy = geoInterface.getDerY(fn.getGlobFacePoints(), INTSTATIC::getYn(i));
    N  = Vx ^ Vy; 
  
    //Sum
    for(UInt j=1; j <= matSize; ++j)
    { mat(j) += bufMat(j) * N.norm2() * INTSTATIC::getWn(i);}
  }
  
  return(mat);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorFLA3d_STC<FUNCTIONAL,TYPE,PRECISION>::
integration(const Teuchos::RCP<FUNCTIONAL> & fn)
{
  typedef intStaticCard<GEOSHAPE2D,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = fn->numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  point3d Vx, Vy, N;
  
  //Integration loop
  for(UInt i=1; i <= (INTSTATIC::N * UInt(intCard.getIsActive())); ++i)
  {
    //Evaluate field
    fn->eval(INTSTATIC::getYn(i), bufMat);
      
    //Evaluate volume element
    Vx = geoInterface.getDerX(fn->getGlobFacePoints(), INTSTATIC::getYn(i));
    Vy = geoInterface.getDerY(fn->getGlobFacePoints(), INTSTATIC::getYn(i));
    N  = Vx ^ Vy; 
  
    //Sum
    for(UInt j=1; j <= matSize; ++j)
    { mat(j) += bufMat(j) * N.norm2() * INTSTATIC::getWn(i); }
  }

  return(mat);
}


#endif
