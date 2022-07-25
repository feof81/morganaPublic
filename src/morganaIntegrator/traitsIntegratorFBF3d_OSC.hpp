/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef TRAITSINTEGRATORFBF3D_OSC_HPP
#define TRAITSINTEGRATORFBF3D_OSC_HPP

#include "morganaFiniteElements.hpp"

#include "geoShapes.h"
#include "geoMapInterface.hpp"
#include "pointElement.hpp"
#include "elCard3d.hpp"

#include "feOsc3d_card.h"

#include "morganaFiniteElements.hpp"
#include "morganaIntegrator.hpp"
#include "intOscSymplexSupport.h"
#include "intStaticCards.h"
#include "intPolicySTD.hpp"


//_________________________________________________________________________________________________
// EMPTY TRAIT
//-------------------------------------------------------------------------------------------------
template<typename FUNCTIONAL, typename GEOSHAPE, intTypes TYPE, UInt PRECISION>
class traitsIntegratorFBF3d_OSC
{ };


//_________________________________________________________________________________________________
// TETRA
//-------------------------------------------------------------------------------------------------
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
class traitsIntegratorFBF3d_OSC<FUNCTIONAL,linearTriangle,TYPE,PRECISION>
{
    /*! @name Typedefs */ //@{
  public:    
    typedef typename FUNCTIONAL::TEST_FECARD   TEST_FECARD;
    typedef typename FUNCTIONAL::TEST_GEOSHAPE TEST_GEOSHAPE;
    typedef typename FUNCTIONAL::TEST_PMAPTYPE TEST_PMAPTYPE;
  
    typedef typename FUNCTIONAL::INTGRID   INTGRID;
    typedef typename INTGRID::GRID_ELMAP  PMAPTYPE;
    typedef TEST_GEOSHAPE                 GEOSHAPE3D;
    typedef typename INTGRID::GEOSHAPE2D  GEOSHAPE2D;
    
    typedef elCard2d<GEOSHAPE2D,PMAPTYPE> ELCARD;
    typedef intPolicySTD                  INTCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const intClass intFlag = FBF3d_OSC;
    static const OPType   opFlag  = fnBF3d;
    
    sVect<pointElement<GEOSHAPE2D> > subElements;
    geoMapInterface<GEOSHAPE2D>   geoInterface2d;
    geoMapInterface<GEOSHAPE3D>   geoInterface3d;
    //@}
    
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    traitsIntegratorFBF3d_OSC();
    
    /*! Constructor */
    traitsIntegratorFBF3d_OSC(const INTCARD & IntCard);
    
    /*! Set integration card */
    void setIntCard(const INTCARD & IntCard);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(FUNCTIONAL & op);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(const Teuchos::RCP<FUNCTIONAL> & op);
    //@}
};


template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
traitsIntegratorFBF3d_OSC<FUNCTIONAL,linearTriangle,TYPE,PRECISION>::
traitsIntegratorFBF3d_OSC()
{
  //Assert
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
  assert(staticAssert<TEST_FECARD::cardLabel == CR_OS_3d>::returnValue);
  
  //Build sub elements
  pointElement<GEOSHAPE2D> elem;
  
  for(UInt i=1; i <= 3; ++i)
  { elem.setPoint(i,linearTriangle::getRefNodes(i)); }
  
  subElements = intOscSymplexSupport::refineTriangle(PRECISION -1,elem);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
traitsIntegratorFBF3d_OSC<FUNCTIONAL,linearTriangle,TYPE,PRECISION>::
traitsIntegratorFBF3d_OSC(const INTCARD & IntCard)
{
  //Assert
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
  assert(staticAssert<TEST_FECARD::cardLabel == CR_OS_3d>::returnValue);
  
  //Build sub elements
  pointElement<GEOSHAPE2D> elem;
  
  for(UInt i=1; i <= 3; ++i)
  { elem.setPoint(i,linearTriangle::getRefNodes(i)); }
  
  subElements = intOscSymplexSupport::refineTriangle(PRECISION -1,elem);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
void
traitsIntegratorFBF3d_OSC<FUNCTIONAL,linearTriangle,TYPE,PRECISION>::
setIntCard(const INTCARD & IntCard)
{
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
sVect<Real>
traitsIntegratorFBF3d_OSC<FUNCTIONAL,linearTriangle,TYPE,PRECISION>::
integration(FUNCTIONAL & op)
{
  //Alloc
  point3d P1, P2, P3, YM, XM;
  komplex C;
  
  //Extracted data
  UInt el3d    = op.getEl3d();
  UInt matSize = op.numIndex_test();
  
  //Card management
  TEST_FECARD   testCard = op.getFeCardL_test();
  sVect<point3d> nodes3d = op.getIntegrationGrid()->getElementNodesL(el3d);
  
  //Compute the global K and X
  sVect<point3d> Ktest(testCard.size());
  
  tensor3d T3d;
  T3d.setCol(1, nodes3d(2) - nodes3d(1));
  T3d.setCol(2, nodes3d(3) - nodes3d(1));
  T3d.setCol(3, nodes3d(4) - nodes3d(1));
  T3d.computeInverse();
  
  for(UInt i=1; i <= testCard.size(); ++i)
  { Ktest(i)  = T3d.firstIndexSaturation(testCard.getH(i)); }
  
  //Compute the integral
  sVect<komplex> matC(matSize);  
  sVect<Real>     mat(matSize);
  
  UInt tot;
  
  for(UInt k=1; k <= subElements.size(); ++k)
  {
    //Local nodes
    P1 = geoInterface2d.getPosition(op.getGlobFacePoints(), subElements(k).getPoint(1));
    P2 = geoInterface2d.getPosition(op.getGlobFacePoints(), subElements(k).getPoint(2));
    P3 = geoInterface2d.getPosition(op.getGlobFacePoints(), subElements(k).getPoint(3));
    
    YM = (subElements(k).getPoint(1) 
        + subElements(k).getPoint(2)
        + subElements(k).getPoint(3)) / 3.0;

    XM = (P1 + P2 + P3) / 3.0;
    
    //Evaluate field
    op.eval(YM,matC);
    
    for(UInt j=1; j <= op.numIndex_test(); ++j)
    {
      C = intOscSymplexSupport::intS0_2d(Ktest(j), P1,P2,P3)
        * komplex::iexp(-point3d::dot(Ktest(j), XM))
        * matC(j);

      mat(j) += komplex::cComp(C,op.I_test);
    }
  }

  return(mat);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
sVect<Real>
traitsIntegratorFBF3d_OSC<FUNCTIONAL,linearTriangle,TYPE,PRECISION>::
integration(const Teuchos::RCP<FUNCTIONAL> & op)
{
  return(integration(*op));
}

#endif
