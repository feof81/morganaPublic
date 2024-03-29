/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef TRAITSINTEGRATORFEL3D_OSC_HPP
#define TRAITSINTEGRATORFEL3D_OSC_HPP

#include "morganaFiniteElements.hpp"

#include "geoShapes.h"
#include "geoMapInterface.hpp"
#include "pointElement.hpp"
#include "elCard3d.hpp"

#include "feOsc3d_card.h"

#include "morganaIntegrator.hpp"
#include "intOscSymplexSupport.h"
#include "intStaticCards.h"
#include "intPolicySTD.hpp"


//_________________________________________________________________________________________________
// EMPTY TRAIT
//-------------------------------------------------------------------------------------------------
template<typename FUNCTIONAL, typename GEOSHAPE, intTypes TYPE, UInt PRECISION>
class traitsIntegratorFEL3d_OSC
{ };


//_________________________________________________________________________________________________
// TETRA
//-------------------------------------------------------------------------------------------------
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
class traitsIntegratorFEL3d_OSC<FUNCTIONAL,linearTetra,TYPE,PRECISION>
{
    /*! @name Typedefs */ //@{
  public:    
    typedef typename FUNCTIONAL::TEST_FECARD    TEST_FECARD;
    typedef typename FUNCTIONAL::TEST_GEOSHAPE  TEST_GEOSHAPE;
    typedef typename FUNCTIONAL::TEST_PMAPTYPE  TEST_PMAPTYPE;
    
    typedef TEST_GEOSHAPE               GEOSHAPE;
    typedef TEST_PMAPTYPE               PMAPTYPE;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef intPolicySTD                INTCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const intClass intFlag = FEL3d_OSC;
    static const OPType   opFlag  = fnEL3d;
    
    sVect<pointElement<GEOSHAPE> > subElements;
    geoMapInterface<GEOSHAPE>     geoInterface;
    //@}
    
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    traitsIntegratorFEL3d_OSC();
    
    /*! Constructor */
    traitsIntegratorFEL3d_OSC(const INTCARD & IntCard);
    
    /*! Set integration card */
    void setIntCard(const INTCARD & IntCard);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(FUNCTIONAL & op);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(const Teuchos::RCP<FUNCTIONAL> & op);
    //@}
};


template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
traitsIntegratorFEL3d_OSC<FUNCTIONAL,linearTetra,TYPE,PRECISION>::
traitsIntegratorFEL3d_OSC()
{
  //Assert
  assert(staticAssert<TEST_GEOSHAPE::nDim    == 3>::returnValue);
  assert(staticAssert<TEST_FECARD::cardLabel == CR_OS_3d>::returnValue);
  
  //Build sub elements
  pointElement<GEOSHAPE> elem;
  
  for(UInt i=1; i <= 4; ++i)
  { elem.setPoint(i,linearTetra::getRefNodes(i)); }
  
  subElements = intOscSymplexSupport::refineTetra2(PRECISION -1,elem);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
traitsIntegratorFEL3d_OSC<FUNCTIONAL,linearTetra,TYPE,PRECISION>::
traitsIntegratorFEL3d_OSC(const INTCARD & IntCard)
{
  //Assert
  assert(staticAssert<TEST_GEOSHAPE::nDim    == 3>::returnValue);
  assert(staticAssert<TEST_FECARD::cardLabel == CR_OS_3d>::returnValue);
  
  //Build sub elements
  pointElement<GEOSHAPE> elem;
  
  for(UInt i=1; i <= 4; ++i)
  { elem.setPoint(i,linearTetra::getRefNodes(i)); }
  
  subElements = intOscSymplexSupport::refineTetra2(PRECISION -1,elem);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
void
traitsIntegratorFEL3d_OSC<FUNCTIONAL,linearTetra,TYPE,PRECISION>::
setIntCard(const INTCARD & IntCard)
{
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
sVect<Real>
traitsIntegratorFEL3d_OSC<FUNCTIONAL,linearTetra,TYPE,PRECISION>::
integration(FUNCTIONAL & op)
{
  //Alloc
  point3d P1, P2, P3, P4, YM, XM;
  komplex C;
    
  UInt matSize = op.numIndex_test();
  sVect<komplex> matC(matSize);  
  sVect<Real>     mat(matSize);
  
  //Card management
  TEST_FECARD testCard = op.getFeCardL_test();
  
  //Compute the global K and X
  sVect<point3d> Ktest(testCard.size());
  sVect<point3d> Xtest(testCard.size());
  sVect<point3d> nodes(op.getElCard().getNodes());
  
  tensor3d T;
  T.setCol(1,nodes(2) - nodes(1));
  T.setCol(2,nodes(3) - nodes(1));
  T.setCol(3,nodes(4) - nodes(1));
  T.computeInverse();
  
  for(UInt i=1; i <= testCard.size(); ++i)
  {
    Ktest(i) = T.firstIndexSaturation(testCard.getH(i));
    Xtest(i) = geoInterface.getPosition(nodes,testCard.getY(i));
  }

  //Compute the integral
  UInt tot;
  
  for(UInt k=1; k <= subElements.size(); ++k)
  {
    //Local nodes
    P1 = geoInterface.getPosition(op.getElCard().getNodes(), subElements(k).getPoint(1));
    P2 = geoInterface.getPosition(op.getElCard().getNodes(), subElements(k).getPoint(2));
    P3 = geoInterface.getPosition(op.getElCard().getNodes(), subElements(k).getPoint(3));
    P4 = geoInterface.getPosition(op.getElCard().getNodes(), subElements(k).getPoint(4));
    
    YM = (subElements(k).getPoint(1) 
        + subElements(k).getPoint(2)
        + subElements(k).getPoint(3)
        + subElements(k).getPoint(4)) / 4.0;

    XM = (P1 + P2 + P3 + P4) / 4.0;
    
    //Evaluate field
    op.eval(YM,matC);
    
    for(UInt i=1; i <= matSize; ++i)
    {
      C = intOscSymplexSupport::intS0_3d(Ktest(i), P1,P2,P3,P4)
        * komplex::iexp(-point3d::dot(Ktest(i), XM))
        * matC(i);

        mat(i) += komplex::cComp(C,op.I_test);
    }
  }

  return(mat);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
sVect<Real>
traitsIntegratorFEL3d_OSC<FUNCTIONAL,linearTetra,TYPE,PRECISION>::
integration(const Teuchos::RCP<FUNCTIONAL> & op)
{
  return(integration(*op));
}


#endif
