/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef TRAITSINTEGRATOROEL3D_OSC_HPP
#define TRAITSINTEGRATOROEL3D_OSC_HPP

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
template<typename OPERATOR, typename GEOSHAPE, intTypes TYPE, UInt PRECISION>
class traitsIntegratorOEL3d_OSC
{ };


//_________________________________________________________________________________________________
// TETRA
//-------------------------------------------------------------------------------------------------
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
class traitsIntegratorOEL3d_OSC<OPERATOR,linearTetra,TYPE,PRECISION>
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
    typedef intPolicySTD                INTCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const intClass intFlag = OEL3d_OSC;
    static const OPType   opFlag  = opEL3d;
    
    sVect<pointElement<GEOSHAPE> > subElements;
    geoMapInterface<GEOSHAPE>     geoInterface;
    //@}
    
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    traitsIntegratorOEL3d_OSC();
    
    /*! Constructor */
    traitsIntegratorOEL3d_OSC(const INTCARD & IntCard);
    
    /*! Set integration card */
    void setIntCard(const INTCARD & IntCard);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(OPERATOR & op);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(const Teuchos::RCP<OPERATOR> & op);
    //@}
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
traitsIntegratorOEL3d_OSC<OPERATOR,linearTetra,TYPE,PRECISION>::
traitsIntegratorOEL3d_OSC()
{
  //Assert
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 3>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName      == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  assert(staticAssert<FIELD_FECARD::cardLabel == CR_OS_3d>::returnValue);
  assert(staticAssert<TEST_FECARD::cardLabel  == CR_OS_3d>::returnValue);
  
  //Build sub elements
  pointElement<GEOSHAPE> elem;
  
  for(UInt i=1; i <= 4; ++i)
  { elem.setPoint(i,linearTetra::getRefNodes(i)); }
  
  subElements = intOscSymplexSupport::refineTetra2(PRECISION -1,elem);
}
    
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
traitsIntegratorOEL3d_OSC<OPERATOR,linearTetra,TYPE,PRECISION>::
traitsIntegratorOEL3d_OSC(const INTCARD & IntCard)
{
  //Assert
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 3>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName      == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  assert(staticAssert<FIELD_FECARD::cardLabel == CR_OS_3d>::returnValue);
  assert(staticAssert<TEST_FECARD::cardLabel  == CR_OS_3d>::returnValue);
  
  //Build sub elements
  pointElement<GEOSHAPE> elem;
  
  for(UInt i=1; i <= 4; ++i)
  { elem.setPoint(i,linearTetra::getRefNodes(i)); }
  
  subElements = intOscSymplexSupport::refineTetra2(PRECISION -1, elem);
}
    
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
void
traitsIntegratorOEL3d_OSC<OPERATOR,linearTetra,TYPE,PRECISION>::
setIntCard(const INTCARD & IntCard)
{
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
sVect<Real>
traitsIntegratorOEL3d_OSC<OPERATOR,linearTetra,TYPE,PRECISION>::
integration(OPERATOR & op)
{
  //Alloc
  point3d P1, P2, P3, P4, YM, XM;
  komplex C;
    
  UInt matSize = op.numIndex_field() * op.numIndex_test();
  sVect<komplex> matC(matSize);  
  sVect<Real>     mat(matSize);
  
  //Card management
  FIELD_FECARD fieldCard = op.getFeCardL_field();
  TEST_FECARD   testCard = op.getFeCardL_test();
  
  //Compute the global K and X
  sVect<point3d> Kfield(fieldCard.size()), Ktest(testCard.size());
  sVect<point3d> Xfield(fieldCard.size()), Xtest(testCard.size());
  sVect<point3d> nodes(op.getElCard().getNodes());
  
  tensor3d T;
  T.setCol(1,nodes(2) - nodes(1));
  T.setCol(2,nodes(3) - nodes(1));
  T.setCol(3,nodes(4) - nodes(1));
  T.computeInverse();
  
  for(UInt i=1; i <= fieldCard.size(); ++i)
  {
    Kfield(i) = T.firstIndexSaturation(fieldCard.getH(i));
    Xfield(i) = geoInterface.getPosition(nodes,fieldCard.getY(i));
  }
  
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
    
    tot = 1;

    for(UInt i=1; i <= op.numIndex_field(); ++i)
    {
      for(UInt j=1; j <= op.numIndex_test(); ++j)
      {
        C = intOscSymplexSupport::intS0_3d(Kfield(i) + Ktest(j), P1,P2,P3,P4)
          * komplex::iexp(-point3d::dot(Kfield(i) + Ktest(j), XM))
          * matC(tot);
	  
	  
	 if((tot == 7) && (k == 7))
	 {
	   cout << "K   : " << std::setprecision(16) << Kfield(i) + Ktest(j);
	   cout << "P1  : " << std::setprecision(16) << P1;
	   cout << "P2  : " << std::setprecision(16) << P2;
	   cout << "P3  : " << std::setprecision(16) << P3;
	   cout << "P4  : " << std::setprecision(16) << P4;
	   cout << "Ex  : " << std::setprecision(16) << intOscSymplexSupport::intS0_3d(Kfield(i) + Ktest(j), P1,P2,P3,P4);
	   cout << C;
	 }
	  
	  
	  /*cout << "i j :" << i << " " << j << endl;
	  
	  cout << "A   : " << intOscSymplexSupport::intS0_3d(Kfield(i) + Ktest(j), P1,P2,P3,P4);
	  cout << "B   : " << komplex::iexp(-point3d::dot(Kfield(i) + Ktest(j), XM));
	  cout << "C   : " << matC(tot);
	  cout << "Tot : " << komplex::cComp(C,op.I_test) << endl;*/

        mat(tot) += komplex::cComp(C,op.I_test);

        ++tot;
      }
    } 
  }

  return(mat);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
sVect<Real>
traitsIntegratorOEL3d_OSC<OPERATOR,linearTetra,TYPE,PRECISION>::
integration(const Teuchos::RCP<OPERATOR> & op)
{
  return(integration(*op));
}


#endif
