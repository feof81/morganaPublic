#ifndef TRAITSINTEGRATORFLA3D_OSC_HPP
#define TRAITSINTEGRATORFLA3D_OSC_HPP

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
class traitsIntegratorFLA3d_OSC
{ };


//_________________________________________________________________________________________________
// TETRA
//-------------------------------------------------------------------------------------------------
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
class traitsIntegratorFLA3d_OSC<FUNCTIONAL,linearTriangle,TYPE,PRECISION>
{
    /*! @name Typedefs */ //@{
  public:   
    typedef typename FUNCTIONAL::TEST_FECARD   TEST_FECARD;
    typedef typename FUNCTIONAL::TEST_GEOSHAPE TEST_GEOSHAPE;
    typedef typename FUNCTIONAL::TEST_PMAPTYPE TEST_PMAPTYPE;
    
    typedef typename FUNCTIONAL::INTGRID  INTGRID;
    typedef typename INTGRID::GRID_ELMAP  PMAPTYPE;
    
    typedef typename INTGRID::GEOSHAPE3D   GEOSHAPE3D;
    typedef typename GEOSHAPE3D::GEOBSHAPE GEOSHAPE2D;
    
    typedef elCard3d<GEOSHAPE3D,PMAPTYPE> ELCARD;
    typedef intPolicySTD                 INTCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const intClass intFlag = FLA3d_OSC;
    static const OPType   opFlag  = fnLA3d;
    
    sVect<pointElement<GEOSHAPE2D> > subElements;
    geoMapInterface<GEOSHAPE2D>   geoInterface2d;
    //@}
    
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    traitsIntegratorFLA3d_OSC();
    
    /*! Constructor */
    traitsIntegratorFLA3d_OSC(const INTCARD & IntCard);
    
    /*! Set integration card */
    void setIntCard(const INTCARD & IntCard);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(FUNCTIONAL & op);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(const Teuchos::RCP<FUNCTIONAL> & op);
    //@}
};


template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
traitsIntegratorFLA3d_OSC<FUNCTIONAL,linearTriangle,TYPE,PRECISION>::
traitsIntegratorFLA3d_OSC()
{
  //Assert
  assert(staticAssert<TEST_GEOSHAPE::nDim    == 2>::returnValue);
  assert(staticAssert<TEST_FECARD::cardLabel == CR_OS_3d>::returnValue);
  
  //Build sub elements
  pointElement<GEOSHAPE2D> elem;
  
  for(UInt i=1; i <= 3; ++i)
  { elem.setPoint(i,linearTriangle::getRefNodes(i)); }
  
  subElements = intOscSymplexSupport::refineTriangle(PRECISION -1,elem);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
traitsIntegratorFLA3d_OSC<FUNCTIONAL,linearTriangle,TYPE,PRECISION>::
traitsIntegratorFLA3d_OSC(const INTCARD & IntCard)
{
  //Assert
  assert(staticAssert<TEST_GEOSHAPE::nDim    == 2>::returnValue);
  assert(staticAssert<TEST_FECARD::cardLabel == CR_OS_3d>::returnValue);
  
  //Build sub elements
  pointElement<GEOSHAPE2D> elem;
  
  for(UInt i=1; i <= 3; ++i)
  { elem.setPoint(i,linearTriangle::getRefNodes(i)); }
  
  subElements = intOscSymplexSupport::refineTriangle(PRECISION -1, elem);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
void
traitsIntegratorFLA3d_OSC<FUNCTIONAL,linearTriangle,TYPE,PRECISION>::
setIntCard(const INTCARD & IntCard)
{
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
sVect<Real>
traitsIntegratorFLA3d_OSC<FUNCTIONAL,linearTriangle,TYPE,PRECISION>::
integration(FUNCTIONAL & fn)
{
  //Alloc
  point3d P1, P2, P3, YM, XM;
  komplex C;
  
  //Extracted data
  UInt el2d    = fn.getEl2d();
  UInt matSize = fn.numIndex_test();
  
  //Card management
  TEST_FECARD   testCard = fn.getFeCardL_test();
  sVect<point3d> nodes2d = fn.getGlobFacePoints();
  
  //Compute the global K and X
  sVect<point3d> Ktest(testCard.size());
  
  tensor3d T2d;
  T2d.setCol(1, nodes2d(2) - nodes2d(1));
  T2d.setCol(2, nodes2d(3) - nodes2d(1));
  T2d.completeThirdColoumn();
  T2d.computeInverse();
    
  for(UInt i=1; i <= testCard.size(); ++i)
  { Ktest(i) = T2d.firstIndexSaturation(testCard.getH(i)); }
  
  //Compute the integral
  sVect<komplex> matC(matSize);  
  sVect<Real>     mat(matSize);
  
  for(UInt k=1; k <= subElements.size(); ++k)
  {
    //Local nodes
    P1 = geoInterface2d.getPosition(fn.getGlobFacePoints(), subElements(k).getPoint(1));
    P2 = geoInterface2d.getPosition(fn.getGlobFacePoints(), subElements(k).getPoint(2));
    P3 = geoInterface2d.getPosition(fn.getGlobFacePoints(), subElements(k).getPoint(3));
    
    YM = (subElements(k).getPoint(1) 
        + subElements(k).getPoint(2)
        + subElements(k).getPoint(3)) / 3.0;

    XM = (P1 + P2 + P3) / 3.0;
    
    //Evaluate field
    fn.eval(YM,matC);
    
    for(UInt j=1; j <= fn.numIndex_test(); ++j)
    {
      C = intOscSymplexSupport::intS0_2d(Ktest(j), P1,P2,P3)
        * komplex::iexp(-point3d::dot(Ktest(j), XM))
        * matC(j);

      mat(j) += komplex::cComp(C,fn.I_test);
    }
  }

  return(mat);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
sVect<Real>
traitsIntegratorFLA3d_OSC<FUNCTIONAL,linearTriangle,TYPE,PRECISION>::
integration(const Teuchos::RCP<FUNCTIONAL> & fn)
{
  return(integration(*fn));
}

#endif
