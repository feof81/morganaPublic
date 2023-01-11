/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FUNCTIONALEL2D_HPP
#define FUNCTIONALEL2D_HPP

#include "elCardFeeder2d.hpp"
#include "morganaFiniteElements.hpp"


/*! Functional for standard integration 2d
<b> Scheme </b>
<ol>
<li> test 2d
<li> coeff 2d
<li> integration: elements 2d
</ol>
*/
template<typename TEST>
class functionalEL2d
{
    /*! @name Typedefs */ //@{
  public:    
    typedef typename TEST::FEINTERFACE        TEST_INTERFACE;
    typedef typename TEST::FECARD             TEST_FECARD;
    typedef typename TEST::FECARDS            TEST_FECARDS;
    typedef typename TEST::OPTIONS            TEST_OPTIONS;
    typedef typename TEST::GEOSHAPE           TEST_GEOSHAPE;
    typedef typename TEST::PMAPTYPE           TEST_PMAPTYPE;
    typedef typename TEST::OUTTYPE            TEST_OUTTYPE;
    typedef typename TEST::FIELD_DOFTYPE      TEST_DOFTYPE;
    typedef typename TEST_INTERFACE::LISTMAP  TEST_LISTMAP;
    
    typedef typename TEST::MESH2D                       MESH2D;
    typedef typename TEST::CONNECT2D                    CONNECT2D;
    typedef elCard2d<TEST_GEOSHAPE,TEST_PMAPTYPE>       ELCARD;
    typedef elCardFeeder2d<TEST_GEOSHAPE,TEST_PMAPTYPE> ELCARDFEEDER;
    
    typedef MESH2D INTGRID;
    typedef MESH2D LOOPGRID;
    //@}
    
    /*! @name Internal flags */ //@{
  public:
    bool startupOk;
    bool commDevLoaded;
    bool geometryLoaded;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    UInt el, I_test, J_test;
    Teuchos::RCP<MESH2D>     grid2d;
    TEST_INTERFACE  testInterface;
    ELCARDFEEDER    elCardFeeder;
    //@}
    
    /*! @name Constructors and setting - EXTERNAL */ //@{
  public:
    functionalEL2d();
    functionalEL2d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d);
    functionalEL2d(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d);
    virtual ~functionalEL2d() {};
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d);
    void setGeometry(MESH2D & Grid2d, CONNECT2D & ConnedGrid2d);
    //@}
    
    /*! @name Set map and startup - EXTERNAL */ //@{
  public:    
    void setFeCardL_test(const UInt & elL, const TEST_FECARD & FeCards);
    void setFeCardG_test(const UInt & elG, const TEST_FECARD & FeCards);
    void setFeCards_test(const TEST_FECARDS & FECards);
    void setOptions_test(const Teuchos::RCP<TEST_OPTIONS> & Options);
    void setOptions_test(TEST_OPTIONS & Options);
    
    void startup();
    void startup(const TEST & CloneTest);
    void startup(const Teuchos::RCP<TEST> & CloneTest);
    //@}
    
     /*! @name LocalVector - INTERFACE */ //@{
  public:
    void setElement(const UInt & El);
    void setIJ(const UInt & II_test, const UInt & JJ_test);
    const Teuchos::RCP<LOOPGRID> & getLoopGrid() const;
    //@}
    
    /*! @name Maps - INTERFACE */ //@{
  public:
    const UInt         & getNumDofsL_test() const;
    const UInt         & getNumDofsG_test() const;
    const UInt         & getSizeListL_test() const;
    const UInt         & getSizeListG_test() const;
    const TEST_LISTMAP & getListMap_test() const;
    //@}
    
    /*! @name Indices - INTERFACE */ //@{
  public:
    UInt numIndex_test();
    void indexL_test(sVect<UInt> & indices);
    void indexG_test(sVect<UInt> & indices);
    //@}
    
    /*! @name Integration - INTERFACE (elements 3d mesh) */ //@{
  public:
    const TEST_FECARD  & getFeCardL_test();
    ELCARD               getElCard() const;
    const Teuchos::RCP<MESH2D>  & getIntegrationGrid() const;
    const UInt         & getElement() const;
    //@}
    
    /*! @name Operator - INTERFACE */ //@{
  public:
    void eval_test(const point3d & Y, sVect<TEST_OUTTYPE> & val);
    void evalGrad_test(const point3d & Y, sVect<TEST_OUTTYPE> & gradX, sVect<TEST_OUTTYPE> & gradY, sVect<TEST_OUTTYPE> & gradZ);
    
    virtual void eval(const point3d & Y, sVect<komplex> & mat) = 0;
    virtual void eval(const point3d & Y, sVect<Real> & mat) = 0;
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS AND SETTING FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename TEST>
functionalEL2d<TEST>::
functionalEL2d()
{
  assert(staticAssert<TEST_GEOSHAPE::nDim == 2>::returnValue);
  
  startupOk      = false;
  commDevLoaded  = false;
  geometryLoaded = false;
  
  el = 0;
  
  I_test   = 0;
  J_test   = 0;
}

template<typename TEST>
functionalEL2d<TEST>::
functionalEL2d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d)
{
  assert(staticAssert<TEST_GEOSHAPE::nDim == 2>::returnValue);
  
  grid2d = Grid2d;
  
  testInterface.setCommunicator(CommDev);
  testInterface.setGeometry(Grid2d,ConnedGrid2d);
  
  elCardFeeder.setGeometry(Grid2d,ConnedGrid2d);
  
  startupOk      = false;
  commDevLoaded  = true;
  geometryLoaded = true;
  
  el = 0;
   
  I_test   = 0;
  J_test   = 0;
}

template<typename TEST>
functionalEL2d<TEST>::
functionalEL2d(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d)
{
  assert(staticAssert<TEST_GEOSHAPE::nDim == 2>::returnValue);
  
  grid2d = Teuchos::rcp(new MESH2D(Grid2d));
  
  testInterface.setGeometry(Grid2d,ConnedGrid2d);
  testInterface.setCommunicator(CommDev);
  
  elCardFeeder.setGeometry(Grid2d,ConnedGrid2d);
  
  startupOk      = false;
  commDevLoaded  = true;
  geometryLoaded = true;
  
  el = 0;
  
  I_test   = 0;
  J_test   = 0;
}

template<typename TEST>
void
functionalEL2d<TEST>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  testInterface.setCommunicator(CommDev);
  
  commDevLoaded = true;
}

template<typename TEST>
void
functionalEL2d<TEST>::
setCommunicator(communicator & CommDev)
{
  testInterface.setCommunicator(CommDev);
  
  commDevLoaded = true;
}

template<typename TEST>
void
functionalEL2d<TEST>::
setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d)
{
  grid2d = Grid2d;
  
  testInterface.setGeometry(Grid2d,ConnedGrid2d);
  elCardFeeder.setGeometry(Grid2d,ConnedGrid2d);
  
  geometryLoaded = true;
}

template<typename TEST>
void
functionalEL2d<TEST>::
setGeometry(MESH2D & Grid2d, CONNECT2D & ConnedGrid2d)
{
  grid2d = Teuchos::rcp(new MESH2D(Grid2d));
  
  testInterface.setGeometry(Grid2d,ConnedGrid2d);
  elCardFeeder.setGeometry(Grid2d,ConnedGrid2d);
  
  geometryLoaded = true;
}

template<typename TEST>
void
functionalEL2d<TEST>::
startup(const TEST & CloneTest)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  startupOk = true;
  
  testInterface.clone(CloneTest.getDofMapper());
}

template<typename TEST>
void
functionalEL2d<TEST>::
startup(const Teuchos::RCP<TEST> & CloneTest)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  startupOk = true;
  
  testInterface.clone(CloneTest->getDofMapper());
}



//_________________________________________________________________________________________________
// SET-MAP AND STARTUP
//-------------------------------------------------------------------------------------------------
template<typename TEST>
void
functionalEL2d<TEST>::
setFeCardL_test(const UInt & elL, const TEST_FECARD & FeCards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCardL(elL,FeCards);
}

template<typename TEST>
void
functionalEL2d<TEST>::
setFeCardG_test(const UInt & elG, const TEST_FECARD & FeCards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCardG(elG,FeCards);
}

template<typename TEST>
void
functionalEL2d<TEST>::
setFeCards_test(const TEST_FECARDS & FECards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCards(FECards);
}

template<typename TEST>
void
functionalEL2d<TEST>::
setOptions_test(const Teuchos::RCP<TEST_OPTIONS> & Options)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setOptions(Options);
}

template<typename TEST>
void
functionalEL2d<TEST>::
setOptions_test(TEST_OPTIONS & Options)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setOptions(Options);
}
    
template<typename TEST>
void
functionalEL2d<TEST>::
startup()
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  startupOk = true;
  
  testInterface.startup();
}


//_________________________________________________________________________________________________
// LOCAL VECTOR
//-------------------------------------------------------------------------------------------------
template<typename TEST>
void
functionalEL2d<TEST>::
setElement(const UInt & El)
{
  assert(startupOk);
  el = El;
}

template<typename TEST>
void
functionalEL2d<TEST>::
setIJ(const UInt & II_test, const UInt & JJ_test)
{
  assert(startupOk);
  
  I_test  = II_test;
  J_test  = JJ_test;
}

template<typename TEST>
const Teuchos::RCP<typename functionalEL2d<TEST>::LOOPGRID> &
functionalEL2d<TEST>::
getLoopGrid() const
{
  assert(startupOk);
  return(grid2d);
}



//_________________________________________________________________________________________________
// MAPS
//-------------------------------------------------------------------------------------------------
template<typename TEST>
const UInt &
functionalEL2d<TEST>::
getNumDofsL_test() const
{
  assert(startupOk);
  return(testInterface.getNumDofsL());
}

template<typename TEST>
const UInt &
functionalEL2d<TEST>::
getNumDofsG_test() const
{
  assert(startupOk);
  return(testInterface.getNumDofsG());
}

template<typename TEST>
const UInt &
functionalEL2d<TEST>::
getSizeListL_test() const
{
  assert(startupOk);
  return(testInterface.getSizeListL());
}

template<typename TEST>
const UInt &
functionalEL2d<TEST>::
getSizeListG_test() const
{
  assert(startupOk);
  return(testInterface.getSizeListG());
}

template<typename TEST>
const typename functionalEL2d<TEST>::TEST_LISTMAP &
functionalEL2d<TEST>::
getListMap_test() const
{
  assert(startupOk);
  return(testInterface.getListMap());
}



//_________________________________________________________________________________________________
// INDICES
//-------------------------------------------------------------------------------------------------
template<typename TEST>
UInt
functionalEL2d<TEST>::
numIndex_test()
{
  assert(startupOk);
  return(testInterface.feNumBasisL(el));
}

template<typename TEST>
void
functionalEL2d<TEST>::
indexL_test(sVect<UInt> & indices)
{
  assert(startupOk);
  return(testInterface.feIndicesLL(el,I_test,J_test,indices));
}

template<typename TEST>
void
functionalEL2d<TEST>::
indexG_test(sVect<UInt> & indices)
{
  assert(startupOk);
  return(testInterface.feIndicesLG(el,I_test,J_test,indices));
}



//_________________________________________________________________________________________________
// INTEGRATION
//-------------------------------------------------------------------------------------------------
template<typename TEST>
const typename functionalEL2d<TEST>::TEST_FECARD &
functionalEL2d<TEST>::
getFeCardL_test()
{
  assert(startupOk);
  return(testInterface.getFeCardL(el));
}

template<typename TEST>
typename functionalEL2d<TEST>::ELCARD
functionalEL2d<TEST>::
getElCard() const
{
  assert(startupOk);
  return(elCardFeeder.getCardLocal(el));
}

template<typename TEST>
const Teuchos::RCP<typename functionalEL2d<TEST>::MESH2D> &
functionalEL2d<TEST>::
getIntegrationGrid() const
{
  assert(startupOk);
  return(grid2d);
}

template<typename TEST>
const UInt &
functionalEL2d<TEST>::
getElement() const
{
  assert(startupOk);
  return(el);
}



//_________________________________________________________________________________________________
// OPERATOR
//-------------------------------------------------------------------------------------------------
template<typename TEST>
void
functionalEL2d<TEST>::
eval_test(const point3d & Y, sVect<TEST_OUTTYPE> & val)
{
  assert(startupOk);
  testInterface.feEvalL(el, I_test, J_test, Y, val);
}

template<typename TEST>
void
functionalEL2d<TEST>::
evalGrad_test(const point3d & Y, sVect<TEST_OUTTYPE> & gradX, sVect<TEST_OUTTYPE> & gradY, sVect<TEST_OUTTYPE> & gradZ)
{
  assert(startupOk);
  testInterface.feGradL(el, I_test, J_test, Y, gradX, gradY, gradZ);
}

#endif
