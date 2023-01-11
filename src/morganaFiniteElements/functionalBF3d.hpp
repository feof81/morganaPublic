/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FUNCTIONALBF3D_HPP
#define FUNCTIONALBF3D_HPP


#include "elCardFeeder3d.hpp"
#include "supportLagrange3d.hpp"
#include "morganaFiniteElements.hpp"


/*! Functional for boundary integration 3d
<b> Scheme </b>
<ol>
<li> test 3d
<li> coeff 2d
<li> integration: boundary faces 3d
</ol>
*/
template<typename TEST>
class functionalBF3d
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
    
    typedef TEST_PMAPTYPE PMAPTYPE;
    
    typedef typename TEST::MESH3D                    MESH3D;
    typedef typename TEST::CONNECT3D                 CONNECT3D;
    typedef typename MESH3D::GEOSHAPE2D              GEOSHAPE2D;
    typedef mesh2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>     MESH2D;
    typedef connect2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>  CONNECT2D;
    typedef typename MESH3D::GEOSHAPE3D              GEOSHAPE3D;
    
    typedef elCard3d<TEST_GEOSHAPE,PMAPTYPE>       ELCARD;
    typedef elCardFeeder3d<TEST_GEOSHAPE,PMAPTYPE> ELCARDFEEDER;
    typedef supportLagrange3d<GEOSHAPE3D,PMAPTYPE> SUPPORTLAGRANGE; 
    
    typedef MESH3D INTGRID;
    typedef MESH3D LOOPGRID;
    //@}
    
    /*! @name Internal flags */ //@{
  public:
    bool startupOk;
    bool commDevLoaded;
    bool geometryLoaded;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    UInt el3d, locFc, numBoundary, I_test, J_test;
    Teuchos::RCP<MESH2D>     grid2d;
    Teuchos::RCP<MESH3D>     grid3d;
    Teuchos::RCP<CONNECT3D>  connectGrid3d;
    Teuchos::RCP<CONNECT2D>  connectGrid2d;
    
    TEST_INTERFACE  testInterface;
    ELCARDFEEDER    elCardFeeder;
    SUPPORTLAGRANGE supportLagrange;
    //@}
    
    /*! @name Constructors and setting - EXTERNAL */ //@{
  public:
    functionalBF3d();
    functionalBF3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    functionalBF3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnectGrid3d, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d);
    virtual ~functionalBF3d() {};
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    void setGeometry(MESH3D & Grid3d, CONNECT3D & ConnectGrid3d, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d);
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
    
    void setLocFace(const UInt & LocFc);
    bool isBoundary() const;
    UInt getNumBoundary() const;
    UInt getEl2d() const;
    UInt getEl3d() const;
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
    
    /*! @name Integration - INTERFACE (boundary faces of the 3d mesh) */ //@{
  public:
    const TEST_FECARD          & getFeCardL_test();
    ELCARD                       getElCard() const;
    const Teuchos::RCP<MESH3D> & getIntegrationGrid() const;
    const UInt                 & getElement() const;
    const sVect<point3d>       & getGlobFacePoints() const;
    //@}
    
    /*! @name Operator - INTERFACE */ //@{
  public:
    point3d mapVolumeY(const point3d & Yf) const;
    point3d mapSurfaceY(const point3d & Yf);
    point3d computeNormal(const point3d & Yf) const;
    
    void eval_test(const point3d & Yf, sVect<TEST_OUTTYPE> & val);
    void evalGrad_test(const point3d & Yf, sVect<TEST_OUTTYPE> & gradX, sVect<TEST_OUTTYPE> & gradY, sVect<TEST_OUTTYPE> & gradZ);
    
    virtual void eval(const point3d & Y, sVect<komplex> & mat) = 0;
    virtual void eval(const point3d & Y, sVect<Real> & mat) = 0;
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS AND SETTING FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename TEST>
functionalBF3d<TEST>::
functionalBF3d()
{
  assert(staticAssert<TEST_GEOSHAPE::nDim == 3>::returnValue);
  
  startupOk      = false;
  commDevLoaded  = false;
  geometryLoaded = false;
  
  el3d  = 0;
  locFc = 1;
  
  I_test  = 0;
  J_test  = 0; 
}

template<typename TEST>
functionalBF3d<TEST>::
functionalBF3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  assert(staticAssert<TEST_GEOSHAPE::nDim == 3>::returnValue);
  
  grid2d        = Grid2d;
  grid3d        = Grid3d;
  connectGrid2d = ConnectGrid2d;
  connectGrid3d = ConnectGrid3d;
  
  testInterface.setCommunicator(CommDev);
  testInterface.setGeometry(Grid3d,ConnectGrid3d);
  
  elCardFeeder.setGeometry(Grid3d,ConnectGrid3d);
  supportLagrange.setGeometry(Grid3d,ConnectGrid3d,Grid2d,ConnectGrid2d);
  
  startupOk      = false;
  commDevLoaded  = true;
  geometryLoaded = true;
  
  el3d  = 0;
  locFc = 1;
  
  I_test  = 0;
  J_test  = 0; 
}

template<typename TEST>
functionalBF3d<TEST>::
functionalBF3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnectGrid3d, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d)
{
  assert(staticAssert<TEST_GEOSHAPE::nDim == 3>::returnValue);
  
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  grid3d        = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
  connectGrid3d = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
  
  testInterface.setCommunicator(CommDev);
  testInterface.setGeometry(Grid3d,ConnectGrid3d);
  
  elCardFeeder.setGeometry(Grid3d,ConnectGrid3d);
  supportLagrange.setGeometry(Grid3d,ConnectGrid3d,Grid2d,ConnectGrid2d);
  
  startupOk      = false;
  commDevLoaded  = true;
  geometryLoaded = true;
  
  el3d  = 0;
  locFc = 1;
  
  I_test  = 0;
  J_test  = 0; 
}

template<typename TEST>
void
functionalBF3d<TEST>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  testInterface.setCommunicator(CommDev);
  
  commDevLoaded = true;
}

template<typename TEST>
void
functionalBF3d<TEST>::
setCommunicator(communicator & CommDev)
{
  testInterface.setCommunicator(CommDev);
  
  commDevLoaded = true;
}

template<typename TEST>
void
functionalBF3d<TEST>::
setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  grid2d        = Grid2d;
  grid3d        = Grid3d;
  connectGrid2d = ConnectGrid2d;
  connectGrid3d = ConnectGrid3d;
  
  testInterface.setGeometry(Grid3d,ConnectGrid3d);
  
  elCardFeeder.setGeometry(Grid3d,ConnectGrid3d);
  supportLagrange.setGeometry(Grid3d,ConnectGrid3d,Grid2d,ConnectGrid2d);
  
  geometryLoaded = true;
}

template<typename TEST>
void
functionalBF3d<TEST>::
setGeometry(MESH3D & Grid3d, CONNECT3D & ConnectGrid3d, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d)
{
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  grid3d        = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
  connectGrid3d = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
  
  testInterface.setGeometry(Grid3d,ConnectGrid3d);
  
  elCardFeeder.setGeometry(Grid3d,ConnectGrid3d);
  supportLagrange.setGeometry(Grid3d,ConnectGrid3d,Grid2d,ConnectGrid2d);
  
  geometryLoaded = true;
}



//_________________________________________________________________________________________________
// SET-MAP AND STARTUP
//-------------------------------------------------------------------------------------------------
template<typename TEST>
void
functionalBF3d<TEST>::
setFeCardL_test(const UInt & elL, const TEST_FECARD & FeCards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCardL(elL,FeCards);
}

template<typename TEST>
void
functionalBF3d<TEST>::
setFeCardG_test(const UInt & elG, const TEST_FECARD & FeCards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCardG(elG,FeCards);
}

template<typename TEST>
void
functionalBF3d<TEST>::
setFeCards_test(const TEST_FECARDS & FECards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCards(FECards);
}

template<typename TEST>
void
functionalBF3d<TEST>::
setOptions_test(const Teuchos::RCP<TEST_OPTIONS> & Options)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setOptions(Options);
}

template<typename TEST>
void
functionalBF3d<TEST>::
setOptions_test(TEST_OPTIONS & Options)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setOptions(Options);
}
    
template<typename TEST>
void
functionalBF3d<TEST>::
startup()
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  startupOk = true;
  
  testInterface.startup();
}

template<typename TEST>
void
functionalBF3d<TEST>::
startup(const TEST & CloneTest)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  startupOk = true;
  
  testInterface.clone(CloneTest.getDofMapper());
}

template<typename TEST>
void
functionalBF3d<TEST>::
startup(const Teuchos::RCP<TEST> & CloneTest)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  startupOk = true;
  
  testInterface.clone(CloneTest->getDofMapper());
}



//_________________________________________________________________________________________________
// LOCAL MATRIX
//-------------------------------------------------------------------------------------------------
template<typename TEST>
void
functionalBF3d<TEST>::
setElement(const UInt & El)
{
  assert(startupOk);
  
  el3d = El;
  supportLagrange.setElement3d(El);
  
  //Count boundary faces
  numBoundary = 0;
  
  for(UInt f=1; f <= GEOSHAPE3D::numFaces; ++f)
  {
    supportLagrange.setLocalFace(f);
    numBoundary += UInt(supportLagrange.isBoundary());
  }
  
  supportLagrange.setLocalFace(locFc);
}

template<typename TEST>
void
functionalBF3d<TEST>::
setIJ(const UInt & II_test, const UInt & JJ_test)
{
  assert(startupOk);
  I_test  = II_test;
  J_test  = JJ_test;
}

template<typename TEST>
const Teuchos::RCP<typename functionalBF3d<TEST>::LOOPGRID> &
functionalBF3d<TEST>::
getLoopGrid() const
{
  assert(startupOk);
  return(grid3d);
}
    
template<typename TEST>
void
functionalBF3d<TEST>::
setLocFace(const UInt & LocFc)
{ 
  assert(startupOk);
  
  locFc = LocFc;
  supportLagrange.setLocalFace(LocFc);
}

template<typename TEST>
bool
functionalBF3d<TEST>::
isBoundary() const
{
  assert(startupOk);
  return(supportLagrange.isBoundary());
}

template<typename TEST>
UInt
functionalBF3d<TEST>::
getNumBoundary() const
{
  assert(startupOk);
  return(numBoundary);
}

template<typename TEST>
UInt
functionalBF3d<TEST>::
getEl2d() const
{
  return(supportLagrange.getElement2d());
}

template<typename TEST>
UInt
functionalBF3d<TEST>::
getEl3d() const
{
  return(el3d);
}




//_________________________________________________________________________________________________
// MAPS
//-------------------------------------------------------------------------------------------------
template<typename TEST>
const UInt &
functionalBF3d<TEST>::
getNumDofsL_test() const
{
  assert(startupOk);
  return(testInterface.getNumDofsL());
}

template<typename TEST>
const UInt &
functionalBF3d<TEST>::
getNumDofsG_test() const
{
  assert(startupOk);
  return(testInterface.getNumDofsG());
}

template<typename TEST>
const UInt &
functionalBF3d<TEST>::
getSizeListL_test() const
{
  assert(startupOk);
  return(testInterface.getSizeListL());
}

template<typename TEST>
const UInt &
functionalBF3d<TEST>::
getSizeListG_test() const
{
  assert(startupOk);
  return(testInterface.getSizeListG());
}

template<typename TEST>
const typename functionalBF3d<TEST>::TEST_LISTMAP &
functionalBF3d<TEST>::
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
functionalBF3d<TEST>::
numIndex_test()
{
  assert(startupOk);
  return( testInterface.feNumBasisL(el3d) );
}

template<typename TEST>
void
functionalBF3d<TEST>::
indexL_test(sVect<UInt> & indices)
{
  assert(startupOk);
  return(testInterface.feIndicesLL(el3d,I_test,J_test,indices));
}

template<typename TEST>
void
functionalBF3d<TEST>::
indexG_test(sVect<UInt> & indices)
{
  assert(startupOk);
  return(testInterface.feIndicesLG(el3d,I_test,J_test,indices));
}



//_________________________________________________________________________________________________
// INTEGRATION
//-------------------------------------------------------------------------------------------------
template<typename TEST>
const typename functionalBF3d<TEST>::TEST_FECARD &
functionalBF3d<TEST>::
getFeCardL_test()
{
  assert(startupOk);
  return(testInterface.getFeCardL(el3d));
}

template<typename TEST>
typename functionalBF3d<TEST>::ELCARD
functionalBF3d<TEST>::
getElCard() const
{
  assert(startupOk);
  return(elCardFeeder.getCardLocal(el3d));
}

template<typename TEST>
const Teuchos::RCP<typename functionalBF3d<TEST>::MESH3D> &
functionalBF3d<TEST>::
getIntegrationGrid() const
{
  assert(startupOk);
  return(grid3d);
}

template<typename TEST>
const UInt &
functionalBF3d<TEST>::
getElement() const
{
  assert(startupOk);
  return(el3d);
}

template<typename TEST>
const sVect<point3d> &
functionalBF3d<TEST>::
getGlobFacePoints() const
{
  assert(startupOk);
  return(supportLagrange.getGlobFacePoints());
}



//_________________________________________________________________________________________________
// OPERATOR
//-------------------------------------------------------------------------------------------------
template<typename TEST>
point3d
functionalBF3d<TEST>::
mapVolumeY(const point3d & Yf) const
{
  assert(startupOk);
  return(supportLagrange.mapVolumeY(Yf));
}

template<typename TEST>
point3d
functionalBF3d<TEST>::
computeNormal(const point3d & Yf) const
{
  assert(startupOk);
  return(supportLagrange.computeNormal(Yf));
}

template<typename TEST>
point3d
functionalBF3d<TEST>::
mapSurfaceY(const point3d & Yf)
{
  assert(startupOk);
  return(supportLagrange.mapSurfaceY(Yf));
}

template<typename TEST>
void
functionalBF3d<TEST>::
eval_test(const point3d & Yf, sVect<TEST_OUTTYPE> & val)
{
  assert(startupOk);
  testInterface.feEvalL(el3d, I_test, J_test, supportLagrange.mapVolumeY(Yf), val);
}

template<typename TEST>
void
functionalBF3d<TEST>::
evalGrad_test(const point3d & Yf, sVect<TEST_OUTTYPE> & gradX, sVect<TEST_OUTTYPE> & gradY, sVect<TEST_OUTTYPE> & gradZ)
{
  assert(startupOk);
  testInterface.feGradL(el3d, I_test, J_test, supportLagrange.mapVolumeY(Yf), gradX, gradY, gradZ);
}


#endif
