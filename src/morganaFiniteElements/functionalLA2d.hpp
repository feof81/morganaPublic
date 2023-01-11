/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FUNCTIONALLA2D_HPP
#define FUNCTIONALLA2D_HPP

#include "elCardFeeder2d.hpp"
#include "supportLagrange2d.hpp"
#include "morganaFiniteElements.hpp"

/*! Functional Lagrange BC imposition 2d
<b> Scheme </b>
<ol>
<li> test 1d
<li> coeff 1d
<li> integration: boundary edges 2d
</ol>
*/
template<typename TEST, typename GEOSHAPE2D>
class functionalLA2d
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
    
    typedef typename TEST::MESH1D                    MESH1D;
    typedef typename TEST::CONNECT1D                 CONNECT1D;
    typedef mesh2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>     MESH2D;
    typedef connect2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>  CONNECT2D;
    
    typedef elCard2d<TEST_GEOSHAPE,PMAPTYPE>       ELCARD;
    typedef elCardFeeder2d<GEOSHAPE2D,PMAPTYPE>    ELCARDFEEDER;
    typedef supportLagrange2d<GEOSHAPE2D,PMAPTYPE> SUPPORTLAGRANGE; 
    
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
    UInt el2d, locEd, numBoundary, I_test, J_test;
    Teuchos::RCP<MESH1D>     grid1d;
    Teuchos::RCP<MESH2D>     grid2d;
    Teuchos::RCP<CONNECT2D>  connectGrid2d;
    Teuchos::RCP<CONNECT1D>  connectGrid1d;
    
    TEST_INTERFACE  testInterface;
    ELCARDFEEDER    elCardFeeder;
    SUPPORTLAGRANGE supportLagrange;
    //@}
    
    /*! @name Constructors and setting - EXTERNAL */ //@{
  public:
    functionalLA2d();
    functionalLA2d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnectGrid1d);
    functionalLA2d(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d, MESH1D & Grid1d, CONNECT1D & ConnectGrid1d);
    virtual ~functionalLA2d() {};
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnectGrid1d);
    void setGeometry(MESH2D & Grid2d, CONNECT2D & ConnectGrid2d, MESH1D & Grid1d, CONNECT1D & ConnectGrid1d);
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
    
    void setLocEdge(const UInt & LocEd);
    bool isBoundary() const;
    UInt getNumBoundary() const;
    UInt getEl1d() const;
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
    const TEST_FECARD           & getFeCardL_test();
    ELCARD                        getElCard() const;
    const Teuchos::RCP<MESH2D>  & getIntegrationGrid() const;
    const UInt                  & getElement() const;
    const sVect<point3d>        & getGlobEdgePoints() const;
    //@}
    
    /*! @name Operator - INTERFACE */ //@{
  public:
    point3d mapVolumeY(const point3d & Yd) const;
    point3d mapSurfaceY(const point3d & Yd) const;
    point3d computeNormal(const point3d & Yd) const;
    
    void eval_test(const point3d & Yd, sVect<TEST_OUTTYPE> & val);
    void evalGrad_test(const point3d & Yd, sVect<TEST_OUTTYPE> & gradX, sVect<TEST_OUTTYPE> & gradY, sVect<TEST_OUTTYPE> & gradZ);
    
    virtual void eval(const point3d & Y, sVect<komplex> & mat) = 0;
    virtual void eval(const point3d & Y, sVect<Real> & mat) = 0;
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS AND SETTING FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename TEST, typename GEOSHAPE2D>
functionalLA2d<TEST,GEOSHAPE2D>::
functionalLA2d()
{
  typedef typename GEOSHAPE2D::GEOBSHAPE GEOTESTING;
  
  assert(staticAssert<TEST_GEOSHAPE::nDim == 1>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::geoName == GEOTESTING::geoName>::returnValue);
  
  startupOk      = false;
  commDevLoaded  = false;
  geometryLoaded = false;
  
  el2d  = 0;
  locEd = 1;
  
  I_test  = 0;
  J_test  = 0; 
}

template<typename TEST, typename GEOSHAPE2D>
functionalLA2d<TEST,GEOSHAPE2D>::
functionalLA2d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnectGrid1d)
{
  typedef typename GEOSHAPE2D::GEOBSHAPE GEOTESTING;
  
  assert(staticAssert<TEST_GEOSHAPE::nDim == 1>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::geoName == GEOTESTING::geoName>::returnValue);
  
  grid1d        = Grid1d;
  grid2d        = Grid2d;
  connectGrid1d = ConnectGrid1d;
  connectGrid2d = ConnectGrid2d;
  
  testInterface.setCommunicator(CommDev);
  testInterface.setGeometry(Grid1d,ConnectGrid1d);
  
  elCardFeeder.setGeometry(Grid2d,ConnectGrid2d);
  supportLagrange.setGeometry(Grid2d,ConnectGrid2d,Grid1d,ConnectGrid1d);
  
  startupOk      = false;
  commDevLoaded  = true;
  geometryLoaded = true;
  
  el2d  = 0;
  locEd = 1;
  
  I_test  = 0;
  J_test  = 0; 
}

template<typename TEST, typename GEOSHAPE2D>
functionalLA2d<TEST,GEOSHAPE2D>::
functionalLA2d(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d, MESH1D & Grid1d, CONNECT1D & ConnectGrid1d)
{
  typedef typename GEOSHAPE2D::GEOBSHAPE GEOTESTING;
  
  assert(staticAssert<TEST_GEOSHAPE::nDim == 1>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::geoName == GEOTESTING::geoName>::returnValue);
  
  grid1d        = Teuchos::rcp(new MESH1D(Grid1d));
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  connectGrid1d = Teuchos::rcp(new CONNECT1D(ConnectGrid1d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
  
  testInterface.setCommunicator(CommDev);
  testInterface.setGeometry(Grid1d,ConnectGrid1d);
  
  elCardFeeder.setGeometry(Grid2d,ConnectGrid2d);
  supportLagrange.setGeometry(Grid2d,ConnectGrid2d,Grid1d,ConnectGrid1d);
  
  startupOk      = false;
  commDevLoaded  = true;
  geometryLoaded = true;
  
  el2d  = 0;
  locEd = 1;
  
  I_test  = 0;
  J_test  = 0; 
}

template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  testInterface.setCommunicator(CommDev);
  
  commDevLoaded = true;
}

template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
setCommunicator(communicator & CommDev)
{
  testInterface.setCommunicator(CommDev);
  
  commDevLoaded = true;
}

template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnectGrid1d)
{
  grid1d        = Grid1d;
  grid2d        = Grid2d;
  connectGrid1d = ConnectGrid1d;
  connectGrid2d = ConnectGrid2d;
  
  testInterface.setGeometry(Grid1d,ConnectGrid1d);
  
  elCardFeeder.setGeometry(Grid2d,ConnectGrid2d);
  supportLagrange.setGeometry(Grid2d,ConnectGrid2d,Grid1d,ConnectGrid1d);
  
  geometryLoaded = true;
}

template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
setGeometry(MESH2D & Grid2d, CONNECT2D & ConnectGrid2d, MESH1D & Grid1d, CONNECT1D & ConnectGrid1d)
{
  grid1d        = Teuchos::rcp(new MESH1D(Grid1d));
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  connectGrid1d = Teuchos::rcp(new CONNECT1D(ConnectGrid1d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
  
  testInterface.setGeometry(Grid1d,ConnectGrid1d);
  
  elCardFeeder.setGeometry(Grid2d,ConnectGrid2d);
  supportLagrange.setGeometry(Grid2d,ConnectGrid2d,Grid1d,ConnectGrid1d);
  
  geometryLoaded = true;
}



//_________________________________________________________________________________________________
// SET-MAP AND STARTUP
//-------------------------------------------------------------------------------------------------
template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
setFeCardL_test(const UInt & elL, const TEST_FECARD & FeCards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCardL(elL,FeCards);
}

template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
setFeCardG_test(const UInt & elG, const TEST_FECARD & FeCards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCardG(elG,FeCards);
}

template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
setFeCards_test(const TEST_FECARDS & FECards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCards(FECards);
}

template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
setOptions_test(const Teuchos::RCP<TEST_OPTIONS> & Options)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setOptions(Options);
}

template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
setOptions_test(TEST_OPTIONS & Options)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setOptions(Options);
}
    
template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
startup()
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  startupOk = true;
  
  testInterface.startup();
}

template<typename TEST, typename GEOSHAPE3D>
void
functionalLA2d<TEST,GEOSHAPE3D>::
startup(const TEST & CloneTest)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  startupOk = true;
  
  testInterface.clone(CloneTest.getDofMapper());
}

template<typename TEST, typename GEOSHAPE3D>
void
functionalLA2d<TEST,GEOSHAPE3D>::
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
template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
setElement(const UInt & El)
{
  assert(startupOk);
  
  el2d = El;
  supportLagrange.setElement2d(El);
  
  //Count boundary faces
  numBoundary = 0;
  
  for(UInt e=1; e <= GEOSHAPE2D::numEdges; ++e)
  {
    supportLagrange.setLocalEdge(e);
    numBoundary += UInt(supportLagrange.isBoundary());
  }
  
  supportLagrange.setLocalEdge(locEd);
}

template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
setIJ(const UInt & II_test, const UInt & JJ_test)
{
  assert(startupOk);
  I_test  = II_test;
  J_test  = JJ_test;
}

template<typename TEST, typename GEOSHAPE2D>
const Teuchos::RCP<typename functionalLA2d<TEST,GEOSHAPE2D>::LOOPGRID> &
functionalLA2d<TEST,GEOSHAPE2D>::
getLoopGrid() const
{
  assert(startupOk);
  return(grid2d);
}

template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
setLocEdge(const UInt & LocEd)
{
  assert(startupOk);
  
  locEd = LocEd;
  supportLagrange.setLocalEdge(LocEd);
}

template<typename TEST, typename GEOSHAPE2D>
bool
functionalLA2d<TEST,GEOSHAPE2D>::
isBoundary() const
{
  assert(startupOk);
  return(supportLagrange.isBoundary());
}

template<typename TEST, typename GEOSHAPE2D>
UInt
functionalLA2d<TEST,GEOSHAPE2D>::
getNumBoundary() const
{
  assert(startupOk);
  return(numBoundary);
}

template<typename TEST, typename GEOSHAPE2D>
UInt
functionalLA2d<TEST,GEOSHAPE2D>::
getEl1d() const
{
  return(supportLagrange.getElement1d());
}



//_________________________________________________________________________________________________
// MAPS
//-------------------------------------------------------------------------------------------------
template<typename TEST, typename GEOSHAPE2D>
const UInt &
functionalLA2d<TEST,GEOSHAPE2D>::
getNumDofsL_test() const
{
  assert(startupOk);
  return(testInterface.getNumDofsL());
}

template<typename TEST, typename GEOSHAPE2D>
const UInt &
functionalLA2d<TEST,GEOSHAPE2D>::
getNumDofsG_test() const
{
  assert(startupOk);
  return(testInterface.getNumDofsG());
}

template<typename TEST, typename GEOSHAPE2D>
const UInt &
functionalLA2d<TEST,GEOSHAPE2D>::
getSizeListL_test() const
{
  assert(startupOk);
  return(testInterface.getSizeListL());
}

template<typename TEST, typename GEOSHAPE2D>
const UInt &
functionalLA2d<TEST,GEOSHAPE2D>::
getSizeListG_test() const
{
  assert(startupOk);
  return(testInterface.getSizeListG());
}

template<typename TEST, typename GEOSHAPE2D>
const typename functionalLA2d<TEST,GEOSHAPE2D>::TEST_LISTMAP &
functionalLA2d<TEST,GEOSHAPE2D>::
getListMap_test() const
{
  assert(startupOk);
  return(testInterface.getListMap());
}



//_________________________________________________________________________________________________
// INDICES
//-------------------------------------------------------------------------------------------------
template<typename TEST, typename GEOSHAPE2D>
UInt
functionalLA2d<TEST,GEOSHAPE2D>::
numIndex_test()
{
  assert(startupOk);
  return( testInterface.feNumBasisL(supportLagrange.getElement1d()) );
}

template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
indexL_test(sVect<UInt> & indices)
{
  assert(startupOk);
  return(testInterface.feIndicesLL(supportLagrange.getElement1d(),I_test,J_test,indices));
}

template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
indexG_test(sVect<UInt> & indices)
{
  assert(startupOk);
  return(testInterface.feIndicesLG(supportLagrange.getElement1d(),I_test,J_test,indices));
}



//_________________________________________________________________________________________________
// INTEGRATION
//-------------------------------------------------------------------------------------------------
template<typename TEST, typename GEOSHAPE2D>
const typename functionalLA2d<TEST,GEOSHAPE2D>::TEST_FECARD &
functionalLA2d<TEST,GEOSHAPE2D>::
getFeCardL_test()
{
  assert(startupOk);
  return(testInterface.getFeCardL(supportLagrange.getElement1d()));
}

template<typename TEST, typename GEOSHAPE2D>
typename functionalLA2d<TEST,GEOSHAPE2D>::ELCARD
functionalLA2d<TEST,GEOSHAPE2D>::
getElCard() const
{
  assert(startupOk);
  return(elCardFeeder.getCardLocal(supportLagrange.getElement1d()));
}

template<typename TEST, typename GEOSHAPE2D>
const Teuchos::RCP<typename functionalLA2d<TEST,GEOSHAPE2D>::LOOPGRID> &
functionalLA2d<TEST,GEOSHAPE2D>::
getIntegrationGrid() const
{
  assert(startupOk);
  return(grid2d);
}

template<typename TEST, typename GEOSHAPE2D>
const UInt &
functionalLA2d<TEST,GEOSHAPE2D>::
getElement() const
{
  assert(startupOk);
  return(el2d);
}

template<typename TEST, typename GEOSHAPE2D>
const sVect<point3d> &
functionalLA2d<TEST,GEOSHAPE2D>::
getGlobEdgePoints() const
{
  assert(startupOk);
  return(supportLagrange.getGlobEdgePoints());
}



//_________________________________________________________________________________________________
// OPERATOR
//-------------------------------------------------------------------------------------------------
template<typename TEST, typename GEOSHAPE2D>
point3d
functionalLA2d<TEST,GEOSHAPE2D>::
mapVolumeY(const point3d & Yd) const
{
  assert(startupOk);
  return(supportLagrange.mapVolumeY(Yd));
}

template<typename TEST, typename GEOSHAPE3D>
point3d
functionalLA2d<TEST,GEOSHAPE3D>::
mapSurfaceY(const point3d & Yd) const
{
  assert(startupOk);
  return(supportLagrange.mapSurfaceY(Yd));
}

template<typename TEST, typename GEOSHAPE2D>
point3d
functionalLA2d<TEST,GEOSHAPE2D>::
computeNormal(const point3d & Yd) const
{
  assert(startupOk);
  return(supportLagrange.computeNormal(Yd));
}

template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
eval_test(const point3d & Yd, sVect<TEST_OUTTYPE> & val)
{
  assert(startupOk);
  testInterface.feEvalL(supportLagrange.getElement1d(), I_test, J_test, supportLagrange.mapSurfaceY(Yd), val);
}

template<typename TEST, typename GEOSHAPE2D>
void
functionalLA2d<TEST,GEOSHAPE2D>::
evalGrad_test(const point3d & Yd, sVect<TEST_OUTTYPE> & gradX, sVect<TEST_OUTTYPE> & gradY, sVect<TEST_OUTTYPE> & gradZ)
{
  assert(startupOk);
  testInterface.feGradL(supportLagrange.getElement1d(), I_test, J_test, supportLagrange.mapSurfaceY(Yd), gradX, gradY, gradZ);
}


#endif
