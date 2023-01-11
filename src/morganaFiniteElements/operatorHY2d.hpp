/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef OPERATORLA2D_HPP
#define OPERATORLA2D_HPP

#include "elCardFeeder2d.hpp"
#include "supportHybrid2d.hpp"
#include "morganaFiniteElements.hpp"


/*! Operarator for hybrid FE 2d
<b> Scheme </b>
<ol>
<li> field 2d
<li> test  2d
<li> integration: internal faces 2d
</ol>
*/
template<typename FIELD, typename TEST>
class operatorHY2d
{
  /*! @name Typedefs */ //@{
  public:
    typedef typename FIELD::FEINTERFACE       FIELD_INTERFACE;
    typedef typename FIELD::FECARD            FIELD_FECARD;
    typedef typename FIELD::FECARDS           FIELD_FECARDS;
    typedef typename FIELD::OPTIONS           FIELD_OPTIONS;
    typedef typename FIELD::GEOSHAPE          FIELD_GEOSHAPE;
    typedef typename FIELD::PMAPTYPE          FIELD_PMAPTYPE;
    typedef typename FIELD::OUTTYPE           FIELD_OUTTYPE;
    typedef typename FIELD::FIELD_DOFTYPE     FIELD_DOFTYPE;
    typedef typename FIELD_INTERFACE::LISTMAP FIELD_LISTMAP;
    
    typedef typename TEST::FEINTERFACE        TEST_INTERFACE;
    typedef typename TEST::FECARD             TEST_FECARD;
    typedef typename TEST::FECARDS            TEST_FECARDS;
    typedef typename TEST::OPTIONS            TEST_OPTIONS;
    typedef typename TEST::GEOSHAPE           TEST_GEOSHAPE;
    typedef typename TEST::PMAPTYPE           TEST_PMAPTYPE;
    typedef typename TEST::OUTTYPE            TEST_OUTTYPE;
    typedef typename TEST::FIELD_DOFTYPE      TEST_DOFTYPE;
    typedef typename TEST_INTERFACE::LISTMAP  TEST_LISTMAP;
    
    typedef typename FIELD::MESH2D      MESH2D;
    typedef typename FIELD::CONNECT2D   CONNECT2D;
    typedef typename FIELD::GEOSHAPE    GEOSHAPE2D;

    typedef TEST_PMAPTYPE PMAPTYPE;
    typedef MESH2D        INTGRID;
    typedef MESH2D        LOOPGRID;
    
    typedef elCard2d<TEST_GEOSHAPE,PMAPTYPE>        ELCARD;
    typedef elCardFeeder2d<TEST_GEOSHAPE,PMAPTYPE>  ELCARDFEEDER;
    typedef supportHybrid2d<GEOSHAPE2D,PMAPTYPE>    SUPPORTHYBRID;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const OPClass opClass = opRowMajor;
    //@}
    
    /*! @name Internal flags */ //@{
  public:
    bool startupOk;
    bool commDevLoaded;
    bool geometryLoaded;
    //@}
    
     /*! @name Internal data */ //@{
  public:
    UInt el2d, locEd, I_test, J_test, I_field, J_field;
    Teuchos::RCP<MESH2D>     grid2d;
    Teuchos::RCP<CONNECT2D>  connectGrid2d;
    FIELD_INTERFACE fieldInterface;
    TEST_INTERFACE  testInterface;
    ELCARDFEEDER    elCardFeeder;
    SUPPORTHYBRID   supportHy;
    //@}
    
    /*! @name Constructors and setting - EXTERNAL */ //@{
  public:
    operatorHY2d();
    operatorHY2d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    operatorHY2d(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d);
    virtual ~operatorHY2d() {};
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    void setGeometry(MESH2D & Grid2d, CONNECT2D & ConnectGrid2d);
    //@}
    
     /*! @name Set map and startup - EXTERNAL */ //@{
  public:
    void setFeCardL_field(const UInt & elL, const FIELD_FECARD & FeCards);
    void setFeCardG_field(const UInt & elG, const FIELD_FECARD & FeCards);
    void setFeCards_field(const FIELD_FECARDS & FECards);
    void setOptions_field(const Teuchos::RCP<FIELD_OPTIONS> & Options);
    void setOptions_field(FIELD_OPTIONS & Options);
    
    void setFeCardL_test(const UInt & elL, const TEST_FECARD & FeCards);
    void setFeCardG_test(const UInt & elG, const TEST_FECARD & FeCards);
    void setFeCards_test(const TEST_FECARDS & FECards);
    void setOptions_test(const Teuchos::RCP<TEST_OPTIONS> & Options);
    void setOptions_test(TEST_OPTIONS & Options);
    
    void startup();
    void startup(const FIELD & CloneField, const TEST & CloneTest);
    void startup(const Teuchos::RCP<FIELD> & CloneField, const Teuchos::RCP<TEST> & CloneTest);
    //@}
    
    /*! @name LocalMatrix - INTERFACE */ //@{
  public:
    void setElement(const UInt & El);
    void setIJ(const UInt & II_test, const UInt & JJ_test, const UInt & II_field, const UInt & JJ_field);
    const Teuchos::RCP<MESH2D> & getLoopGrid() const;
    
    void setLocEdge(const UInt & LocEd);
    //@}
    
    /*! @name Maps - INTERFACE */ //@{
  public:
    const UInt          & getNumDofsL_field() const;
    const UInt          & getNumDofsG_field() const;
    const UInt          & getSizeListL_field() const;
    const UInt          & getSizeListG_field() const;
    const FIELD_LISTMAP & getListMap_field() const;
    
    const UInt         & getNumDofsL_test() const;
    const UInt         & getNumDofsG_test() const;
    const UInt         & getSizeListL_test() const;
    const UInt         & getSizeListG_test() const;
    const TEST_LISTMAP & getListMap_test() const;
    //@}
    
    /*! @name Indices - INTERFACE */ //@{
  public:
    UInt numIndex_field();
    void indexL_field(sVect<UInt> & indices);
    void indexG_field(sVect<UInt> & indices);
    
    UInt numIndex_test();
    void indexL_test(sVect<UInt> & indices);
    void indexG_test(sVect<UInt> & indices);
    //@}
    
    /*! @name Integration - INTERFACE (boundary faces of the 3d mesh) */ //@{
  public:
    const FIELD_FECARD   & getFeCardL_field();
    const TEST_FECARD    & getFeCardL_test();
    ELCARD                 getElCard() const;
    const Teuchos::RCP<INTGRID>   & getIntegrationGrid() const;
    const UInt           & getElement() const;
    const sVect<point3d> & getGlobEdgePoints() const;
    //@}
    
    /*! @name Operator - INTERFACE */ //@{
  public:
    point3d mapVolumeY(const point3d & Ys) const;
    point3d computeNormal(const point3d & Ys) const;
    
    void eval_field(const point3d & Y, sVect<FIELD_OUTTYPE> & val);
    void evalGrad_field(const point3d & Y, sVect<FIELD_OUTTYPE> & gradX, sVect<FIELD_OUTTYPE> & gradY, sVect<FIELD_OUTTYPE> & gradZ);
    
    void eval_test(const point3d & Y, sVect<TEST_OUTTYPE> & val);
    void evalGrad_test(const point3d & Y, sVect<TEST_OUTTYPE> & gradX, sVect<TEST_OUTTYPE> & gradY, sVect<TEST_OUTTYPE> & gradZ);
    
    virtual void eval(const point3d & Ys, sVect<komplex> & mat) = 0;
    virtual void eval(const point3d & Ys, sVect<Real> & mat) = 0;
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS AND SETTINGS
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST>
operatorHY2d<FIELD,TEST>::
operatorHY2d()
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 2>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 2>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName  == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  startupOk      = false;
  commDevLoaded  = false;
  geometryLoaded = false;
  
  el2d  = 0;
  locEd = 1;
  
  I_test  = 0;
  J_test  = 0; 
  I_field = 0;
  J_field = 0;
}

template<typename FIELD, typename TEST>
operatorHY2d<FIELD,TEST>::
operatorHY2d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 2>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 2>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName  == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  grid2d        = Grid2d;
  connectGrid2d = ConnectGrid2d;
  
  fieldInterface.setCommunicator(CommDev);
  fieldInterface.setGeometry(Grid2d,ConnectGrid2d);
  
  testInterface.setCommunicator(CommDev);
  testInterface.setGeometry(Grid2d,ConnectGrid2d);
  
  elCardFeeder.setGeometry(Grid2d,ConnectGrid2d);
  supportHy.setGeometry(Grid2d,ConnectGrid2d);
  
  startupOk      = false;
  commDevLoaded  = true;
  geometryLoaded = true;
  
  el2d  = 0;
  locEd = 1;
  
  I_test  = 0;
  J_test  = 0; 
  I_field = 0;
  J_field = 0;
}

template<typename FIELD, typename TEST>
operatorHY2d<FIELD,TEST>::
operatorHY2d(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d)
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 2>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 2>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName  == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
  
  fieldInterface.setCommunicator(CommDev);
  fieldInterface.setGeometry(Grid2d,ConnectGrid2d);
  
  testInterface.setCommunicator(CommDev);
  testInterface.setGeometry(Grid2d,ConnectGrid2d);
  
  elCardFeeder.setGeometry(Grid2d,ConnectGrid2d);
  supportHy.setGeometry(Grid2d,ConnectGrid2d);
  
  startupOk      = false;
  commDevLoaded  = true;
  geometryLoaded = true;
  
  el2d  = 0;
  locEd = 1;
  
  I_test  = 0;
  J_test  = 0; 
  I_field = 0;
  J_field = 0;
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  fieldInterface.setCommunicator(CommDev);
  testInterface.setCommunicator(CommDev);
  
  commDevLoaded = true;
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setCommunicator(communicator & CommDev)
{
  fieldInterface.setCommunicator(CommDev);
  testInterface.setCommunicator(CommDev);
  
  commDevLoaded = true;
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  grid2d        = Grid2d;
  connectGrid2d = ConnectGrid2d;
  
  fieldInterface.setGeometry(Grid2d,ConnectGrid2d);
  testInterface.setGeometry(Grid2d,ConnectGrid2d);
  
  elCardFeeder.setGeometry(Grid2d,ConnectGrid2d);
  supportHy.setGeometry(Grid2d,ConnectGrid2d);
  
  geometryLoaded = true;
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setGeometry(MESH2D & Grid2d, CONNECT2D & ConnectGrid2d)
{
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
  
  fieldInterface.setGeometry(Grid2d,ConnectGrid2d);
  testInterface.setGeometry(Grid2d,ConnectGrid2d);
  
  elCardFeeder.setGeometry(Grid2d,ConnectGrid2d);
  supportHy.setGeometry(Grid2d,ConnectGrid2d);
  
  geometryLoaded = true;
}



//_________________________________________________________________________________________________
// SET UP AND STARTUP
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setFeCardL_field(const UInt & elL, const FIELD_FECARD & FeCards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  fieldInterface.setFeCardL(elL,FeCards);
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setFeCardG_field(const UInt & elG, const FIELD_FECARD & FeCards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  fieldInterface.setFeCardG(elG,FeCards);
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setFeCards_field(const FIELD_FECARDS & FECards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  fieldInterface.setFeCards(FECards);
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setOptions_field(const Teuchos::RCP<FIELD_OPTIONS> & Options)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  fieldInterface.setOptions(Options);
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setOptions_field(FIELD_OPTIONS & Options)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  fieldInterface.setOptions(Options);
}


template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setFeCardL_test(const UInt & elL, const TEST_FECARD & FeCards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCardL(elL,FeCards);
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setFeCardG_test(const UInt & elG, const TEST_FECARD & FeCards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCardG(elG,FeCards);
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setFeCards_test(const TEST_FECARDS & FECards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCards(FECards);
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setOptions_test(const Teuchos::RCP<TEST_OPTIONS> & Options)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setOptions(Options);
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setOptions_test(TEST_OPTIONS & Options)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setOptions(Options);
}
    
template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
startup()
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  startupOk = true;
  
  fieldInterface.startup();
  testInterface.startup();
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
startup(const FIELD & CloneField, const TEST & CloneTest)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  startupOk = true;
  
  fieldInterface.clone(CloneField.getDofMapper());
  testInterface.clone(CloneTest.getDofMapper());
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
startup(const Teuchos::RCP<FIELD> & CloneField, const Teuchos::RCP<TEST> & CloneTest)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  startupOk = true;
  
  fieldInterface.clone(CloneField->getDofMapper());
  testInterface.clone(CloneTest->getDofMapper());
}



//_________________________________________________________________________________________________
// LOCAL MATRIX
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setElement(const UInt & El)
{
  assert(startupOk);
  
  el2d = El;
  supportHy.setElement2d(El);
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setIJ(const UInt & II_test, const UInt & JJ_test, const UInt & II_field, const UInt & JJ_field)
{
  assert(startupOk);
  I_test  = II_test;
  J_test  = JJ_test;
  I_field = II_field;
  J_field = JJ_field;
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
setLocEdge(const UInt & LocEd)
{
  assert(startupOk);
  
  locEd = LocEd;
  supportHy.setLocalEdge(LocEd);
}

template<typename FIELD, typename TEST>
const Teuchos::RCP<typename operatorHY2d<FIELD,TEST>::MESH2D> &
operatorHY2d<FIELD,TEST>::
getLoopGrid() const
{
  return(grid2d);
}



//_________________________________________________________________________________________________
// MAPS
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST>
const UInt &
operatorHY2d<FIELD,TEST>::
getNumDofsL_field() const
{
  assert(startupOk);
  return(fieldInterface.getNumDofsL());
}

template<typename FIELD, typename TEST>
const UInt &
operatorHY2d<FIELD,TEST>::
getNumDofsG_field() const
{
  assert(startupOk);
  return(fieldInterface.getNumDofsG());
}

template<typename FIELD, typename TEST>
const UInt &
operatorHY2d<FIELD,TEST>::
getSizeListL_field() const
{
  assert(startupOk);
  return(fieldInterface.getSizeListL());
}

template<typename FIELD, typename TEST>
const UInt &
operatorHY2d<FIELD,TEST>::
getSizeListG_field() const
{
  assert(startupOk);
  return(fieldInterface.getSizeListG());
}

template<typename FIELD, typename TEST>
const typename operatorHY2d<FIELD,TEST>::FIELD_LISTMAP &
operatorHY2d<FIELD,TEST>::
getListMap_field() const
{
  assert(startupOk);
  return(fieldInterface.getListMap());
}
    

template<typename FIELD, typename TEST>
const UInt &
operatorHY2d<FIELD,TEST>::
getNumDofsL_test() const
{
  assert(startupOk);
  return(testInterface.getNumDofsL());
}

template<typename FIELD, typename TEST>
const UInt &
operatorHY2d<FIELD,TEST>::
getNumDofsG_test() const
{
  assert(startupOk);
  return(testInterface.getNumDofsG());
}

template<typename FIELD, typename TEST>
const UInt &
operatorHY2d<FIELD,TEST>::
getSizeListL_test() const
{
  assert(startupOk);
  return(testInterface.getSizeListL());
}

template<typename FIELD, typename TEST>
const UInt &
operatorHY2d<FIELD,TEST>::
getSizeListG_test() const
{
  assert(startupOk);
  return(testInterface.getSizeListG());
}

template<typename FIELD, typename TEST>
const typename operatorHY2d<FIELD,TEST>::TEST_LISTMAP &
operatorHY2d<FIELD,TEST>::
getListMap_test() const
{
  assert(startupOk);
  return(testInterface.getListMap());
}



//_________________________________________________________________________________________________
// INDICES
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST>
UInt
operatorHY2d<FIELD,TEST>::
numIndex_field()
{
  assert(startupOk);
  return(fieldInterface.feNumBasisL(el2d));
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
indexL_field(sVect<UInt> & indices)
{
  assert(startupOk);
  return(fieldInterface.feIndicesLL(el2d,I_field,J_field,indices));
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
indexG_field(sVect<UInt> & indices)
{
  assert(startupOk);
  return(fieldInterface.feIndicesLG(el2d,I_field,J_field,indices));
}
    

template<typename FIELD, typename TEST>
UInt
operatorHY2d<FIELD,TEST>::
numIndex_test()
{  
  assert(startupOk);
  return( testInterface.feNumBasisL(el2d) );
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
indexL_test(sVect<UInt> & indices)
{
  assert(startupOk);
  return(testInterface.feIndicesLL(el2d,I_test,J_test,indices));
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
indexG_test(sVect<UInt> & indices)
{
  assert(startupOk);
  return(testInterface.feIndicesLG(el2d,I_test,J_test,indices));
}



//_________________________________________________________________________________________________
// INTEGRATION
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST>
const typename operatorHY2d<FIELD,TEST>::FIELD_FECARD &
operatorHY2d<FIELD,TEST>::
getFeCardL_field()
{
  assert(startupOk);
  return(fieldInterface.getFeCardL(el2d));
}

template<typename FIELD, typename TEST>
const typename operatorHY2d<FIELD,TEST>::TEST_FECARD &
operatorHY2d<FIELD,TEST>::
getFeCardL_test()
{
  assert(startupOk);
  return(testInterface.getFeCardL(el2d));
}

template<typename FIELD, typename TEST>
typename operatorHY2d<FIELD,TEST>::ELCARD
operatorHY2d<FIELD,TEST>::
getElCard() const
{
  assert(startupOk);
  return(elCardFeeder.getCardLocal(el2d));
}

template<typename FIELD, typename TEST>
const Teuchos::RCP<typename operatorHY2d<FIELD,TEST>::INTGRID> &
operatorHY2d<FIELD,TEST>::
getIntegrationGrid() const
{
  assert(startupOk);
  return(grid2d);
}

template<typename FIELD, typename TEST>
const UInt &
operatorHY2d<FIELD,TEST>::
getElement() const
{
  assert(startupOk);
  return(el2d);
}

template<typename FIELD, typename TEST>
const sVect<point3d> &
operatorHY2d<FIELD,TEST>::
getGlobEdgePoints() const
{
  assert(startupOk);
  return(supportHy.getGlobEdgePoints());
}



//_________________________________________________________________________________________________
// OPERATOR
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST>
point3d
operatorHY2d<FIELD,TEST>::
mapVolumeY(const point3d & Ys) const
{
  assert(startupOk);
  return(supportHy.mapVolumeY(Ys));
}

template<typename FIELD, typename TEST>
point3d
operatorHY2d<FIELD,TEST>::
computeNormal(const point3d & Ys) const
{
  assert(startupOk);
  return(supportHy.computeNormal(Ys));
}
    
template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
eval_field(const point3d & Ys, sVect<FIELD_OUTTYPE> & val)
{
  assert(startupOk);
  fieldInterface.feEvalL(el2d, I_field, J_field, supportHy.mapVolumeY(Ys), val);
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
evalGrad_field(const point3d & Ys, sVect<FIELD_OUTTYPE> & gradX, sVect<FIELD_OUTTYPE> & gradY, sVect<FIELD_OUTTYPE> & gradZ)
{
  assert(startupOk);
  fieldInterface.feGradL(el2d, I_field, J_field, supportHy.mapVolumeY(Ys), gradX, gradY, gradZ);
}
    
template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
eval_test(const point3d & Ys, sVect<TEST_OUTTYPE> & val)
{
  assert(startupOk);
  testInterface.feEvalL(el2d, I_test, J_test, supportHy.mapVolumeY(Ys), val);
}

template<typename FIELD, typename TEST>
void
operatorHY2d<FIELD,TEST>::
evalGrad_test(const point3d & Ys, sVect<TEST_OUTTYPE> & gradX, sVect<TEST_OUTTYPE> & gradY, sVect<TEST_OUTTYPE> & gradZ)
{
  assert(startupOk);
  testInterface.feGradL(el2d, I_test, J_test, supportHy.mapVolumeY(Ys), gradX, gradY, gradZ);
}

#endif
