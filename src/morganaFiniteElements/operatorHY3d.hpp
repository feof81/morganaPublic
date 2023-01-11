/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef OPERATORHY3D_HPP
#define OPERATORHY3D_HPP

#include "elCardFeeder2d.hpp"
#include "supportHybrid3d.hpp"
#include "morganaFiniteElements.hpp"


/*! Operarator for hybrid FE 3d
<b> Scheme </b>
<ol>
<li> field 3d
<li> test  3d
<li> integration: internal faces 3d
</ol>
*/
template<typename FIELD, typename TEST>
class operatorHY3d
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
    
    typedef typename FIELD::MESH3D         MESH3D;
    typedef typename FIELD::CONNECT3D      CONNECT3D;
    typedef typename FIELD::GEOSHAPE       GEOSHAPE3D;
    typedef typename GEOSHAPE3D::GEOBSHAPE GEOSHAPE2D;

    typedef TEST_PMAPTYPE PMAPTYPE;
    typedef MESH3D        INTGRID;
    typedef MESH3D        LOOPGRID;
    
    typedef elCard3d<TEST_GEOSHAPE,PMAPTYPE>        ELCARD;
    typedef elCardFeeder3d<TEST_GEOSHAPE,PMAPTYPE>  ELCARDFEEDER;
    typedef supportHybrid3d<GEOSHAPE3D,PMAPTYPE>    SUPPORTHYBRID;
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
    UInt el3d, locFc, I_test, J_test, I_field, J_field;
    Teuchos::RCP<MESH3D>     grid3d;
    Teuchos::RCP<CONNECT3D>  connectGrid3d;
    FIELD_INTERFACE fieldInterface;
    TEST_INTERFACE  testInterface;
    ELCARDFEEDER    elCardFeeder;
    SUPPORTHYBRID   supportHy;
    //@}
    
    /*! @name Constructors and setting - EXTERNAL */ //@{
  public:
    operatorHY3d();
    operatorHY3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d);
    operatorHY3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnectGrid3d);
    virtual ~operatorHY3d() {};
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d);
    void setGeometry(MESH3D & Grid3d, CONNECT3D & ConnectGrid3d);
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
    const Teuchos::RCP<MESH3D> & getLoopGrid() const;
    
    void setLocFace(const UInt & LocFc);
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
    const FIELD_FECARD          & getFeCardL_field();
    const TEST_FECARD           & getFeCardL_test();
    ELCARD                        getElCard() const;
    const Teuchos::RCP<INTGRID> & getIntegrationGrid() const;
    const UInt                  & getElement() const;
    const sVect<point3d>        & getGlobFacePoints() const;
    //@}
    
    /*! @name Operator - INTERFACE */ //@{
  public:
    point3d mapVolumeY(const point3d & Yf) const;
    point3d computeNormal(const point3d & Yf) const;
    
    void eval_field(const point3d & Y, sVect<FIELD_OUTTYPE> & val);
    void evalGrad_field(const point3d & Y, sVect<FIELD_OUTTYPE> & gradX, sVect<FIELD_OUTTYPE> & gradY, sVect<FIELD_OUTTYPE> & gradZ);
    
    void eval_test(const point3d & Y, sVect<TEST_OUTTYPE> & val);
    void evalGrad_test(const point3d & Y, sVect<TEST_OUTTYPE> & gradX, sVect<TEST_OUTTYPE> & gradY, sVect<TEST_OUTTYPE> & gradZ);
    
    virtual void eval(const point3d & Yf, sVect<komplex> & mat) = 0;
    virtual void eval(const point3d & Yf, sVect<Real> & mat) = 0;
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS AND SETTINGS
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST>
operatorHY3d<FIELD,TEST>::
operatorHY3d()
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 3>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName  == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  startupOk      = false;
  commDevLoaded  = false;
  geometryLoaded = false;
  
  el3d  = 0;
  locFc = 1;
  
  I_test  = 0;
  J_test  = 0; 
  I_field = 0;
  J_field = 0;
}

template<typename FIELD, typename TEST>
operatorHY3d<FIELD,TEST>::
operatorHY3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d)
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 3>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName  == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  grid3d        = Grid3d;
  connectGrid3d = ConnectGrid3d;
  
  fieldInterface.setCommunicator(CommDev);
  fieldInterface.setGeometry(Grid3d,ConnectGrid3d);
  
  testInterface.setCommunicator(CommDev);
  testInterface.setGeometry(Grid3d,ConnectGrid3d);
  
  elCardFeeder.setGeometry(Grid3d,ConnectGrid3d);
  supportHy.setGeometry(Grid3d,ConnectGrid3d);
  
  startupOk      = false;
  commDevLoaded  = true;
  geometryLoaded = true;
  
  el3d  = 0;
  locFc = 1;
  
  I_test  = 0;
  J_test  = 0; 
  I_field = 0;
  J_field = 0;
}

template<typename FIELD, typename TEST>
operatorHY3d<FIELD,TEST>::
operatorHY3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnectGrid3d)
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 3>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName  == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  grid3d        = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid3d = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
  
  fieldInterface.setCommunicator(CommDev);
  fieldInterface.setGeometry(Grid3d,ConnectGrid3d);
  
  testInterface.setCommunicator(CommDev);
  testInterface.setGeometry(Grid3d,ConnectGrid3d);
  
  elCardFeeder.setGeometry(Grid3d,ConnectGrid3d);
  supportHy.setGeometry(Grid3d,ConnectGrid3d);
  
  startupOk      = false;
  commDevLoaded  = true;
  geometryLoaded = true;
  
  el3d  = 0;
  locFc = 1;
  
  I_test  = 0;
  J_test  = 0; 
  I_field = 0;
  J_field = 0;
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  fieldInterface.setCommunicator(CommDev);
  testInterface.setCommunicator(CommDev);
  
  commDevLoaded = true;
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
setCommunicator(communicator & CommDev)
{
  fieldInterface.setCommunicator(CommDev);
  testInterface.setCommunicator(CommDev);
  
  commDevLoaded = true;
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d)
{
  grid3d        = Grid3d;
  connectGrid3d = ConnectGrid3d;
  
  fieldInterface.setGeometry(Grid3d,ConnectGrid3d);
  testInterface.setGeometry(Grid3d,ConnectGrid3d);
  
  elCardFeeder.setGeometry(Grid3d,ConnectGrid3d);
  supportHy.setGeometry(Grid3d,ConnectGrid3d);
  
  geometryLoaded = true;
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
setGeometry(MESH3D & Grid3d, CONNECT3D & ConnectGrid3d)
{
  grid3d        = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid3d = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
  
  fieldInterface.setGeometry(Grid3d,ConnectGrid3d);
  testInterface.setGeometry(Grid3d,ConnectGrid3d);
  
  elCardFeeder.setGeometry(Grid3d,ConnectGrid3d);
  supportHy.setGeometry(Grid3d,ConnectGrid3d);
  
  geometryLoaded = true;
}



//_________________________________________________________________________________________________
// SET UP AND STARTUP
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
setFeCardL_field(const UInt & elL, const FIELD_FECARD & FeCards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  fieldInterface.setFeCardL(elL,FeCards);
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
setFeCardG_field(const UInt & elG, const FIELD_FECARD & FeCards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  fieldInterface.setFeCardG(elG,FeCards);
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
setFeCards_field(const FIELD_FECARDS & FECards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  fieldInterface.setFeCards(FECards);
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
setOptions_field(const Teuchos::RCP<FIELD_OPTIONS> & Options)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  fieldInterface.setOptions(Options);
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
setOptions_field(FIELD_OPTIONS & Options)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  fieldInterface.setOptions(Options);
}


template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
setFeCardL_test(const UInt & elL, const TEST_FECARD & FeCards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCardL(elL,FeCards);
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
setFeCardG_test(const UInt & elG, const TEST_FECARD & FeCards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCardG(elG,FeCards);
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
setFeCards_test(const TEST_FECARDS & FECards)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setFeCards(FECards);
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
setOptions_test(const Teuchos::RCP<TEST_OPTIONS> & Options)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setOptions(Options);
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
setOptions_test(TEST_OPTIONS & Options)
{
  assert(commDevLoaded);
  assert(geometryLoaded);
  
  testInterface.setOptions(Options);
}
    
template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
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
operatorHY3d<FIELD,TEST>::
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
operatorHY3d<FIELD,TEST>::
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
operatorHY3d<FIELD,TEST>::
setElement(const UInt & El)
{
  assert(startupOk);
  
  el3d = El;
  supportHy.setElement3d(El);
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
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
operatorHY3d<FIELD,TEST>::
setLocFace(const UInt & LocFc)
{
  assert(startupOk);
  
  locFc = LocFc;
  supportHy.setLocalFace(LocFc);
}

template<typename FIELD, typename TEST>
const Teuchos::RCP<typename operatorHY3d<FIELD,TEST>::MESH3D> &
operatorHY3d<FIELD,TEST>::
getLoopGrid() const
{
  return(grid3d);
}



//_________________________________________________________________________________________________
// MAPS
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST>
const UInt &
operatorHY3d<FIELD,TEST>::
getNumDofsL_field() const
{
  assert(startupOk);
  return(fieldInterface.getNumDofsL());
}

template<typename FIELD, typename TEST>
const UInt &
operatorHY3d<FIELD,TEST>::
getNumDofsG_field() const
{
  assert(startupOk);
  return(fieldInterface.getNumDofsG());
}

template<typename FIELD, typename TEST>
const UInt &
operatorHY3d<FIELD,TEST>::
getSizeListL_field() const
{
  assert(startupOk);
  return(fieldInterface.getSizeListL());
}

template<typename FIELD, typename TEST>
const UInt &
operatorHY3d<FIELD,TEST>::
getSizeListG_field() const
{
  assert(startupOk);
  return(fieldInterface.getSizeListG());
}

template<typename FIELD, typename TEST>
const typename operatorHY3d<FIELD,TEST>::FIELD_LISTMAP &
operatorHY3d<FIELD,TEST>::
getListMap_field() const
{
  assert(startupOk);
  return(fieldInterface.getListMap());
}
    

template<typename FIELD, typename TEST>
const UInt &
operatorHY3d<FIELD,TEST>::
getNumDofsL_test() const
{
  assert(startupOk);
  return(testInterface.getNumDofsL());
}

template<typename FIELD, typename TEST>
const UInt &
operatorHY3d<FIELD,TEST>::
getNumDofsG_test() const
{
  assert(startupOk);
  return(testInterface.getNumDofsG());
}

template<typename FIELD, typename TEST>
const UInt &
operatorHY3d<FIELD,TEST>::
getSizeListL_test() const
{
  assert(startupOk);
  return(testInterface.getSizeListL());
}

template<typename FIELD, typename TEST>
const UInt &
operatorHY3d<FIELD,TEST>::
getSizeListG_test() const
{
  assert(startupOk);
  return(testInterface.getSizeListG());
}

template<typename FIELD, typename TEST>
const typename operatorHY3d<FIELD,TEST>::TEST_LISTMAP &
operatorHY3d<FIELD,TEST>::
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
operatorHY3d<FIELD,TEST>::
numIndex_field()
{
  assert(startupOk);
  return(fieldInterface.feNumBasisL(el3d));
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
indexL_field(sVect<UInt> & indices)
{
  assert(startupOk);
  return(fieldInterface.feIndicesLL(el3d,I_field,J_field,indices));
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
indexG_field(sVect<UInt> & indices)
{
  assert(startupOk);
  return(fieldInterface.feIndicesLG(el3d,I_field,J_field,indices));
}
    

template<typename FIELD, typename TEST>
UInt
operatorHY3d<FIELD,TEST>::
numIndex_test()
{  
  assert(startupOk);
  return( testInterface.feNumBasisL(el3d) );
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
indexL_test(sVect<UInt> & indices)
{
  assert(startupOk);
  return(testInterface.feIndicesLL(el3d,I_test,J_test,indices));
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
indexG_test(sVect<UInt> & indices)
{
  assert(startupOk);
  return(testInterface.feIndicesLG(el3d,I_test,J_test,indices));
}



//_________________________________________________________________________________________________
// INTEGRATION
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST>
const typename operatorHY3d<FIELD,TEST>::FIELD_FECARD &
operatorHY3d<FIELD,TEST>::
getFeCardL_field()
{
  assert(startupOk);
  return(fieldInterface.getFeCardL(el3d));
}

template<typename FIELD, typename TEST>
const typename operatorHY3d<FIELD,TEST>::TEST_FECARD &
operatorHY3d<FIELD,TEST>::
getFeCardL_test()
{
  assert(startupOk);
  return(testInterface.getFeCardL(el3d));
}

template<typename FIELD, typename TEST>
typename operatorHY3d<FIELD,TEST>::ELCARD
operatorHY3d<FIELD,TEST>::
getElCard() const
{
  assert(startupOk);
  return(elCardFeeder.getCardLocal(el3d));
}

template<typename FIELD, typename TEST>
const Teuchos::RCP<typename operatorHY3d<FIELD,TEST>::INTGRID> &
operatorHY3d<FIELD,TEST>::
getIntegrationGrid() const
{
  assert(startupOk);
  return(grid3d);
}

template<typename FIELD, typename TEST>
const UInt &
operatorHY3d<FIELD,TEST>::
getElement() const
{
  assert(startupOk);
  return(el3d);
}

template<typename FIELD, typename TEST>
const sVect<point3d> &
operatorHY3d<FIELD,TEST>::
getGlobFacePoints() const
{
  assert(startupOk);
  return(supportHy.getGlobFacePoints());
}



//_________________________________________________________________________________________________
// OPERATOR
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST>
point3d
operatorHY3d<FIELD,TEST>::
mapVolumeY(const point3d & Yf) const
{
  assert(startupOk);
  return(supportHy.mapVolumeY(Yf));
}

template<typename FIELD, typename TEST>
point3d
operatorHY3d<FIELD,TEST>::
computeNormal(const point3d & Yf) const
{
  assert(startupOk);
  return(supportHy.computeNormal(Yf));
}
    
template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
eval_field(const point3d & Yf, sVect<FIELD_OUTTYPE> & val)
{
  assert(startupOk);
  fieldInterface.feEvalL(el3d, I_field, J_field, supportHy.mapVolumeY(Yf), val);
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
evalGrad_field(const point3d & Yf, sVect<FIELD_OUTTYPE> & gradX, sVect<FIELD_OUTTYPE> & gradY, sVect<FIELD_OUTTYPE> & gradZ)
{
  assert(startupOk);
  fieldInterface.feGradL(el3d, I_field, J_field, supportHy.mapVolumeY(Yf), gradX, gradY, gradZ);
}
    
template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
eval_test(const point3d & Yf, sVect<TEST_OUTTYPE> & val)
{
  assert(startupOk);
  testInterface.feEvalL(el3d, I_test, J_test, supportHy.mapVolumeY(Yf), val);
}

template<typename FIELD, typename TEST>
void
operatorHY3d<FIELD,TEST>::
evalGrad_test(const point3d & Yf, sVect<TEST_OUTTYPE> & gradX, sVect<TEST_OUTTYPE> & gradY, sVect<TEST_OUTTYPE> & gradZ)
{
  assert(startupOk);
  testInterface.feGradL(el3d, I_test, J_test, supportHy.mapVolumeY(Yf), gradX, gradY, gradZ);
}


#endif
