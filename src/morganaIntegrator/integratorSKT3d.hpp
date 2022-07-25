/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef NODEFEEDERSKT3D_HH_
#define NODEFEEDERSKT3D_HH_

#include "glRule.h"

#include "elCard3d.hpp"
#include "intStaticCards.hpp"
#include "geoMapInterface.hpp"

#include "morganaIntegrator.hpp"


/*! Integration dummy card */
class intCardSKT3d
{ };


/*! Standard 3d integrator */
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
class integratorSKT3d
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
    typedef intCardSKT3d                INTCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const intClass intFlag = intSKT3d;
    
    geoMapInterface<GEOSHAPE> geoInterface;
    //@}
    
    /*! @name Functions */ //@{
  public:
    integratorSKT3d();
    
    /*! Constructor */
    integratorSKT3d(const INTCARD & IntCard);
    
    /*! Set integration card */
    void setIntCard(const INTCARD & IntCard);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(const OPERATOR & op);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(const RCP<OPERATOR> & op);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
integratorSKT3d<OPERATOR,TYPE,PRECISION>::
integratorSKT3d()
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 3>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName      == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
integratorSKT3d<OPERATOR,TYPE,PRECISION>::
integratorSKT3d(const INTCARD & IntCard)
{
  assert(staticAssert<FIELD_GEOSHAPE::nDim == 3>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
  assert(staticAssert<FIELD_GEOSHAPE::geoName      == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<FIELD_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
void
integratorSKT3d<OPERATOR,TYPE,PRECISION>::
setIntCard(const INTCARD & IntCard)
{
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorSKT3d<OPERATOR,TYPE,PRECISION>::
integration(const OPERATOR & op)
{
  //Degree determination
  UInt rx = PRECISION * op.getFeCardL_field().getRx() * op.getFeCardL_test().getRx() * UInt(op.getFeCardL_field().getIsActive()) * UInt(op.getFeCardL_test().getIsActive());
  UInt ry = PRECISION * op.getFeCardL_field().getRy() * op.getFeCardL_test().getRy() * UInt(op.getFeCardL_field().getIsActive()) * UInt(op.getFeCardL_test().getIsActive());
  UInt rz = PRECISION * op.getFeCardL_field().getRz() * op.getFeCardL_test().getRz() * UInt(op.getFeCardL_field().getIsActive()) * UInt(op.getFeCardL_test().getIsActive());
  
  //Set toll
  Real toll = 1e-10;
  
  //Nodes
  sVect<Real>  xgl,  ygl,  zgl;
  sVect<Real> wxgl, wygl, wzgl;
  
  glRule ruleX( UInt((rx + 1) / 2.0) + 1 ); 
  glRule ruleY( UInt((ry + 1) / 2.0) + 1 );
  glRule ruleZ( UInt((rz + 1) / 2.0) + 1 );
  
  ruleX.compute(xgl,wxgl,toll);
  ruleY.compute(ygl,wygl,toll);
  ruleZ.compute(zgl,wzgl,toll);
  
  //Alloc
  UInt matSize = op.numIndex_field() * op.numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  
  Real jac;
  point3d Y;
  
  mat.resize(matSize);
  
  //Integration loop
  for(UInt ix = 1; ix <= xgl.size(); ++ix)
  {
    Y.setX(xgl(ix));
    
    for(UInt iy = 1; iy <= ygl.size(); ++iy)
    {
      Y.setY(ygl(iy));
      
      for(UInt iz = 1; iz <= zgl.size(); ++iz)
      {
	Y.setZ(zgl(iz));
	
	//Evaluate field
        op.eval(Y, bufMat);
	
	//Evaluate jacobian
        jac = geoInterface.getGradientDet(op.getElCard().getNodes(), Y);
	
	//Sum
        for(UInt j=1; j <= matSize; ++j)
        {
          mat(j) += bufMat(j) * jac * wxgl(ix) * wygl(iy) * wzgl(iz);
        }
      }
    }
  }
  
  return(mat);
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorSKT3d<OPERATOR,TYPE,PRECISION>::
integration(const RCP<OPERATOR> & op)
{
  //Degree determination
  UInt rx = PRECISION * op->getFeCardL_field().getRx() * op->getFeCardL_test().getRx() * UInt(op->getFeCardL_field().getIsActive()) * UInt(op->getFeCardL_test().getIsActive());
  UInt ry = PRECISION * op->getFeCardL_field().getRy() * op->getFeCardL_test().getRy() * UInt(op->getFeCardL_field().getIsActive()) * UInt(op->getFeCardL_test().getIsActive());
  UInt rz = PRECISION * op->getFeCardL_field().getRz() * op->getFeCardL_test().getRz() * UInt(op->getFeCardL_field().getIsActive()) * UInt(op->getFeCardL_test().getIsActive());
  
  //Set toll
  Real toll = 1e-10;
  
  //Nodes
  sVect<Real>  xgl,  ygl,  zgl;
  sVect<Real> wxgl, wygl, wzgl;
  
  glRule ruleX( UInt((rx + 1) / 2.0) + 1 ); 
  glRule ruleY( UInt((ry + 1) / 2.0) + 1 );
  glRule ruleZ( UInt((rz + 1) / 2.0) + 1 );
  
  ruleX.compute(xgl,wxgl,toll);
  ruleY.compute(ygl,wygl,toll);
  ruleZ.compute(zgl,wzgl,toll);
  
  //Alloc
  UInt matSize = op->numIndex_field() * op->numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  
  Real jac;
  point3d Y;
  
  mat.resize(matSize);
  
  //Integration loop
  for(UInt ix = 1; ix <= xgl.size(); ++ix)
  {
    Y.setX(xgl(ix));
    
    for(UInt iy = 1; iy <= ygl.size(); ++iy)
    {
      Y.setY(ygl(iy));
      
      for(UInt iz = 1; iz <= zgl.size(); ++iz)
      {
	Y.setZ(zgl(iz));
	
	//Evaluate field
        op->eval(Y, bufMat);
	
	//Evaluate jacobian
        jac = geoInterface.getGradientDet(op->getElCard().getNodes(), Y);
	
	//Sum
        for(UInt j=1; j <= matSize; ++j)
        {
          mat(j) += bufMat(j) * jac * wxgl(ix) * wygl(iy) * wzgl(iz);
        }
      }
    }
  }
 
  return(mat);
}

#endif
