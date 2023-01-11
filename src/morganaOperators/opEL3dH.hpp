/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef OPEL3DH_HPP
#define OPEL3DH_HPP

#include "morganaFiniteElements.hpp"
#include "operatorEL3d.hpp"
#include "../morganaDofs/komplex.h"


/*! Operator 3d : c * U * V */
template<typename FIELD, typename TEST, typename COEFF>
class opEL3dH : public operatorEL3d<FIELD,TEST>
{
    /*! @name Typedefs */ //@{
  public:
    typedef operatorEL3d<FIELD,TEST>  OPERATOR;
    
    typedef typename OPERATOR::FIELD_OUTTYPE FIELD_OUTTYPE;
    typedef typename OPERATOR::TEST_OUTTYPE  TEST_OUTTYPE;
    
    typedef typename OPERATOR::MESH3D     MESH3D;
    typedef typename OPERATOR::CONNECT3D  CONNECT3D;
    
    typedef typename TEST::GEOSHAPE   TEST_GEOSHAPE;
    typedef typename FIELD::GEOSHAPE  FIELD_GEOSHAPE;
    typedef typename COEFF::GEOSHAPE  COEFF_GEOSHAPE;
    typedef typename COEFF::OUTTYPE   COEFF_OUTTYPE;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    bool coeffOk;
    Teuchos::RCP<COEFF> coeff;
    //@}
    
    /*! @name Functions */ //@{
  public:
    opEL3dH();
    opEL3dH(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    opEL3dH(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
    void setCoeff(Teuchos::RCP<COEFF> & Coeff);
    void setCoeff(COEFF & Coeff);
    void eval(const point3d & Y, sVect<komplex> & mat);
    void eval(const point3d & Y, sVect<Real> & mat);
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST, typename COEFF>
opEL3dH<FIELD,TEST,COEFF>::
opEL3dH() : operatorEL3d<FIELD,TEST>()
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  assert(staticAssert<COEFF_GEOSHAPE::geoName == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<COEFF_GEOSHAPE::geoName == FIELD_GEOSHAPE::geoName>::returnValue); 
  
  coeffOk = false;
}

template<typename FIELD, typename TEST, typename COEFF>
opEL3dH<FIELD,TEST,COEFF>::
opEL3dH(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d) : operatorEL3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  assert(staticAssert<COEFF_GEOSHAPE::geoName == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<COEFF_GEOSHAPE::geoName == FIELD_GEOSHAPE::geoName>::returnValue); 
  
  coeffOk = false;
}

template<typename FIELD, typename TEST, typename COEFF>
opEL3dH<FIELD,TEST,COEFF>::
opEL3dH(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d) : operatorEL3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  assert(staticAssert<COEFF_GEOSHAPE::geoName == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<COEFF_GEOSHAPE::geoName == FIELD_GEOSHAPE::geoName>::returnValue); 
  
  coeffOk = false;
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL3dH<FIELD,TEST,COEFF>::
setCoeff(Teuchos::RCP<COEFF> & Coeff)
{
  coeffOk = true;
  coeff   = Coeff;
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL3dH<FIELD,TEST,COEFF>::
setCoeff(COEFF & Coeff)
{
  coeffOk = true;
  coeff   = Teuchos::rcpFromRef(Coeff);
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL3dH<FIELD,TEST,COEFF>::
eval(const point3d & Y, sVect<Real> & mat)
{
  typedef operatorEL3d<FIELD,TEST>  OPERATOR;
  
  assert(coeffOk);
  assert(mat.size() == (OPERATOR::numIndex_field() * OPERATOR::numIndex_test()) );
  
  //Field eval
  sVect<FIELD_OUTTYPE> eval_field(OPERATOR::numIndex_field());
  OPERATOR::eval_field(Y,eval_field);
  
  //Test eval
  sVect<TEST_OUTTYPE> eval_test(OPERATOR::numIndex_test());
  OPERATOR::eval_test(Y,eval_test);
  
  //Coeff eval
  Real val;
  coeff->evalL(OPERATOR::el,Y,val);
  
  //Matrixing
  UInt tot = 1;
  
  for(UInt i=1; i <= OPERATOR::numIndex_field(); ++i)
  {
    for(UInt j=1; j <= OPERATOR::numIndex_test(); ++j)
    {
      mat(tot) = point3d::dot( eval_field(i), eval_test(j) ) * val;
      ++tot;
    }
  }
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL3dH<FIELD,TEST,COEFF>::
eval(const point3d & Y, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Y,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}

#endif
