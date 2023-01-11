/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef OPEL3DN_HPP
#define OPEL3DN_HPP

#include "morganaTypes.hpp"
#include "morganaFiniteElements.hpp"
#include "operatorEL3d.hpp"
#include "../morganaDofs/tensor2d.h"
#include "../morganaDofs/komplex.h"


/*! Operator 3d : C * grad(U) * grad(U)
C : 2x2 tensor
U : 2 point solution vector
V : 2 point test function */
template<typename FIELD, typename TEST, typename COEFF>
class opEL3dN : public operatorEL3d<FIELD,TEST>
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
    opEL3dN();
    opEL3dN(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    opEL3dN(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
    void setCoeff(const Teuchos::RCP<COEFF> & Coeff);
    void setCoeff(const COEFF & Coeff);
    void eval(const point3d & Y, sVect<komplex> & mat);
    void eval(const point3d & Y, sVect<Real> & mat);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST, typename COEFF>
opEL3dN<FIELD,TEST,COEFF>::
opEL3dN() : operatorEL3d<FIELD,TEST>()
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint2d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint2d>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeTensor2d>::returnValue);
  
  assert(staticAssert<COEFF_GEOSHAPE::geoName == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<COEFF_GEOSHAPE::geoName == FIELD_GEOSHAPE::geoName>::returnValue); 
  
  coeffOk = false;
}

template<typename FIELD, typename TEST, typename COEFF>
opEL3dN<FIELD,TEST,COEFF>::
opEL3dN(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d) : operatorEL3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint2d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint2d>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeTensor2d>::returnValue);
  assert(staticAssert<COEFF_GEOSHAPE::geoName == TEST_GEOSHAPE::geoName>::returnValue);
  
  coeffOk = false;
}

template<typename FIELD, typename TEST, typename COEFF>
opEL3dN<FIELD,TEST,COEFF>::
opEL3dN(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d) : operatorEL3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint2d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint2d>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeTensor2d>::returnValue);
  
  assert(staticAssert<COEFF_GEOSHAPE::geoName == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<COEFF_GEOSHAPE::geoName == FIELD_GEOSHAPE::geoName>::returnValue); 
  
  coeffOk = false;
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL3dN<FIELD,TEST,COEFF>::
setCoeff(const Teuchos::RCP<COEFF> & Coeff)
{
  coeffOk = true;
  coeff   = Coeff;
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL3dN<FIELD,TEST,COEFF>::
setCoeff(const COEFF & Coeff)
{
  coeffOk = true;
  coeff   = Teuchos::rcpFromRef(Coeff);
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL3dN<FIELD,TEST,COEFF>::
eval(const point3d & Y, sVect<Real> & mat)
{
  typedef operatorEL3d<FIELD,TEST>  OPERATOR;
  
  assert(coeffOk);
  assert(mat.size() == (OPERATOR::numIndex_field() * OPERATOR::numIndex_test()) );
  
  //Field eval
  sVect<FIELD_OUTTYPE> gradX_field(OPERATOR::numIndex_field());
  sVect<FIELD_OUTTYPE> gradY_field(OPERATOR::numIndex_field());
  sVect<FIELD_OUTTYPE> gradZ_field(OPERATOR::numIndex_field());
  OPERATOR::evalGrad_field(Y,gradX_field,gradY_field,gradZ_field);
  
  //Test eval
  sVect<TEST_OUTTYPE> gradX_test(OPERATOR::numIndex_test());
  sVect<TEST_OUTTYPE> gradY_test(OPERATOR::numIndex_test());
  sVect<TEST_OUTTYPE> gradZ_test(OPERATOR::numIndex_test());
  OPERATOR::evalGrad_test(Y,gradX_test,gradY_test,gradZ_test);
  
  //Coeff eval
  tensor2d val;
  coeff->evalL(OPERATOR::el,Y,val);
  
  //Matrixing
  UInt tot = 1;
  
  for(UInt i=1; i <= OPERATOR::numIndex_field(); ++i)
  {
    for(UInt j=1; j <= OPERATOR::numIndex_test(); ++j)
    {
      mat(tot) = gradX_field(i).getX() * val.getIJ(1,1) * gradX_test(j).getX() +
                 gradX_field(i).getX() * val.getIJ(1,2) * gradX_test(j).getY() +
                 gradX_field(i).getY() * val.getIJ(2,1) * gradX_test(j).getX() +
                 gradX_field(i).getY() * val.getIJ(2,2) * gradX_test(j).getY() +
      
                 gradY_field(i).getX() * val.getIJ(1,1) * gradY_test(j).getX() +
                 gradY_field(i).getX() * val.getIJ(1,2) * gradY_test(j).getY() +
                 gradY_field(i).getY() * val.getIJ(2,1) * gradY_test(j).getX() +
                 gradY_field(i).getY() * val.getIJ(2,2) * gradY_test(j).getY() + 
 
                 gradZ_field(i).getX() * val.getIJ(1,1) * gradZ_test(j).getX() +
                 gradZ_field(i).getX() * val.getIJ(1,2) * gradZ_test(j).getY() +
                 gradZ_field(i).getY() * val.getIJ(2,1) * gradZ_test(j).getX() +
                 gradZ_field(i).getY() * val.getIJ(2,2) * gradZ_test(j).getY();
      ++tot;
    }
  }
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL3dN<FIELD,TEST,COEFF>::
eval(const point3d & Y, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Y,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}

#endif
