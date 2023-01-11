/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef OPEL3DO_HPP
#define OPEL3DO_HPP

#include "traitsBasic.h"
#include "morganaFiniteElements.hpp"
#include "operatorEL3d.hpp"
#include "../morganaDofs/komplex.h"


/*! Operator 3d : lam * div(U) * div(V) + mu * (grad(U) + grad(U)^T) : grad(V) */
template<typename FIELD, typename TEST, typename COEFF>
class opEL3dO : public operatorEL3d<FIELD,TEST>
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
    Teuchos::RCP<COEFF> coeffLam;
    Teuchos::RCP<COEFF> coeffMu;
    //@}
    
    /*! @name Functions */ //@{
  public:
    opEL3dO();
    opEL3dO(const Teuchos::RCP<communicator> & CommDev,
            const Teuchos::RCP<MESH3D>       & Grid3d,
            const Teuchos::RCP<CONNECT3D>    & ConnedGrid3d);
    
    opEL3dO(communicator & CommDev,
            MESH3D       & Grid3d,
            CONNECT3D    & ConnedGrid3d);
    
    void setCoeff(const Teuchos::RCP<COEFF> & CoeffLam,
                  const Teuchos::RCP<COEFF> & CoeffMu);
    
    void setCoeff(const COEFF & CoeffLam,
                  const COEFF & CoeffMu);
    
    void eval(const point3d & Y, sVect<komplex> & mat);
    void eval(const point3d & Y, sVect<Real> & mat);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST, typename COEFF>
opEL3dO<FIELD,TEST,COEFF>::
opEL3dO()
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  assert(staticAssert<COEFF_GEOSHAPE::geoName == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<COEFF_GEOSHAPE::geoName == FIELD_GEOSHAPE::geoName>::returnValue); 
  
  coeffOk = false;
}

template<typename FIELD, typename TEST, typename COEFF>
opEL3dO<FIELD,TEST,COEFF>::
opEL3dO(const Teuchos::RCP<communicator> & CommDev,
        const Teuchos::RCP<MESH3D>       & Grid3d,
        const Teuchos::RCP<CONNECT3D>    & ConnedGrid3d) : operatorEL3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  assert(staticAssert<COEFF_GEOSHAPE::geoName == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<COEFF_GEOSHAPE::geoName == FIELD_GEOSHAPE::geoName>::returnValue); 
  
  coeffOk = false;
}

template<typename FIELD, typename TEST, typename COEFF>
opEL3dO<FIELD,TEST,COEFF>::
opEL3dO(communicator & CommDev,
        MESH3D       & Grid3d,
        CONNECT3D    & ConnedGrid3d) : operatorEL3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d)
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
opEL3dO<FIELD,TEST,COEFF>::
setCoeff(const Teuchos::RCP<COEFF> & CoeffLam,
         const Teuchos::RCP<COEFF> & CoeffMu)
{
  coeffOk  = true;
  coeffLam = CoeffLam;
  coeffMu  = CoeffMu;
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL3dO<FIELD,TEST,COEFF>::
setCoeff(const COEFF & CoeffLam,
         const COEFF & CoeffMu)
{
  coeffOk  = true;
  coeffLam = Teuchos::rcpFromRef(CoeffLam);
  coeffMu  = Teuchos::rcpFromRef(CoeffMu);
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL3dO<FIELD,TEST,COEFF>::
eval(const point3d & Y, sVect<Real> & mat)
{
  typedef operatorEL3d<FIELD,TEST> OPERATOR;
  
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
  
  //Set tensors
  tensor3d gradU, gradUT, gradV;
  
  //Coeff eval
  Real valLam, valMu;
  coeffLam->evalL(OPERATOR::el,Y,valLam);
  coeffMu->evalL(OPERATOR::el,Y,valMu);
  
  //Matrixing
  UInt tot = 1;
  
  for(UInt i=1; i <= OPERATOR::numIndex_field(); ++i)
  {
    for(UInt j=1; j <= OPERATOR::numIndex_test(); ++j)
    {
      gradU.setRow(1,gradX_field(i));
      gradU.setRow(2,gradY_field(i));
      gradU.setRow(3,gradZ_field(i));
  
      gradUT.setCol(1,gradX_field(i));
      gradUT.setCol(2,gradY_field(i));
      gradUT.setCol(3,gradZ_field(i));
  
      gradV.setRow(1,gradX_test(j));
      gradV.setRow(2,gradY_test(j));
      gradV.setRow(3,gradZ_test(j));
      
      mat(tot) = valLam * (gradU.getFirstInvariant()  * gradV.getFirstInvariant()) +
                 valMu  * (gradV.scalarProduct(gradU) + gradV.scalarProduct(gradUT));
      ++tot;
    }
  }
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL3dO<FIELD,TEST,COEFF>::
eval(const point3d & Y, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Y,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}

#endif
