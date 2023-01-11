/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef OPEL2DA_HPP
#define OPEL2DA_HPP

#include "morganaFiniteElements.hpp"
#include "operatorEL2d.hpp"
#include "../morganaDofs/komplex.h"


/*! Operator 3d : c * grad(u) * grad(v) */
template<typename FIELD, typename TEST, typename COEFF>
class opEL2dA : public operatorEL2d<FIELD,TEST>
{
    /*! @name Typedefs */ //@{
  public:
    typedef operatorEL2d<FIELD,TEST>  OPERATOR;
    
    typedef typename OPERATOR::FIELD_OUTTYPE FIELD_OUTTYPE;
    typedef typename OPERATOR::TEST_OUTTYPE  TEST_OUTTYPE;
    
    typedef typename OPERATOR::MESH2D     MESH2D;
    typedef typename OPERATOR::CONNECT2D  CONNECT2D;
    
    typedef typename TEST::GEOSHAPE   TEST_GEOSHAPE;
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
    opEL2dA();
    opEL2dA(const Teuchos::RCP<communicator> & CommDev,
	    const Teuchos::RCP<MESH2D>       & Grid2d,
	    const Teuchos::RCP<CONNECT2D>    & ConnedGrid2d);
    
    opEL2dA(communicator & CommDev,
	    MESH2D       & Grid2d,
	    CONNECT2D    & ConnedGrid2d);
    
    void setCoeff(const Teuchos::RCP<COEFF> & Coeff);
    void setCoeff(COEFF & Coeff);
    
    void eval(const point3d & Y, sVect<komplex> & mat);
    void eval(const point3d & Y, sVect<Real> & mat);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST, typename COEFF>
opEL2dA<FIELD,TEST,COEFF>::
opEL2dA() : operatorEL2d<FIELD,TEST>()
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<COEFF_GEOSHAPE::geoName == TEST_GEOSHAPE::geoName>::returnValue); 
  
  coeffOk = false;
}

template<typename FIELD, typename TEST, typename COEFF>
opEL2dA<FIELD,TEST,COEFF>::
opEL2dA(const Teuchos::RCP<communicator> & CommDev,
	const Teuchos::RCP<MESH2D>       & Grid2d,
	const Teuchos::RCP<CONNECT2D>    & ConnedGrid2d) : operatorEL2d<FIELD,TEST>(CommDev,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<COEFF_GEOSHAPE::geoName == TEST_GEOSHAPE::geoName>::returnValue);
  
  coeffOk = false;
}

template<typename FIELD, typename TEST, typename COEFF>
opEL2dA<FIELD,TEST,COEFF>::
opEL2dA(communicator & CommDev,
	MESH2D       & Grid2d,
	CONNECT2D    & ConnedGrid2d) : operatorEL2d<FIELD,TEST>(CommDev,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<COEFF_GEOSHAPE::geoName == TEST_GEOSHAPE::geoName>::returnValue);
  
  coeffOk = false;
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL2dA<FIELD,TEST,COEFF>::
setCoeff(const Teuchos::RCP<COEFF> & Coeff)
{
  coeffOk = true;
  coeff   = Coeff;
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL2dA<FIELD,TEST,COEFF>::
setCoeff(COEFF & Coeff)
{
  coeffOk = true;
  coeff   = Teuchos::rcpFromRef(Coeff);
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL2dA<FIELD,TEST,COEFF>::
eval(const point3d & Y, sVect<Real> & mat)
{
  typedef operatorEL2d<FIELD,TEST>  OPERATOR;
  
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
  Real val;
  coeff->evalL(OPERATOR::el,Y,val);
  
  //Matrixing
  UInt tot = 1;
  
  for(UInt i=1; i <= OPERATOR::numIndex_field(); ++i)
  {
    for(UInt j=1; j <= OPERATOR::numIndex_test(); ++j)
    {
      mat(tot) = ( gradX_field(i) * gradX_test(j) +
                   gradY_field(i) * gradY_test(j) + 
		   gradZ_field(i) * gradZ_test(j) ) * val;
      ++tot;
    }
  }
}

template<typename FIELD, typename TEST, typename COEFF>
void
opEL2dA<FIELD,TEST,COEFF>::
eval(const point3d & Y, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Y,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}


#endif
