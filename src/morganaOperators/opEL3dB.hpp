/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef OPEL3DB_HPP
#define OPEL3DB_HPP

#include "morganaFiniteElements.hpp"
#include "operatorEL3d.hpp"
#include "../morganaDofs/komplex.h"


/*! Operator 3d : - div(U) * v */
template<typename FIELD, typename TEST>
class opEL3dB : public operatorEL3d<FIELD,TEST>
{
  /*! @name Typedefs */ //@{
  public:
    typedef operatorEL3d<FIELD,TEST>  OPERATOR;
    
    typedef typename OPERATOR::FIELD_OUTTYPE FIELD_OUTTYPE;
    typedef typename OPERATOR::TEST_OUTTYPE  TEST_OUTTYPE;
    
    typedef typename OPERATOR::MESH3D     MESH3D;
    typedef typename OPERATOR::CONNECT3D  CONNECT3D;
    //@}
    
    /*! @name Functions */ //@{
  public:
    opEL3dB();
    opEL3dB(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    opEL3dB(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
    void eval(const point3d & Y, sVect<komplex> & mat);
    void eval(const point3d & Y, sVect<Real> & mat);
    //@}
};



template<typename FIELD, typename TEST>
opEL3dB<FIELD,TEST>::
opEL3dB() : operatorEL3d<FIELD,TEST>()
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST>
opEL3dB<FIELD,TEST>::
opEL3dB(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d) : operatorEL3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST>
opEL3dB<FIELD,TEST>::
opEL3dB(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d) : operatorEL3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST>
void
opEL3dB<FIELD,TEST>::
eval(const point3d & Y, sVect<Real> & mat)
{
  typedef operatorEL3d<FIELD,TEST>  OPERATOR;
  assert(mat.size() == (OPERATOR::numIndex_field() * OPERATOR::numIndex_test()) );
  
  //Field eval
  sVect<FIELD_OUTTYPE> gradX_field(OPERATOR::numIndex_field());
  sVect<FIELD_OUTTYPE> gradY_field(OPERATOR::numIndex_field());
  sVect<FIELD_OUTTYPE> gradZ_field(OPERATOR::numIndex_field());
  OPERATOR::evalGrad_field(Y,gradX_field, gradY_field, gradZ_field);
  
  //Test eval
  sVect<TEST_OUTTYPE> val_test(OPERATOR::numIndex_test());
  OPERATOR::eval_test(Y,val_test);
  
  //Matrixing
  UInt tot = 1;
  
  for(UInt i=1; i <= OPERATOR::numIndex_field(); ++i)
  {
    for(UInt j=1; j <= OPERATOR::numIndex_test(); ++j)
    {     
      mat(tot) = - ( gradX_field(i).getX() + gradY_field(i).getY() + gradZ_field(i).getZ() ) * val_test(j);
      ++tot;
    }
  }
}

template<typename FIELD, typename TEST>
void
opEL3dB<FIELD,TEST>::
eval(const point3d & Y, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Y,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}

#endif
