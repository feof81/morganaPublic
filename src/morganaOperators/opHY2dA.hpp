/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef OPHY3DA_HPP
#define OPHY3DA_HPP

#include "morganaFiniteElements.hpp"
#include "operatorHY2d.hpp"
#include "../morganaDofs/komplex.h"


/*! Operator hybrid 2d : (U * N) * lam */
template<typename FIELD, typename TEST>
class opHY2dA : public operatorHY2d<FIELD,TEST>
{
  /*! @name Typedefs */ //@{
  public:
    typedef operatorHY2d<FIELD,TEST>  OPERATOR;
    
    typedef typename OPERATOR::FIELD_OUTTYPE FIELD_OUTTYPE;
    typedef typename OPERATOR::TEST_OUTTYPE  TEST_OUTTYPE;
    
    typedef typename OPERATOR::MESH2D     MESH2D;
    typedef typename OPERATOR::CONNECT2D  CONNECT2D;
    //@}
    
    /*! @name Functions */ //@{
  public:
    opHY2dA();
    opHY2dA(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d);
    opHY2dA(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d);
    void eval(const point3d & Ys, sVect<komplex> & mat);
    void eval(const point3d & Ys, sVect<Real> & mat);
    //@}
};


template<typename FIELD, typename TEST>
opHY2dA<FIELD,TEST>::
opHY2dA() : operatorHY2d<FIELD,TEST>()
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST>
opHY2dA<FIELD,TEST>::
opHY2dA(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d) : operatorHY2d<FIELD,TEST>(CommDev,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST>
opHY2dA<FIELD,TEST>::
opHY2dA(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d) : operatorHY2d<FIELD,TEST>(CommDev,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST>
void
opHY2dA<FIELD,TEST>::
eval(const point3d & Ys, sVect<Real> & mat)
{
  typedef operatorHY2d<FIELD,TEST>  OPERATOR;
  
  assert(mat.size() == (OPERATOR::numIndex_field() * OPERATOR::numIndex_test()) );
  
  //Field eval
  sVect<FIELD_OUTTYPE> eval_field(OPERATOR::numIndex_field());
  OPERATOR::eval_field(Ys,eval_field);
  
  //Test eval
  sVect<TEST_OUTTYPE> eval_test(OPERATOR::numIndex_test());
  OPERATOR::eval_test(Ys,eval_test);
  
  //Normal eval
  point3d N = OPERATOR::computeNormal(Ys);
  
  //Matrixing
  UInt tot = 1;
  
  for(UInt j=1; j <= OPERATOR::numIndex_test(); ++j)
  {
    for(UInt i=1; i <= OPERATOR::numIndex_field(); ++i)
    {
      mat(tot) = ( point3d::dot(eval_field(i), N) * eval_test(j) );
      ++tot;
    }
  }
}

template<typename FIELD, typename TEST>
void
opHY2dA<FIELD,TEST>::
eval(const point3d & Ys, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Ys,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}

#endif
