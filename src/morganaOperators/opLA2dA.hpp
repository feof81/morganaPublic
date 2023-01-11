/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef OPLA2DA_HPP
#define OPLA2DA_HPP

#include "morganaFiniteElements.hpp"
#include "operatorLA2d.hpp"
#include "../morganaDofs/komplex.h"


/*! Operator 2d : u * q */
template<typename FIELD, typename TEST>
class opLA2dA : public operatorLA2d<FIELD,TEST>
{
    /*! @name Typedefs */ //@{
  public:
    typedef operatorLA2d<FIELD,TEST>  OPERATOR;
    
    typedef typename OPERATOR::FIELD_OUTTYPE FIELD_OUTTYPE;
    typedef typename OPERATOR::TEST_OUTTYPE  TEST_OUTTYPE;
    
    typedef typename OPERATOR::MESH1D     MESH1D;
    typedef typename OPERATOR::MESH2D     MESH2D;
    typedef typename OPERATOR::CONNECT1D  CONNECT1D;
    typedef typename OPERATOR::CONNECT2D  CONNECT2D;
    
    typedef typename TEST::GEOSHAPE   GEOSHAPE2D;
    typedef typename TEST::GEOSHAPE   TEST_GEOSHAPE;
    typedef typename TEST::PMAPTYPE   TEST_PMAPTYPE;
    
    typedef typename FIELD::GEOSHAPE  GEOSHAPE3D;
    typedef typename FIELD::GEOSHAPE  FIELD_GEOSHAPE;
    typedef typename FIELD::PMAPTYPE  FIELD_PMAPTYPE;
    //@}
    
    /*! @name Functions */ //@{
  public:
    opLA2dA();
    opLA2dA(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d);
    opLA2dA(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d, MESH1D & Grid1d, CONNECT1D & ConnedGrid1d);
    void eval(const point3d & Yd, sVect<komplex> & mat);
    void eval(const point3d & Yd, sVect<Real> & mat);
    //@}
};


template<typename FIELD, typename TEST>
opLA2dA<FIELD,TEST>::
opLA2dA() : operatorLA2d<FIELD,TEST>()
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST>
opLA2dA<FIELD,TEST>::
opLA2dA(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d) : operatorLA2d<FIELD,TEST>(CommDev,Grid2d,ConnedGrid2d,Grid1d,ConnedGrid1d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST>
opLA2dA<FIELD,TEST>::
opLA2dA(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d, MESH1D & Grid1d, CONNECT1D & ConnedGrid1d) : operatorLA2d<FIELD,TEST>(CommDev,Grid2d,ConnedGrid2d,Grid1d,ConnedGrid1d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST>
void
opLA2dA<FIELD,TEST>::
eval(const point3d & Yd, sVect<Real> & mat)
{
  typedef operatorLA2d<FIELD,TEST>  OPERATOR;
  
  assert(mat.size() == (OPERATOR::numIndex_field() * OPERATOR::numIndex_test()) );
  
  //Field eval
  sVect<FIELD_OUTTYPE> eval_field(OPERATOR::numIndex_field());
  OPERATOR::eval_field(Yd,eval_field);
  
  //Test eval
  sVect<TEST_OUTTYPE> eval_test(OPERATOR::numIndex_test());
  OPERATOR::eval_test(Yd,eval_test);
  
  //Matrixing
  UInt tot = 1;
  
  for(UInt j=1; j <= OPERATOR::numIndex_test(); ++j)
  {
    for(UInt i=1; i <= OPERATOR::numIndex_field(); ++i)
    {  
      mat(tot) = ( eval_field(i) * eval_test(j) );
      ++tot;
    }
  }
}

template<typename FIELD, typename TEST>
void
opLA2dA<FIELD,TEST>::
eval(const point3d & Yd, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Yd,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}


#endif
