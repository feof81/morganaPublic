/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef OPEL2DB_HPP
#define OPEL2DB_HPP

#include "morganaFiniteElements.hpp"
#include "operatorEL2d.hpp"
#include "../morganaDofs/komplex.h"


/*! Operator 2d : - div(U) * v */
template<typename FIELD, typename TEST>
class opEL2dB : public operatorEL2d<FIELD,TEST>
{
    /*! @name Typedefs */ //@{
  public:
    typedef operatorEL2d<FIELD,TEST>  OPERATOR;
    
    typedef typename OPERATOR::FIELD_OUTTYPE FIELD_OUTTYPE;
    typedef typename OPERATOR::TEST_OUTTYPE  TEST_OUTTYPE;
    
    typedef typename OPERATOR::MESH2D     MESH2D;
    typedef typename OPERATOR::CONNECT2D  CONNECT2D;
    //@}
    
    /*! @name Functions */ //@{
  public:
    opEL2dB();
    opEL2dB(const Teuchos::RCP<communicator> & CommDev,
            const Teuchos::RCP<MESH2D>       & Grid2d,
            const Teuchos::RCP<CONNECT2D>    & ConnectGrid2d);
    
    opEL2dB(communicator & CommDev,
            MESH2D       & Grid2d,
            CONNECT2D    & ConnectGrid2d);
    
    void eval(const point3d & Y, sVect<komplex> & mat);
    void eval(const point3d & Y, sVect<Real> & mat);
    //@}
};



template<typename FIELD, typename TEST>
opEL2dB<FIELD,TEST>::
opEL2dB() : operatorEL2d<FIELD,TEST>()
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST>
opEL2dB<FIELD,TEST>::
opEL2dB(const Teuchos::RCP<communicator> & CommDev,
        const Teuchos::RCP<MESH2D>       & Grid2d,
        const Teuchos::RCP<CONNECT2D>    & ConnectGrid2d) 
      : operatorEL2d<FIELD,TEST>(CommDev,Grid2d,ConnectGrid2d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST>
opEL2dB<FIELD,TEST>::
opEL2dB(communicator & CommDev,
        MESH2D       & Grid2d,
        CONNECT2D    & ConnectGrid2d)
      : operatorEL2d<FIELD,TEST>(CommDev,Grid2d,ConnectGrid2d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST>
void
opEL2dB<FIELD,TEST>::
eval(const point3d & Y, sVect<Real> & mat)
{
  typedef operatorEL2d<FIELD,TEST>  OPERATOR;
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
opEL2dB<FIELD,TEST>::
eval(const point3d & Y, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Y,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}

#endif
