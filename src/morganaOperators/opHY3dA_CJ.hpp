/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef OPHY3DA_CJ_HPP
#define OPHY3DA_CJ_HPP

#include "morganaFiniteElements.hpp"
#include "operatorHY3d.hpp"
#include "../morganaDofs/komplex.h"


/*! Operator hybrid 3d complex : U * lam */
template<typename FIELD, typename TEST>
class opHY3dA_CJ : public operatorHY3d<FIELD,TEST>
{
    /*! @name Typedefs */ //@{
  public:
    typedef operatorHY3d<FIELD,TEST>  OPERATOR;
    
    typedef typename OPERATOR::FIELD_OUTTYPE FIELD_OUTTYPE;
    typedef typename OPERATOR::TEST_OUTTYPE  TEST_OUTTYPE;
    
    typedef typename OPERATOR::MESH3D     MESH3D;
    typedef typename OPERATOR::CONNECT3D  CONNECT3D;
    
    typedef typename TEST::FIELD_DOFTYPE  TEST_DOFTYPE;
    typedef typename FIELD::FIELD_DOFTYPE FIELD_DOFTYPE;
    //@}
    
    /*! @name Functions */ //@{
  public:
    opHY3dA_CJ();
    opHY3dA_CJ(const Teuchos::RCP<communicator> & CommDev,
               const Teuchos::RCP<MESH3D>       & Grid3d,
               const Teuchos::RCP<CONNECT3D>    & ConnedGrid3);
    
    opHY3dA_CJ(communicator & CommDev,
               MESH3D       & Grid3d,
               CONNECT3D    & ConnedGrid3d);
    
    void eval(const point3d & Yf, sVect<komplex> & mat);
    void eval(const point3d & Yf, sVect<Real> & mat);
    //@}
};

template<typename FIELD, typename TEST>
opHY3dA_CJ<FIELD,TEST>::
opHY3dA_CJ() : operatorHY3d<FIELD,TEST>()
{
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<FIELD_DOFTYPE>::myType  == typeKomplex>::returnValue);
}

template<typename FIELD, typename TEST>
opHY3dA_CJ<FIELD,TEST>::
opHY3dA_CJ(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d) : operatorHY3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<FIELD_DOFTYPE>::myType  == typeKomplex>::returnValue);
}

template<typename FIELD, typename TEST>
opHY3dA_CJ<FIELD,TEST>::
opHY3dA_CJ(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d) : operatorHY3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<FIELD_DOFTYPE>::myType  == typeKomplex>::returnValue);
}

template<typename FIELD, typename TEST>
void
opHY3dA_CJ<FIELD,TEST>::
eval(const point3d & Yf, sVect<komplex> & mat)
{
  typedef operatorHY3d<FIELD,TEST>  OPERATOR;
  
  assert(mat.size() == (OPERATOR::numIndex_field() * OPERATOR::numIndex_test()) );
  
  //Field eval
  sVect<FIELD_OUTTYPE> eval_field(OPERATOR::numIndex_field());
  OPERATOR::eval_field(Yf,eval_field);
  
  //Test eval
  sVect<TEST_OUTTYPE> eval_test(OPERATOR::numIndex_test());
  OPERATOR::eval_test(Yf,eval_test);
  
  //Matrixing
  UInt tot = 1;
  komplex V(Real(OPERATOR::I_test == 1), -Real(OPERATOR::I_test == 2));
  
  for(UInt j=1; j <= OPERATOR::numIndex_test(); ++j)
  {
    for(UInt i=1; i <= OPERATOR::numIndex_field(); ++i)
    {
      mat(tot) = komplex::conj(V * eval_field(i) * eval_test(j));
      ++tot;
    }
  }
}

template<typename FIELD, typename TEST>
void
opHY3dA_CJ<FIELD,TEST>::
eval(const point3d & Yf, sVect<Real> & mat)
{
  sVect<komplex> matC(mat.size());
  eval(Yf,matC);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex::cComp(matC(i), OPERATOR::I_test); }
}

#endif
