/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef OPBB2DA_HPP
#define OPBB2DA_HPP

#include "morganaFiniteElements.hpp"
#include "operatorBB2d.hpp"
#include "../morganaDofs/komplex.h"


/*! Operator 2d : u * q */
template<typename FIELD, typename TEST>
class opBB2dA : public operatorBB2d<FIELD,TEST>
{
    /*! @name Typedefs */ //@{
  public:
    typedef operatorBB2d<FIELD,TEST>  OPERATOR;
    
    typedef typename OPERATOR::FIELD_OUTTYPE FIELD_OUTTYPE;
    typedef typename OPERATOR::TEST_OUTTYPE  TEST_OUTTYPE;
    
    typedef typename OPERATOR::MESH1D     MESH1D;
    typedef typename OPERATOR::MESH2D     MESH2D;
    typedef typename OPERATOR::CONNECT1D  CONNECT1D;
    typedef typename OPERATOR::CONNECT2D  CONNECT2D;
    //@}
    
    /*! @name Functions */ //@{
  public:
    Real coeff;
    set<UInt> geoIds1d;
    Teuchos::RCP<MESH1D> grid1d;
    //@}
    
    /*! @name Functions */ //@{
  public:
    opBB2dA();
    opBB2dA(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d);
    opBB2dA(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d, MESH1D & Grid1d, CONNECT1D & ConnedGrid1d);
    void setCoeff(const Real & Coeff);
    void setGeoIds(const set<UInt> & GeoIds);
    void eval(const point3d & Yf, sVect<komplex> & mat);
    void eval(const point3d & Yf, sVect<Real> & mat);
    //@}
};


template<typename FIELD, typename TEST>
opBB2dA<FIELD,TEST>::
opBB2dA() : operatorBB2d<FIELD,TEST>()
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
  
  coeff = 1.0;
}

template<typename FIELD, typename TEST>
opBB2dA<FIELD,TEST>::
opBB2dA(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d) 
: operatorBB2d<FIELD,TEST>(CommDev,Grid2d,ConnedGrid2d,Grid1d,ConnedGrid1d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
  
  coeff  = 1.0;
  grid1d = Grid1d;
}

template<typename FIELD, typename TEST>
opBB2dA<FIELD,TEST>::
opBB2dA(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d, MESH1D & Grid1d, CONNECT1D & ConnedGrid1d) 
: operatorBB2d<FIELD,TEST>(CommDev,Grid2d,ConnedGrid2d,Grid1d,ConnedGrid1d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
  
  coeff  = 1.0;
  grid1d = Teuchos::rcpFromRef(Grid1d);
}

template<typename FIELD, typename TEST>
void
opBB2dA<FIELD,TEST>::
setCoeff(const Real & Coeff)
{
  coeff = Coeff;
}

template<typename FIELD, typename TEST>
void
opBB2dA<FIELD,TEST>::
setGeoIds(const set<UInt> & GeoIds)
{
  geoIds1d = GeoIds;
}

template<typename FIELD, typename TEST>
void
opBB2dA<FIELD,TEST>::
eval(const point3d & Yf, sVect<Real> & mat)
{
  typedef operatorBB2d<FIELD,TEST>  OPERATOR;
  
  assert(mat.size() == (OPERATOR::numIndex_field() * OPERATOR::numIndex_test()) );
  
  //Field eval
  sVect<FIELD_OUTTYPE> eval_field(OPERATOR::numIndex_field());
  OPERATOR::eval_field(Yf,eval_field);
  
  //Test eval
  sVect<TEST_OUTTYPE> eval_test(OPERATOR::numIndex_test());
  OPERATOR::eval_test(Yf,eval_test);
  
  //Element active
  UInt el1d   = OPERATOR::getEl1d();
  UInt geoId  = grid1d->getElementL(el1d).getGeoId();
  Real active = Real(geoIds1d.count(geoId) == 1);
  
  //Matrixing
  UInt tot = 1;
  
  for(UInt j=1; j <= OPERATOR::numIndex_test(); ++j)
  {
    for(UInt i=1; i <= OPERATOR::numIndex_field(); ++i)
    {  
      mat(tot) = active * coeff * ( eval_field(i) * eval_test(j) );
      ++tot;
    }
  }
}

template<typename FIELD, typename TEST>
void
opBB2dA<FIELD,TEST>::
eval(const point3d & Yf, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Yf,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}


#endif
