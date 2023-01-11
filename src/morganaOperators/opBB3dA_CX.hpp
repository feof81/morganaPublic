/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef OPBB3DA_CX_HPP
#define OPBB3DA_CX_HPP

#include "morganaFiniteElements.hpp"
#include "operatorBB3d.hpp"
#include "../morganaDofs/komplex.h"


/*! Operator 3d : coeff * u * v, complex version*/
template<typename FIELD, typename TEST>
class opBB3dA_CX : public operatorBB3d<FIELD,TEST>
{
    /*! @name Typedefs */ //@{
  public:
    typedef operatorBB3d<FIELD,TEST>  OPERATOR;
    
    typedef typename OPERATOR::FIELD_OUTTYPE FIELD_OUTTYPE;
    typedef typename OPERATOR::TEST_OUTTYPE  TEST_OUTTYPE;
    
    typedef typename OPERATOR::MESH2D     MESH2D;
    typedef typename OPERATOR::MESH3D     MESH3D;
    typedef typename OPERATOR::CONNECT2D  CONNECT2D;
    typedef typename OPERATOR::CONNECT3D  CONNECT3D;
    
    typedef typename TEST::FIELD_DOFTYPE  TEST_DOFTYPE;
    typedef typename FIELD::FIELD_DOFTYPE FIELD_DOFTYPE;
    //@}
    
    /*! @name Functions */ //@{
  public:
    Real coeff;
    set<UInt> geoIds2d;
    Teuchos::RCP<MESH2D> grid2d;
    //@}
    
    /*! @name Functions */ //@{
  public:
    opBB3dA_CX();
    opBB3dA_CX(const Teuchos::RCP<communicator> & CommDev,
            const Teuchos::RCP<MESH3D>       & Grid3d,
            const Teuchos::RCP<CONNECT3D>    & ConnedGrid3d,
            const Teuchos::RCP<MESH2D>       & Grid2d,
            const Teuchos::RCP<CONNECT2D>    & ConnedGrid2d);
    
    opBB3dA_CX(communicator & CommDev,
            MESH3D       & Grid3d,
            CONNECT3D    & ConnedGrid3d,
            MESH2D       & Grid2d,
            CONNECT2D    & ConnedGrid2d);
    
    void setCoeff(const Real & Coeff);
    void setGeoIds(const set<UInt> & GeoIds);
    void eval(const point3d & Yf, sVect<komplex> & mat);
    void eval(const point3d & Yf, sVect<Real> & mat);
    //@}
};


template<typename FIELD, typename TEST>
opBB3dA_CX<FIELD,TEST>::
opBB3dA_CX() : operatorBB3d<FIELD,TEST>()
{
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<FIELD_DOFTYPE>::myType == typeKomplex>::returnValue);
  
  coeff = 1.0;
}

template<typename FIELD, typename TEST>
opBB3dA_CX<FIELD,TEST>::
opBB3dA_CX(const Teuchos::RCP<communicator> & CommDev,
	const Teuchos::RCP<MESH3D>       & Grid3d,
	const Teuchos::RCP<CONNECT3D>    & ConnedGrid3d,
	const Teuchos::RCP<MESH2D>       & Grid2d,
	const Teuchos::RCP<CONNECT2D>    & ConnedGrid2d) 
: operatorBB3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<FIELD_DOFTYPE>::myType == typeKomplex>::returnValue);
  
  coeff  = 1.0;
  grid2d = Grid2d;
}

template<typename FIELD, typename TEST>
opBB3dA_CX<FIELD,TEST>::
opBB3dA_CX(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d) 
: operatorBB3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<FIELD_DOFTYPE>::myType == typeKomplex>::returnValue);
  
  coeff  = 1.0;
  grid2d = Teuchos::rcpFromRef(Grid2d);
}

template<typename FIELD, typename TEST>
void
opBB3dA_CX<FIELD,TEST>::
setCoeff(const Real & Coeff)
{
  coeff = Coeff;
}

template<typename FIELD, typename TEST>
void
opBB3dA_CX<FIELD,TEST>::
setGeoIds(const set<UInt> & GeoIds)
{
  geoIds2d = GeoIds;
}

template<typename FIELD, typename TEST>
void
opBB3dA_CX<FIELD,TEST>::
eval(const point3d & Yf, sVect<komplex> & mat)
{
  typedef operatorBB3d<FIELD,TEST>          OPERATOR;
  typedef traitsMultiply<komplex,komplex> MULTIPLIER;
  
  assert(mat.size() == (OPERATOR::numIndex_field() * OPERATOR::numIndex_test()) );
  
  //Field eval
  sVect<FIELD_OUTTYPE> eval_field(OPERATOR::numIndex_field());
  OPERATOR::eval_field(Yf,eval_field);
  
  //Test eval
  sVect<TEST_OUTTYPE> eval_test(OPERATOR::numIndex_test());
  OPERATOR::eval_test(Yf,eval_test);
  
  //Element active
  UInt el2d   = OPERATOR::getEl2d();
  UInt geoId  = grid2d->getElementL(el2d).getGeoId();
  Real active = Real(geoIds2d.count(geoId) == 1);
  
  //Matrixing
  UInt tot = 1;
  komplex V(Real(OPERATOR::I_test == 1), -Real(OPERATOR::I_test == 2));
  
  for(UInt j=1; j <= OPERATOR::numIndex_test(); ++j)
  {
    for(UInt i=1; i <= OPERATOR::numIndex_field(); ++i)
    {  
      mat(tot) = MULTIPLIER::multiply(eval_field(i), eval_test(j) * V) * active * coeff;
      ++tot;
    }
  }
}

template<typename FIELD, typename TEST>
void
opBB3dA_CX<FIELD,TEST>::
eval(const point3d & Yf, sVect<Real> & mat)
{
  sVect<komplex> matC(mat.size());
  eval(Yf,matC);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex::cComp(matC(i),OPERATOR::I_test); }
}

#endif
