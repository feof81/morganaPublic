/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FNLA2DA_HPP
#define FNLA2DA_HPP
#include "../morganaFiniteElements/functionalLA2d.hpp"
#include "../morganaDofs/morganaTypes.hpp"
#include "../morganaDofs/point3d.h"
#include "../morganaDofs/komplex.h"

template<typename TEST, typename GEOSHAPE2D, typename COEFF>
class fnLA2dA : public functionalLA2d<TEST,GEOSHAPE2D>
{
    /*! @name Typedefs */ //@{
  public:
    typedef functionalLA2d<TEST,GEOSHAPE2D>    FUNCTIONAL;
    
    typedef typename FUNCTIONAL::TEST_OUTTYPE  TEST_OUTTYPE;
    typedef typename FUNCTIONAL::MESH1D        MESH1D;
    typedef typename FUNCTIONAL::MESH2D        MESH2D;
    typedef typename FUNCTIONAL::CONNECT1D     CONNECT1D;
    typedef typename FUNCTIONAL::CONNECT2D     CONNECT2D;
    
    typedef typename TEST::PMAPTYPE   TEST_PMAPTYPE;
    typedef typename TEST::GEOSHAPE   TEST_GEOSHAPE;
    
    typedef typename COEFF::GEOSHAPE  COEFF_GEOSHAPE;
    typedef typename COEFF::OUTTYPE   COEFF_OUTTYPE;
    typedef typename COEFF::PMAPTYPE  COEFF_PMAPTYPE;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    bool coeffOk;
    Teuchos::RCP<COEFF> coeff;
    //@}
    
    /*! @name Functions */ //@{
  public:
    fnLA2dA();
    fnLA2dA(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnectGrid1d);
    fnLA2dA(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d, MESH1D & Grid1d, CONNECT1D & ConnectGrid1d);
    void setCoeff(const Teuchos::RCP<COEFF> & Coeff);
    void setCoeff(COEFF & Coeff);
    void eval(const point3d & Yd, sVect<komplex> & mat);
    void eval(const point3d & Yd, sVect<Real> & mat);
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TEST, typename GEOSHAPE2D, typename COEFF>
fnLA2dA<TEST,GEOSHAPE2D,COEFF>::
fnLA2dA() : functionalLA2d<TEST,GEOSHAPE2D>()
{
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 1>::returnValue);
  
  assert(staticAssert<COEFF_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::geoName       == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typePoint3d>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename GEOSHAPE2D, typename COEFF>
fnLA2dA<TEST,GEOSHAPE2D,COEFF>::
fnLA2dA(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnectGrid1d) :
functionalLA2d<TEST,GEOSHAPE2D>(CommDev,Grid2d,ConnectGrid2d,Grid1d,ConnectGrid1d)
{
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 1>::returnValue);
  
  assert(staticAssert<COEFF_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::geoName       == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typePoint3d>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename GEOSHAPE2D, typename COEFF>
fnLA2dA<TEST,GEOSHAPE2D,COEFF>::
fnLA2dA(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d, MESH1D & Grid1d, CONNECT1D & ConnectGrid1d) :
functionalLA2d<TEST,GEOSHAPE2D>(CommDev,Grid2d,ConnectGrid2d,Grid1d,ConnectGrid1d)
{
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 1>::returnValue);
  
  assert(staticAssert<COEFF_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::geoName       == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typePoint3d>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename GEOSHAPE2D, typename COEFF>
void
fnLA2dA<TEST,GEOSHAPE2D,COEFF>::
setCoeff(const Teuchos::RCP<COEFF> & Coeff)
{
  coeffOk = true;
  coeff   = Coeff;
}

template<typename TEST, typename GEOSHAPE2D, typename COEFF>
void
fnLA2dA<TEST,GEOSHAPE2D,COEFF>::
setCoeff(COEFF & Coeff)
{
  coeffOk = true;
  coeff   = Teuchos::rcpFromRef(Coeff);
}

template<typename TEST, typename GEOSHAPE2D, typename COEFF>
void
fnLA2dA<TEST,GEOSHAPE2D,COEFF>::
eval(const point3d & Yd, sVect<Real> & mat)
{
  typedef functionalLA2d<TEST,GEOSHAPE2D> FUNCTIONAL;
  
  assert(coeffOk);
  assert(mat.size() == FUNCTIONAL::numIndex_test() );
  
  //Field eval
  sVect<TEST_OUTTYPE> eval_field(FUNCTIONAL::numIndex_test());
  FUNCTIONAL::eval_test(Yd, eval_field);
  
  //Normal evaluation
  point3d N = FUNCTIONAL::computeNormal(Yd);
  
  //Coeff eval
  point3d val;
  coeff->evalL(FUNCTIONAL::getEl1d(), FUNCTIONAL::mapSurfaceY(Yd), val);
  
  //Vectoring
  for(UInt j=1; j <= FUNCTIONAL::numIndex_test(); ++j)
  {
    mat(j) = point3d::dot(val, N) * eval_field(j) ;
  }
}

template<typename TEST, typename GEOSHAPE2D, typename COEFF>
void
fnLA2dA<TEST,GEOSHAPE2D,COEFF>::
eval(const point3d & Yd, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Yd,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}


#endif
