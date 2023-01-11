/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FNBF2DA_HPP
#define FNBF2DA_HPP
#include "../morganaFiniteElements/functionalBF2d.hpp"
#include "../morganaDofs/morganaTypes.hpp"
#include "../morganaDofs/point3d.h"
#include "../morganaDofs/komplex.h"


template<typename TEST, typename COEFF>
class fnBF2dA : public functionalBF2d<TEST>
{
    /*! @name Typedefs */ //@{
  public:
    typedef functionalBF2d<TEST>  FUNCTIONAL;
    
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
    fnBF2dA();
    fnBF2dA(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnectGrid1d);
    fnBF2dA(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d, MESH1D & Grid1d, CONNECT1D & ConnectGrid1d);
    void setCoeff(const Teuchos::RCP<COEFF> & Coeff);
    void setCoeff(COEFF & Coeff);
    void eval(const point3d & Ys, sVect<komplex> & mat);
    void eval(const point3d & Ys, sVect<Real> & mat);
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TEST, typename COEFF>
fnBF2dA<TEST,COEFF>::
fnBF2dA() : functionalBF2d<TEST>()
{
  typedef typename TEST_GEOSHAPE::GEOBSHAPE GEOPROOF;
  
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 1>::returnValue);
  
  assert(staticAssert<COEFF_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  assert(staticAssert<GEOPROOF::geoName            == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
fnBF2dA<TEST,COEFF>::
fnBF2dA(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnectGrid1d) :
functionalBF2d<TEST>(CommDev,Grid2d,ConnectGrid2d,Grid1d,ConnectGrid1d)
{
  typedef typename TEST_GEOSHAPE::GEOBSHAPE GEOPROOF;
  
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 1>::returnValue);
  
  assert(staticAssert<COEFF_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  assert(staticAssert<GEOPROOF::geoName            == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
fnBF2dA<TEST,COEFF>::
fnBF2dA(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d, MESH1D & Grid1d, CONNECT1D & ConnectGrid1d) : 
functionalBF2d<TEST>(CommDev,Grid2d,ConnectGrid2d,Grid1d,ConnectGrid1d)
{
  typedef typename TEST_GEOSHAPE::GEOBSHAPE GEOPROOF;
  
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 1>::returnValue);
  
  assert(staticAssert<COEFF_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  assert(staticAssert<GEOPROOF::geoName            == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
void
fnBF2dA<TEST,COEFF>::
setCoeff(const Teuchos::RCP<COEFF> & Coeff)
{
  coeffOk = true;
  coeff   = Coeff;
}

template<typename TEST, typename COEFF>
void
fnBF2dA<TEST,COEFF>::
setCoeff(COEFF & Coeff)
{
  coeffOk = true;
  coeff   = Teuchos::rcpFromRef(Coeff);
}

template<typename TEST, typename COEFF>
void
fnBF2dA<TEST,COEFF>::
eval(const point3d & Ys, sVect<Real> & mat)
{
  typedef functionalBF2d<TEST> FUNCTIONAL;
  
  assert(coeffOk);
  assert(mat.size() == FUNCTIONAL::numIndex_test() );
  
  //Field eval
  sVect<TEST_OUTTYPE> eval_field(FUNCTIONAL::numIndex_test());
  FUNCTIONAL::eval_test(Ys, eval_field);
  
  //Normal evaluation
  point3d N = FUNCTIONAL::computeNormal(Ys);
  
  //Coeff eval
  Real val;
  coeff->evalL(FUNCTIONAL::getEl1d(), Ys, val);
  
  //Vectoring
  for(UInt j=1; j <= FUNCTIONAL::numIndex_test(); ++j)
  {
    mat(j) = - (point3d::dot(eval_field(j), N) * val);
  }
}

template<typename TEST, typename COEFF>
void
fnBF2dA<TEST,COEFF>::
eval(const point3d & Ys, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Ys,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}

#endif
