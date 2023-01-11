/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FNBF3DD_HPP
#define FNBF3DD_HPP
#include "../morganaFiniteElements/functionalBF3d.hpp"
#include "../morganaDofs/morganaTypes.hpp"
#include "../morganaDofs/point3d.h"
#include "../morganaDofs/komplex.h"


/*! Functional 3d: (Coeff * N) v */
template<typename TEST, typename COEFF>
class fnBF3dD : public functionalBF3d<TEST>
{
  /*! @name Typedefs */ //@{
  public:
    typedef functionalBF3d<TEST>  FUNCTIONAL;
    
    typedef typename FUNCTIONAL::TEST_OUTTYPE  TEST_OUTTYPE;
    typedef typename FUNCTIONAL::MESH2D        MESH2D;
    typedef typename FUNCTIONAL::MESH3D        MESH3D;
    typedef typename FUNCTIONAL::CONNECT2D     CONNECT2D;
    typedef typename FUNCTIONAL::CONNECT3D     CONNECT3D;
    
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
    fnBF3dD();
    fnBF3dD(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    fnBF3dD(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnectGrid3d, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d);
    void setCoeff(const Teuchos::RCP<COEFF> & Coeff);
    void setCoeff(COEFF & Coeff);
    void eval(const point3d & Yf, sVect<komplex> & mat);
    void eval(const point3d & Yf, sVect<Real> & mat);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TEST, typename COEFF>
fnBF3dD<TEST,COEFF>::
fnBF3dD()
{
  typedef typename TEST_GEOSHAPE::GEOBSHAPE TEST_GEOBSHAPE;
  
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 2>::returnValue);
  
  assert(staticAssert<TEST_GEOBSHAPE::geoName == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<COEFF_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typePoint3d>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
fnBF3dD<TEST,COEFF>::
fnBF3dD(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d) :
functionalBF3d<TEST>(CommDev,Grid3d,ConnectGrid3d,Grid2d,ConnectGrid2d)
{
  typedef typename TEST_GEOSHAPE::GEOBSHAPE TEST_GEOBSHAPE;
  
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 2>::returnValue);
  
  assert(staticAssert<TEST_GEOBSHAPE::geoName == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<COEFF_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typePoint3d>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
fnBF3dD<TEST,COEFF>::
fnBF3dD(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnectGrid3d, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d) : 
functionalBF3d<TEST>(CommDev,Grid3d,ConnectGrid3d,Grid2d,ConnectGrid2d)
{
  typedef typename TEST_GEOSHAPE::GEOBSHAPE TEST_GEOBSHAPE;
  
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 2>::returnValue);
  
  assert(staticAssert<TEST_GEOBSHAPE::geoName == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<COEFF_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typePoint3d>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
void
fnBF3dD<TEST,COEFF>::
setCoeff(const Teuchos::RCP<COEFF> & Coeff)
{
  coeffOk = true;
  coeff   = Coeff;
}

template<typename TEST, typename COEFF>
void
fnBF3dD<TEST,COEFF>::
setCoeff(COEFF & Coeff)
{
  coeffOk = true;
  coeff   = Teuchos::rcpFromRef(Coeff);
}

template<typename TEST, typename COEFF>
void
fnBF3dD<TEST,COEFF>::
eval(const point3d & Yf, sVect<Real> & mat)
{
  typedef functionalBF3d<TEST> FUNCTIONAL;
  
  assert(coeffOk);
  assert(mat.size() == FUNCTIONAL::numIndex_test() );
  
  //Normal evaluation
  point3d N = FUNCTIONAL::computeNormal(Yf);
  
  //Field eval
  sVect<TEST_OUTTYPE> eval_field(FUNCTIONAL::numIndex_test());
  FUNCTIONAL::eval_test(Yf, eval_field);
  
  //Coeff
  point3d F;
  coeff->evalL(FUNCTIONAL::getEl2d(), Yf, F);
  
  //Vectoring
  for(UInt j=1; j <= FUNCTIONAL::numIndex_test(); ++j)
  { mat(j) = point3d::dot(F, N) * eval_field(j); }
}

template<typename TEST, typename COEFF>
void
fnBF3dD<TEST,COEFF>::
eval(const point3d & Yf, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Yf,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}

#endif
