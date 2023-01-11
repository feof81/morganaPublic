/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FNBF3DA_CX_HPP
#define FNBF3DA_CX_HPP

#include "../morganaFiniteElements/functionalBF3d.hpp"
#include "../morganaDofs/morganaTypes.hpp"
#include "../morganaDofs/point3d.h"
#include "../morganaDofs/komplex.h"


/*! Functional 3d: - v * coeff */
template<typename TEST, typename COEFF>
class fnBF3dA_CX : public functionalBF3d<TEST>
{
    /*! @name Typedefs */ //@{
  public:
    typedef functionalBF3d<TEST>  FUNCTIONAL;
    
    typedef typename FUNCTIONAL::TEST_OUTTYPE  TEST_OUTTYPE;
    typedef typename FUNCTIONAL::MESH2D        MESH2D;
    typedef typename FUNCTIONAL::MESH3D        MESH3D;
    typedef typename FUNCTIONAL::CONNECT2D     CONNECT2D;
    typedef typename FUNCTIONAL::CONNECT3D     CONNECT3D;
    
    typedef typename TEST::PMAPTYPE      TEST_PMAPTYPE;
    typedef typename TEST::GEOSHAPE      TEST_GEOSHAPE;
    typedef typename TEST::FIELD_DOFTYPE TEST_DOFTYPE;
    
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
    fnBF3dA_CX();
    fnBF3dA_CX(const Teuchos::RCP<communicator> & CommDev,
	       const Teuchos::RCP<MESH3D>       & Grid3d,
	       const Teuchos::RCP<CONNECT3D>    & ConnectGrid3d,
	       const Teuchos::RCP<MESH2D>       & Grid2d,
	       const Teuchos::RCP<CONNECT2D>    & ConnectGrid2d);
    
    fnBF3dA_CX(communicator & CommDev,
	       MESH3D       & Grid3d,
	       CONNECT3D    & ConnectGrid3d,
	       MESH2D       & Grid2d,
	       CONNECT2D    & ConnectGrid2d);
    
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
fnBF3dA_CX<TEST,COEFF>::
fnBF3dA_CX() : functionalBF3d<TEST>()
{
  typedef typename TEST_GEOSHAPE::GEOBSHAPE GEOPROOF;
  
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 2>::returnValue);
  
  assert(staticAssert<COEFF_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  assert(staticAssert<GEOPROOF::geoName            == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
fnBF3dA_CX<TEST,COEFF>::
fnBF3dA_CX(const Teuchos::RCP<communicator> & CommDev,
	   const Teuchos::RCP<MESH3D>       & Grid3d,
	   const Teuchos::RCP<CONNECT3D>    & ConnectGrid3d,
	   const Teuchos::RCP<MESH2D>       & Grid2d,
	   const Teuchos::RCP<CONNECT2D>    & ConnectGrid2d) :
functionalBF3d<TEST>(CommDev,Grid3d,ConnectGrid3d,Grid2d,ConnectGrid2d)
{
  typedef typename TEST_GEOSHAPE::GEOBSHAPE GEOPROOF;
  
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 2>::returnValue);
  
  assert(staticAssert<COEFF_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  assert(staticAssert<GEOPROOF::geoName            == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
fnBF3dA_CX<TEST,COEFF>::
fnBF3dA_CX(communicator & CommDev,
	   MESH3D       & Grid3d,
	   CONNECT3D    & ConnectGrid3d,
	   MESH2D       & Grid2d,
	   CONNECT2D    & ConnectGrid2d) : 
functionalBF3d<TEST>(CommDev,Grid3d,ConnectGrid3d,Grid2d,ConnectGrid2d)
{
  typedef typename TEST_GEOSHAPE::GEOBSHAPE GEOPROOF;
  
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 2>::returnValue);
  
  assert(staticAssert<COEFF_PMAPTYPE::parallelType == TEST_PMAPTYPE::parallelType>::returnValue);
  assert(staticAssert<GEOPROOF::geoName            == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
void
fnBF3dA_CX<TEST,COEFF>::
setCoeff(const Teuchos::RCP<COEFF> & Coeff)
{
  coeffOk = true;
  coeff   = Coeff;
}

template<typename TEST, typename COEFF>
void
fnBF3dA_CX<TEST,COEFF>::
setCoeff(COEFF & Coeff)
{
  coeffOk = true;
  coeff   = Teuchos::rcpFromRef(Coeff);
}

template<typename TEST, typename COEFF>
void
fnBF3dA_CX<TEST,COEFF>::
eval(const point3d & Yf, sVect<komplex> & mat)
{
  typedef functionalBF3d<TEST> FUNCTIONAL;
  
  assert(coeffOk);
  assert(mat.size() == FUNCTIONAL::numIndex_test() );
  
  //Field eval
  sVect<TEST_OUTTYPE> eval_field(FUNCTIONAL::numIndex_test());
  FUNCTIONAL::eval_test(Yf, eval_field);  
  
  //Coeff eval
  Real val;
  coeff->evalL(FUNCTIONAL::getEl2d(), FUNCTIONAL::mapSurfaceY(Yf), val);
  komplex V(Real(FUNCTIONAL::I_test == 1), -Real(FUNCTIONAL::I_test == 2));
  
  //Vectoring
  for(UInt j=1; j <= FUNCTIONAL::numIndex_test(); ++j)
  {
    mat(j) = V * eval_field(j) * (val * (-1.0));
  }
}

template<typename TEST, typename COEFF>
void
fnBF3dA_CX<TEST,COEFF>::
eval(const point3d & Yf, sVect<Real> & mat)
{
  sVect<komplex> matC(mat.size());
  eval(Yf,matC);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex::cComp(matC(i),FUNCTIONAL::I_test); }
}

#endif
