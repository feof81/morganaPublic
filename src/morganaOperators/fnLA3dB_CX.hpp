/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FNLA3DB_CX_HPP
#define FNLA3DB_CX_HPP

#include "../morganaFiniteElements/functionalLA3d.hpp"
#include "../morganaDofs/morganaTypes.hpp"
#include "../morganaDofs/point3d.h"
#include "../morganaDofs/komplex.h"


/*!  Functional: u * v  */
template<typename TEST, typename GEOSHAPE3D, typename COEFF>
class fnLA3dB_CX : public functionalLA3d<TEST,GEOSHAPE3D>
{
    /*! @name Typedefs */ //@{
  public:
    typedef functionalLA3d<TEST,GEOSHAPE3D>    FUNCTIONAL;
    
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
    fnLA3dB_CX();
    fnLA3dB_CX(const Teuchos::RCP<communicator> & CommDev,
               const Teuchos::RCP<MESH3D>       & Grid3d,
               const Teuchos::RCP<CONNECT3D>    & ConnectGrid3d,
               const Teuchos::RCP<MESH2D>       & Grid2d,
               const Teuchos::RCP<CONNECT2D>    & ConnectGrid2d);
    
    fnLA3dB_CX(communicator & CommDev,
               MESH3D       & Grid3d,
               CONNECT3D    & ConnectGrid3d,
               MESH2D       & Grid2d,
               CONNECT2D    & ConnectGrid2d);
    
    void setCoeff(const Teuchos::RCP<COEFF> & Coeff);
    void setCoeff(COEFF & Coeff);
    void eval(const point3d & Yf, sVect<komplex> & mat);
    void eval(const point3d & Yf, sVect<Real>    & mat);
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TEST, typename GEOSHAPE3D, typename COEFF>
fnLA3dB_CX<TEST,GEOSHAPE3D,COEFF>::
fnLA3dB_CX() : functionalLA3d<TEST,GEOSHAPE3D>()
{
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 2>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::geoName == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename GEOSHAPE3D, typename COEFF>
fnLA3dB_CX<TEST,GEOSHAPE3D,COEFF>::
fnLA3dB_CX(const Teuchos::RCP<communicator> & CommDev,
	   const Teuchos::RCP<MESH3D>       & Grid3d,
	   const Teuchos::RCP<CONNECT3D>    & ConnectGrid3d,
	   const Teuchos::RCP<MESH2D>       & Grid2d,
	   const Teuchos::RCP<CONNECT2D>    & ConnectGrid2d) :
functionalLA3d<TEST,GEOSHAPE3D>(CommDev,Grid3d,ConnectGrid3d,Grid2d,ConnectGrid2d)
{
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 2>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::geoName == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename GEOSHAPE3D, typename COEFF>
fnLA3dB_CX<TEST,GEOSHAPE3D,COEFF>::
fnLA3dB_CX(communicator & CommDev,
	   MESH3D       & Grid3d,
	   CONNECT3D    & ConnectGrid3d,
	   MESH2D       & Grid2d,
	   CONNECT2D    & ConnectGrid2d) :
functionalLA3d<TEST,GEOSHAPE3D>(CommDev,Grid3d,ConnectGrid3d,Grid2d,ConnectGrid2d)
{
  assert(staticAssert<COEFF_GEOSHAPE::nDim == 2>::returnValue);
  assert(staticAssert<TEST_GEOSHAPE::geoName == COEFF_GEOSHAPE::geoName>::returnValue);
  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename GEOSHAPE3D, typename COEFF>
void
fnLA3dB_CX<TEST,GEOSHAPE3D,COEFF>::
setCoeff(const Teuchos::RCP<COEFF> & Coeff)
{
  coeffOk = true;
  coeff   = Coeff;
}

template<typename TEST, typename GEOSHAPE3D, typename COEFF>
void
fnLA3dB_CX<TEST,GEOSHAPE3D,COEFF>::
setCoeff(COEFF & Coeff)
{
  coeffOk = true;
  coeff   = Teuchos::rcpFromRef(Coeff);
}

template<typename TEST, typename GEOSHAPE3D, typename COEFF>
void
fnLA3dB_CX<TEST,GEOSHAPE3D,COEFF>::
eval(const point3d & Yf, sVect<komplex> & mat)
{
  typedef functionalLA3d<TEST,GEOSHAPE3D> FUNCTIONAL;
  
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
    mat(j) = V * eval_field(j) * val;
  }
}

template<typename TEST, typename GEOSHAPE3D, typename COEFF>
void
fnLA3dB_CX<TEST,GEOSHAPE3D,COEFF>::
eval(const point3d & Yf, sVect<Real> & mat)
{
  sVect<komplex> Cmat(mat.size());
  eval(Yf,Cmat);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex::cComp(Cmat(i),FUNCTIONAL::I_test); }
}


#endif
