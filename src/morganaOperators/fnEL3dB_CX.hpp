/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FNEL3DB_CX_HPP
#define FNEL3DB_CX_HPP

#include "morganaFiniteElements.hpp"
#include "functionalEL3d.hpp"
#include "../morganaDofs/komplex.h"


/*! Functional 3d: f * v complex test function complex source term */
template<typename TEST, typename COEFF>
class fnEL3dB_CX : public functionalEL3d<TEST>
{
    /*! @name Typedefs */ //@{
  public:
    typedef functionalEL3d<TEST>  FUNCTIONAL;
    
    typedef typename FUNCTIONAL::TEST_OUTTYPE  TEST_OUTTYPE;
    typedef typename FUNCTIONAL::MESH3D        MESH3D;
    typedef typename FUNCTIONAL::CONNECT3D     CONNECT3D;
    
    typedef typename TEST::GEOSHAPE      TEST_GEOSHAPE;
    typedef typename TEST::FIELD_DOFTYPE TEST_DOFTYPE;
    
    typedef typename COEFF::GEOSHAPE COEFF_GEOSHAPE;
    typedef typename COEFF::OUTTYPE  COEFF_OUTTYPE;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    bool coeffOk;
    Teuchos::RCP<COEFF> coeff;
    //@}
    
    /*! @name Functions */ //@{
  public:
    fnEL3dB_CX();
    fnEL3dB_CX(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    fnEL3dB_CX(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
    void setCoeff(const Teuchos::RCP<COEFF> & Coeff);
    void setCoeff(COEFF & Coeff);
    void eval(const point3d & Y, sVect<komplex> & mat);
    void eval(const point3d & Y, sVect<Real> & mat);
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TEST, typename COEFF>
fnEL3dB_CX<TEST,COEFF>::
fnEL3dB_CX() : functionalEL3d<TEST>()
{
  assert(staticAssert<COEFF_GEOSHAPE::geoName == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeKomplex>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
fnEL3dB_CX<TEST,COEFF>::
fnEL3dB_CX(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d) : functionalEL3d<TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<COEFF_GEOSHAPE::geoName == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeKomplex>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
fnEL3dB_CX<TEST,COEFF>::
fnEL3dB_CX(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d) : functionalEL3d<TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<COEFF_GEOSHAPE::geoName == TEST_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeKomplex>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
void
fnEL3dB_CX<TEST,COEFF>::
setCoeff(const Teuchos::RCP<COEFF> & Coeff)
{
  coeffOk = true;
  coeff   = Coeff;
}

template<typename TEST, typename COEFF>
void
fnEL3dB_CX<TEST,COEFF>::
setCoeff(COEFF & Coeff)
{
  coeffOk = true;
  coeff   = Teuchos::rcpFromRef(Coeff);
}

template<typename TEST, typename COEFF>
void
fnEL3dB_CX<TEST,COEFF>::
eval(const point3d & Y, sVect<komplex> & mat)
{
  //Assert
  assert(coeffOk);
  assert(mat.size() == FUNCTIONAL::numIndex_test() );
  
  //Typedefs
  typedef functionalEL3d<TEST>            FUNCTIONAL;
  typedef traitsMultiply<komplex,komplex> MULTIPLIER;
  
  //Field eval
  sVect<TEST_OUTTYPE> eval_field(FUNCTIONAL::numIndex_test());
  FUNCTIONAL::eval_test(Y,eval_field);
  
  //Coeff eval
  komplex val;
  coeff->evalL(FUNCTIONAL::el,Y,val);
  
  //Vectoring
  komplex V(Real(FUNCTIONAL::I_test == 1), -Real(FUNCTIONAL::I_test == 2));
  
  for(UInt j=1; j <= FUNCTIONAL::numIndex_test(); ++j)
  {
    mat(j) = MULTIPLIER::multiply(eval_field(j), V * val);
  }
}

template<typename TEST, typename COEFF>
void
fnEL3dB_CX<TEST,COEFF>::
eval(const point3d & Y, sVect<Real> & mat)
{
  sVect<komplex> Cmat(mat.size());
  eval(Y,Cmat);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex::cComp(Cmat(i),FUNCTIONAL::I_test); }
}

#endif
