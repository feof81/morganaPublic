/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FNHY3DA_CX_HPP
#define FNHY3DA_CX_HPP

#include "morganaFiniteElements.hpp"
#include "functionalHY3d.hpp"
#include "../morganaDofs/komplex.h"


/*! Functional hybrid 3d, complex version: v * c*/
template<typename TEST, typename COEFF>
class fnHY3dA_CX : public functionalHY3d<TEST>
{
  /*! @name Typedefs */ //@{
  public:
    typedef functionalHY3d<TEST>  FUNCTIONAL;
    
    typedef typename FUNCTIONAL::TEST_OUTTYPE  TEST_OUTTYPE;
    typedef typename FUNCTIONAL::MESH3D        MESH3D;
    typedef typename FUNCTIONAL::CONNECT3D     CONNECT3D;
    
    typedef typename TEST::GEOSHAPE      TEST_GEOSHAPE;
    typedef typename TEST::FIELD_DOFTYPE TEST_DOFTYPE;
    typedef typename COEFF::GEOSHAPE  COEFF_GEOSHAPE;
    typedef typename COEFF::OUTTYPE   COEFF_OUTTYPE;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    bool coeffOk;
    Teuchos::RCP<COEFF> coeff;
    //@}
    
    /*! @name Functions */ //@{
  public:
    fnHY3dA_CX();
    fnHY3dA_CX(const Teuchos::RCP<communicator> & CommDev,
               const Teuchos::RCP<MESH3D>       & Grid3d,
               const Teuchos::RCP<CONNECT3D>    & ConnedGrid3d);
    
    fnHY3dA_CX(communicator & CommDev,
               MESH3D       & Grid3d,
               CONNECT3D    & ConnedGrid3d);
    
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
fnHY3dA_CX<TEST,COEFF>::
fnHY3dA_CX() : functionalHY3d<TEST>()
{
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
fnHY3dA_CX<TEST,COEFF>::
fnHY3dA_CX(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d) : functionalHY3d<TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
fnHY3dA_CX<TEST,COEFF>::
fnHY3dA_CX(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d) : functionalHY3d<TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<traitsBasic<TEST_DOFTYPE>::myType  == typeKomplex>::returnValue);
  assert(staticAssert<traitsBasic<COEFF_OUTTYPE>::myType == typeReal>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
void
fnHY3dA_CX<TEST,COEFF>::
setCoeff(const Teuchos::RCP<COEFF> & Coeff)
{
  coeffOk = true;
  coeff   = Coeff;
}

template<typename TEST, typename COEFF>
void
fnHY3dA_CX<TEST,COEFF>::
setCoeff(COEFF & Coeff)
{
  coeffOk = true;
  coeff   = Teuchos::rcpFromRef(Coeff);
}

template<typename TEST, typename COEFF>
void
fnHY3dA_CX<TEST,COEFF>::
eval(const point3d & Yf, sVect<Real> & mat)
{
  typedef functionalHY3d<TEST> FUNCTIONAL;
  
  assert(mat.size() == FUNCTIONAL::numIndex_test() );
  
  //Test eval
  sVect<TEST_OUTTYPE> eval_test(FUNCTIONAL::numIndex_test());
  FUNCTIONAL::eval_test(Yf,eval_test);
  
  //Normal and Yv eval
  point3d Yv = FUNCTIONAL::mapVolumeY(Yf);
  
  //Coeff eval
  Real val;
  coeff->evalL(FUNCTIONAL::el3d,Yv,val);
  komplex V(Real(FUNCTIONAL::I_test == 1), -Real(FUNCTIONAL::I_test == 2));
  
  //Vectoring
  for(UInt j=1; j <= FUNCTIONAL::numIndex_test(); ++j)
  { mat(j) = komplex::cComp(0.5 * eval_test(j) * V * (-val), FUNCTIONAL::I_test); }
}

template<typename TEST, typename COEFF>
void
fnHY3dA_CX<TEST,COEFF>::
eval(const point3d & Yf, sVect<komplex> & mat)
{
  typedef functionalHY3d<TEST> FUNCTIONAL;
  
  assert(mat.size() == FUNCTIONAL::numIndex_test() );
  
  //Test eval
  sVect<TEST_OUTTYPE> eval_test(FUNCTIONAL::numIndex_test());
  FUNCTIONAL::eval_test(Yf,eval_test);
  
  //Normal and Yv eval
  point3d Yv = FUNCTIONAL::mapVolumeY(Yf);
  
  //Coeff eval
  Real val;
  coeff->evalL(FUNCTIONAL::el3d,Yv,val);
  komplex V(Real(FUNCTIONAL::I_test == 1), -Real(FUNCTIONAL::I_test == 2));
  
  //Vectoring
  for(UInt j=1; j <= FUNCTIONAL::numIndex_test(); ++j)
  { mat(j) = 0.5 * eval_test(j) * V * (-val); }
}

#endif
