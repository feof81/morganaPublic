/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FNEL2DC_HPP
#define FNEL2DC_HPP

#include "traitsMultiply.h"
#include "morganaFiniteElements.hpp"
#include "functionalEL2d.hpp"
#include "../morganaDofs/komplex.h"


/*! Functional 2d f * v (where * is the scalar product between f and v, f and v has the same outtype) */
template<typename TEST, typename COEFF>
class fnEL2dD : public functionalEL2d<TEST>
{
  /*! @name Typedefs */ //@{
  public:
    typedef functionalEL2d<TEST>  FUNCTIONAL;
    
    typedef typename FUNCTIONAL::TEST_OUTTYPE  TEST_OUTTYPE;
    typedef typename FUNCTIONAL::MESH2D        MESH2D;
    typedef typename FUNCTIONAL::CONNECT2D     CONNECT2D;
    
    typedef typename TEST::GEOSHAPE   TEST_GEOSHAPE;
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
    fnEL2dD();
    fnEL2dD(const Teuchos::RCP<communicator> & CommDev,
            const Teuchos::RCP<MESH2D>       & Grid2d,
            const Teuchos::RCP<CONNECT2D>    & ConnedGrid2d);
    
    fnEL2dD(communicator & CommDev,
            MESH2D       & Grid2d,
            CONNECT2D    & ConnedGrid2d);
    
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
fnEL2dD<TEST,COEFF>::
fnEL2dD() : functionalEL2d<TEST>()
{
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType == traitsBasic<COEFF_OUTTYPE>::myType>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
fnEL2dD<TEST,COEFF>::
fnEL2dD(const Teuchos::RCP<communicator> & CommDev,
        const Teuchos::RCP<MESH2D>       & Grid2d,
        const Teuchos::RCP<CONNECT2D>    & ConnedGrid2d) : functionalEL2d<TEST>(CommDev,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType == traitsBasic<COEFF_OUTTYPE>::myType>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
fnEL2dD<TEST,COEFF>::
fnEL2dD(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d) : functionalEL2d<TEST>(CommDev,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType == traitsBasic<COEFF_OUTTYPE>::myType>::returnValue);
  
  coeffOk = false;
}

template<typename TEST, typename COEFF>
void
fnEL2dD<TEST,COEFF>::
setCoeff(const Teuchos::RCP<COEFF> & Coeff)
{
  coeffOk = true;
  coeff   = Coeff;
}

template<typename TEST, typename COEFF>
void
fnEL2dD<TEST,COEFF>::
setCoeff(COEFF & Coeff)
{
  coeffOk = true;
  coeff   = Teuchos::rcpFromRef(Coeff);
}

template<typename TEST, typename COEFF>
void
fnEL2dD<TEST,COEFF>::
eval(const point3d & Y, sVect<Real> & mat)
{
  typedef functionalEL2d<TEST> FUNCTIONAL;
  
  assert(coeffOk);
  assert(mat.size() == FUNCTIONAL::numIndex_test() );
  
  //Field eval
  sVect<TEST_OUTTYPE> eval_test(FUNCTIONAL::numIndex_test());
  FUNCTIONAL::eval_test(Y,eval_test);
  
  //Coeff eval
  COEFF_OUTTYPE val;
  coeff->evalL(FUNCTIONAL::el,Y,val);

  //Vectoring
  for(UInt j=1; j <= FUNCTIONAL::numIndex_test(); ++j)
  { mat(j) = traitsMultiply<TEST_OUTTYPE,COEFF_OUTTYPE>::scalarProductA(eval_test(j) , val); }
}

template<typename TEST, typename COEFF>
void
fnEL2dD<TEST,COEFF>::
eval(const point3d & Y, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Y,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}


#endif
