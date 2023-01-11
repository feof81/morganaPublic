/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FNEL3DB_HPP
#define FNEL3DB_HPP

#include "morganaFiniteElements.hpp"
#include "functionalEL3d.hpp"
#include "../morganaDofs/komplex.h"


/*! Functional 3d:  v */
template<typename TEST>
class fnEL3dB : public functionalEL3d<TEST>
{
  /*! @name Typedefs */ //@{
  public:
    typedef functionalEL3d<TEST>  FUNCTIONAL;
    
    typedef typename FUNCTIONAL::TEST_OUTTYPE  TEST_OUTTYPE;
    typedef typename FUNCTIONAL::MESH3D        MESH3D;
    typedef typename FUNCTIONAL::CONNECT3D     CONNECT3D;
    //@}
    
    /*! @name Functions */ //@{
  public:
    fnEL3dB();
    fnEL3dB(const Teuchos::RCP<communicator> & CommDev,
            const Teuchos::RCP<MESH3D>       & Grid3d,
            const Teuchos::RCP<CONNECT3D>    & ConnedGrid3d);
    
    fnEL3dB(communicator & CommDev,
            MESH3D       & Grid3d,
            CONNECT3D    & ConnedGrid3d);
    
    void eval(const point3d & Y, sVect<komplex> & mat);
    void eval(const point3d & Y, sVect<Real> & mat);
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TEST>
fnEL3dB<TEST>::
fnEL3dB() : functionalEL3d<TEST>()
{
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename TEST>
fnEL3dB<TEST>::
fnEL3dB(const Teuchos::RCP<communicator> & CommDev,
	const Teuchos::RCP<MESH3D>       & Grid3d,
	const Teuchos::RCP<CONNECT3D>    & ConnedGrid3d) : functionalEL3d<TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename TEST>
fnEL3dB<TEST>::
fnEL3dB(communicator & CommDev,
        MESH3D       & Grid3d,
        CONNECT3D    & ConnedGrid3d) : functionalEL3d<TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename TEST>
void
fnEL3dB<TEST>::
eval(const point3d & Y, sVect<Real> & mat)
{
  typedef functionalEL3d<TEST> FUNCTIONAL;
  
  assert(mat.size() == FUNCTIONAL::numIndex_test() );
  
  //Field eval
  sVect<TEST_OUTTYPE> eval_field(FUNCTIONAL::numIndex_test());
  FUNCTIONAL::eval_test(Y,eval_field);
  
  //Vectoring
  for(UInt j=1; j <= FUNCTIONAL::numIndex_test(); ++j)
  {
    mat(j) = eval_field(j);
  }
}

template<typename TEST>
void
fnEL3dB<TEST>::
eval(const point3d & Y, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Y,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}

#endif
