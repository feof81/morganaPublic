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

#include "morganaFiniteElements.hpp"
#include "functionalEL2d.hpp"
#include "../morganaDofs/komplex.h"


template<typename TEST>
class fnEL2dC : public functionalEL2d<TEST>
{
  /*! @name Typedefs */ //@{
  public:
    typedef functionalEL2d<TEST>  FUNCTIONAL;
    
    typedef typename FUNCTIONAL::TEST_OUTTYPE  TEST_OUTTYPE;
    typedef typename FUNCTIONAL::MESH2D        MESH2D;
    typedef typename FUNCTIONAL::CONNECT2D     CONNECT2D;
    
    typedef typename TEST::GEOSHAPE   TEST_GEOSHAPE;
    //@}
    
    /*! @name Functions */ //@{
  public:
    fnEL2dC();
    fnEL2dC(const Teuchos::RCP<communicator> & CommDev,
            const Teuchos::RCP<MESH2D>       & Grid2d,
            const Teuchos::RCP<CONNECT2D>    & ConnedGrid2d);
    
    fnEL2dC(communicator & CommDev,
            MESH2D       & Grid2d,
            CONNECT2D    & ConnedGrid2d);
    
    void eval(const point3d & Y, sVect<komplex> & mat);
    void eval(const point3d & Y, sVect<Real> & mat);
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TEST>
fnEL2dC<TEST>::
fnEL2dC() : functionalEL2d<TEST>()
{
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename TEST>
fnEL2dC<TEST>::
fnEL2dC(const Teuchos::RCP<communicator> & CommDev,
        const Teuchos::RCP<MESH2D>       & Grid2d,
        const Teuchos::RCP<CONNECT2D>    & ConnedGrid2d) : functionalEL2d<TEST>(CommDev,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename TEST>
fnEL2dC<TEST>::
fnEL2dC(communicator & CommDev,
        MESH2D       & Grid2d,
        CONNECT2D    & ConnedGrid2d) : functionalEL2d<TEST>(CommDev,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename TEST>
void
fnEL2dC<TEST>::
eval(const point3d & Y, sVect<Real> & mat)
{
  typedef functionalEL2d<TEST> FUNCTIONAL;
  
  assert(mat.size() == FUNCTIONAL::numIndex_test() );
  
  //Field eval
  sVect<TEST_OUTTYPE> eval_field(FUNCTIONAL::numIndex_test());
  FUNCTIONAL::eval_test(Y,eval_field);
  
  //Vectoring
  for(UInt j=1; j <= FUNCTIONAL::numIndex_test(); ++j)
  { mat(j) = eval_field(j); }
}

template<typename TEST>
void
fnEL2dC<TEST>::
eval(const point3d & Y, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Y,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}

#endif
