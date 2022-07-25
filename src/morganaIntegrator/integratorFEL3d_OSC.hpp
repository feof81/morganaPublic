/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTEGRATORFEL3D_OSC_HPP
#define INTEGRATORFEL3D_OSC_HPP

#include "traitsIntegratorFEL3d_OSC.hpp"


//_________________________________________________________________________________________________
// INTERFACE
//-------------------------------------------------------------------------------------------------
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
class integratorFEL3d_OSC : public traitsIntegratorFEL3d_OSC<FUNCTIONAL,typename FUNCTIONAL::TEST_GEOSHAPE,TYPE,PRECISION>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename FUNCTIONAL::TEST_GEOSHAPE                            GEOSHAPE;
    typedef traitsIntegratorFEL3d_OSC<FUNCTIONAL,GEOSHAPE,TYPE,PRECISION> OSCTRAIT;
    typedef typename OSCTRAIT::INTCARD                                     INTCARD;
    //@}
  
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    integratorFEL3d_OSC();
    
    /*! Constructor */
    integratorFEL3d_OSC(const INTCARD & IntCard);
    //@}
};


template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
integratorFEL3d_OSC<FUNCTIONAL,TYPE,PRECISION>::
integratorFEL3d_OSC() : traitsIntegratorFEL3d_OSC<FUNCTIONAL,GEOSHAPE,TYPE,PRECISION>()
{
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
integratorFEL3d_OSC<FUNCTIONAL,TYPE,PRECISION>::
integratorFEL3d_OSC(const INTCARD & IntCard) : traitsIntegratorFEL3d_OSC<FUNCTIONAL,GEOSHAPE,TYPE,PRECISION>(IntCard)
{
}

#endif
