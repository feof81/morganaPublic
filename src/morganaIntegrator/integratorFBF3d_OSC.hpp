/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTEGRATORFBF3D_OSC_HPP
#define INTEGRATORFBF3D_OSC_HPP

#include "traitsIntegratorFBF3d_OSC.hpp"


//_________________________________________________________________________________________________
// INTERFACE
//-------------------------------------------------------------------------------------------------
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
class integratorFBF3d_OSC : public traitsIntegratorFBF3d_OSC<FUNCTIONAL,typename FUNCTIONAL::GEOSHAPE2D,TYPE,PRECISION>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename FUNCTIONAL::GEOSHAPE2D                               GEOSHAPE;
    typedef traitsIntegratorFBF3d_OSC<FUNCTIONAL,GEOSHAPE,TYPE,PRECISION> OSCTRAIT;
    typedef typename OSCTRAIT::INTCARD                                    INTCARD;
    //@}
  
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    integratorFBF3d_OSC();
    
    /*! Constructor */
    integratorFBF3d_OSC(const INTCARD & IntCard);
    //@}
};


template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
integratorFBF3d_OSC<FUNCTIONAL,TYPE,PRECISION>::
integratorFBF3d_OSC() : traitsIntegratorFBF3d_OSC<FUNCTIONAL,GEOSHAPE,TYPE,PRECISION>()
{
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
integratorFBF3d_OSC<FUNCTIONAL,TYPE,PRECISION>::
integratorFBF3d_OSC(const INTCARD & IntCard) : traitsIntegratorFBF3d_OSC<FUNCTIONAL,GEOSHAPE,TYPE,PRECISION>(IntCard)
{
}

#endif
