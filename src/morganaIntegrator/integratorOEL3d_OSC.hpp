/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTEGRATOROEL3D_OSC_HPP
#define INTEGRATOROEL3D_OSC_HPP

#include "traitsIntegratorOEL3d_OSC.hpp"


//_________________________________________________________________________________________________
// INTERFACE
//-------------------------------------------------------------------------------------------------
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
class integratorOEL3d_OSC : public traitsIntegratorOEL3d_OSC<OPERATOR,typename OPERATOR::FIELD_GEOSHAPE,TYPE,PRECISION>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename OPERATOR::FIELD_GEOSHAPE                           GEOSHAPE;
    typedef traitsIntegratorOEL3d_OSC<OPERATOR,GEOSHAPE,TYPE,PRECISION> OSCTRAIT;
    typedef typename OSCTRAIT::INTCARD                                  INTCARD;
    //@}
  
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    integratorOEL3d_OSC();
    
    /*! Constructor */
    integratorOEL3d_OSC(const INTCARD & IntCard);
    //@}
};


template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
integratorOEL3d_OSC<OPERATOR,TYPE,PRECISION>::
integratorOEL3d_OSC() : traitsIntegratorOEL3d_OSC<OPERATOR,GEOSHAPE,TYPE,PRECISION>()
{
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
integratorOEL3d_OSC<OPERATOR,TYPE,PRECISION>::
integratorOEL3d_OSC(const INTCARD & IntCard) : traitsIntegratorOEL3d_OSC<OPERATOR,GEOSHAPE,TYPE,PRECISION>(IntCard)
{
}

#endif
