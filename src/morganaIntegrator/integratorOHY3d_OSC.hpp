/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTEGRATOROHY3D_OSC_HPP
#define INTEGRATOROHY3D_OSC_HPP

#include "traitsIntegratorOHY3d_OSC.hpp"


//_________________________________________________________________________________________________
// INTERFACE
//-------------------------------------------------------------------------------------------------
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
class integratorOHY3d_OSC : public traitsIntegratorOHY3d_OSC<OPERATOR,typename OPERATOR::GEOSHAPE2D,TYPE,PRECISION>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename OPERATOR::GEOSHAPE2D                                 GEOSHAPE2D;
    typedef traitsIntegratorOHY3d_OSC<OPERATOR,GEOSHAPE2D,TYPE,PRECISION> OSCTRAIT;
    typedef typename OSCTRAIT::INTCARD                                    INTCARD;
    //@}
  
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    integratorOHY3d_OSC();
    
    /*! Constructor */
    integratorOHY3d_OSC(const INTCARD & IntCard);
    //@}
};


template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
integratorOHY3d_OSC<OPERATOR,TYPE,PRECISION>::
integratorOHY3d_OSC() : traitsIntegratorOHY3d_OSC<OPERATOR,GEOSHAPE2D,TYPE,PRECISION>()
{
}

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
integratorOHY3d_OSC<OPERATOR,TYPE,PRECISION>::
integratorOHY3d_OSC(const INTCARD & IntCard) : traitsIntegratorOHY3d_OSC<OPERATOR,GEOSHAPE2D,TYPE,PRECISION>(IntCard)
{
}

#endif
