/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef INTTABLEHY3D_HPP
#define INTTABLEHY3D_HPP

#include "feHY3d.hpp"
#include "feHY3d_extern.hpp"
#include "feOscCr3d.hpp"
#include "feOscHYB3d.hpp"
#include "feRt0Loc3d.hpp"
#include "morganaFields.hpp"
#include "morganaIntegrator.hpp"

template<FEBaseLabel, FEBaseLabel > struct intTableHY3d;  //! Generic class

template<> struct intTableHY3d< BL_RT0LOC_3d, BL_HY_3d >     { static const intClass intFlag = OHY3d_STD; };  //! RT0_LOC - HY
template<> struct intTableHY3d< BL_OSCRr_3d,  BL_OSCHYB_3d > { static const intClass intFlag = OHY3d_OSC; };  //! OSC_CR  - OSC_HY

#endif
