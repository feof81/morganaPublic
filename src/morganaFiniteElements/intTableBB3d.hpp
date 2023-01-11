/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef INTTABLEBB3D_HPP
#define INTTABLEBB3D_HPP

#include "fePr2d.hpp"
#include "feQr2d.hpp"
#include "fePr3d.hpp"
#include "feOscPr3d.hpp"
#include "feQr3d.hpp"
#include "morganaFields.hpp"
#include "morganaIntegrator.hpp"


template<FEBaseLabel, FEBaseLabel > struct intTableBB3d;  //! Generic class

template<> struct intTableBB3d< BL_Pr_3d,     BL_Pr_3d >     { static const intClass intFlag = OBB3d_STD; };   //! Pr3d - Pr3d
template<> struct intTableBB3d< BL_RT0LT_3d,  BL_RT0LT_3d >  { static const intClass intFlag = OBB3d_STD; };   //! RT0  - RT0
template<> struct intTableBB3d< BL_RT0LOC_3d, BL_RT0LOC_3d > { static const intClass intFlag = OBB3d_STD; };   //! RT0-loc  - RT0-loc
template<> struct intTableBB3d< BL_OSwPr_3d,  BL_OSwPr_3d >  { static const intClass intFlag = OBB3d_OSC; };   //! OscPr3d - OscPr3d
template<> struct intTableBB3d< BL_OSCRr_3d,  BL_OSCRr_3d >  { static const intClass intFlag = OBB3d_OSC; };   //! OscCr3d - OscCr3d

#endif
