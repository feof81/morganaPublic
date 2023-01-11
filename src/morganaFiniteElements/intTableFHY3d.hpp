/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef INTTABLEFHY3D_HPP
#define INTTABLEFHY3D_HPP

#include "feHY3d.hpp"
#include "feHY3d_extern.hpp"
#include "feHYB3d.hpp"
#include "feHYB3d_extern.hpp"
#include "feOscHYB3d.hpp"
#include "feOscHYB3d_extern.hpp"

#include "morganaFields.hpp"
#include "morganaIntegrator.hpp"

template<FEBaseLabel > struct intTableFHY3d;  //! Generic class

template<> struct intTableFHY3d< BL_HY_3d >      { static const intClass intFlag = FHY3d_STD; };   //! HY
template<> struct intTableFHY3d< BL_HYB_3d >     { static const intClass intFlag = FHY3d_STD; };   //! HYB
template<> struct intTableFHY3d< BL_OSCHYB_3d >  { static const intClass intFlag = FHY3d_OSC; };   //! OSC-HYB

#endif
