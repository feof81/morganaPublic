/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef INTTABLEFLA3D_HPP
#define INTTABLEFLA3D_HPP

#include "fePr2d.hpp"
#include "feOscPr2d.hpp"
#include "feQr2d.hpp"
#include "feRt0LT2d.hpp"
#include "morganaFields.hpp"
#include "morganaIntegrator.hpp"


template<FEBaseLabel> struct intTableFLA3d;  //! Generic class

template<> struct intTableFLA3d< BL_Pr_2d    >  { static const intClass intFlag = FLA3d_STD; };   //! Pr
template<> struct intTableFLA3d< BL_OSwPr_2d >  { static const intClass intFlag = FLA3d_OSC; };   //! OSC_Pr3d
template<> struct intTableFLA3d< BL_Qr_2d    >  { static const intClass intFlag = FLA3d_STD; };   //! Qr
template<> struct intTableFLA3d< BL_RT0LT_2d >  { static const intClass intFlag = FLA3d_STD; };   //! RT0 


#endif
