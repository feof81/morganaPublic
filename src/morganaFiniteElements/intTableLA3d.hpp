/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef INTTABLELA3D_HPP
#define INTTABLELA3D_HPP

#include "fePr2d.hpp"
#include "feOscPr2d.hpp"
#include "feQr2d.hpp"
#include "feRt0LT2d.hpp"
#include "fePr3d.hpp"
#include "feOscPr3d.hpp"
#include "feQr3d.hpp"
#include "feSpectralLH2d.hpp"
#include "feSpectralLH3d.hpp"
#include "morganaFields.hpp"
#include "morganaIntegrator.hpp"


template<FEBaseLabel, FEBaseLabel > struct intTableLA3d;  //! Generic class

template<> struct intTableLA3d< BL_Pr_3d,     BL_Pr_2d    >  { static const intClass intFlag = OLA3d_STD; };   //! Pr3d     - Pr2d
template<> struct intTableLA3d< BL_RT0LT_3d,  BL_Pr_2d    >  { static const intClass intFlag = OLA3d_STD; };   //! RT0      - Pr2d
template<> struct intTableLA3d< BL_Qr_3d,     BL_Qr_2d    >  { static const intClass intFlag = OLA3d_STD; };   //! Qr3d     - Qr2d
template<> struct intTableLA3d< BL_RT0LOC_3d, BL_Pr_2d    >  { static const intClass intFlag = OLA3d_STD; };   //! RT0_LOC  - Pr2d
template<> struct intTableLA3d< BL_RT0LOC_3d, BL_RT0LT_2d >  { static const intClass intFlag = OLA3d_STD; };   //! RT0_LOC  - RT0LT_2d
template<> struct intTableLA3d< BL_RT0LT_3d,  BL_RT0LT_2d >  { static const intClass intFlag = OLA3d_STD; };   //! RT0LT_3d - RT0LT_2d
template<> struct intTableLA3d< BL_OSwPr_3d,  BL_OSwPr_2d >  { static const intClass intFlag = OLA3d_OSC; };   //! OSC_Pr3d - OSC_Pr2d
template<> struct intTableLA3d< BL_OSCRr_3d,  BL_OSwPr_2d >  { static const intClass intFlag = OLA3d_OSC; };   //! OSC_Cr3d - OSC_Pr2d

#endif
