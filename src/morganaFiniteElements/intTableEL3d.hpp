/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef INTTABLEEL3D_HPP
#define INTTABLEEL3D_HPP

#include "fePr3d.hpp"
#include "feQr3d.hpp"
#include "feOscPr3d.hpp"
#include "feSpectralLH3d.hpp"
#include "morganaFields.hpp"
#include "morganaIntegrator.hpp"


template<FEBaseLabel, FEBaseLabel > struct intTableEL3d;  //! Generic class

template<> struct intTableEL3d< BL_Pr_3d,      BL_Pr_3d     >  { static const intClass intFlag = OEL3d_STD; };   //! Pr - Pr
template<> struct intTableEL3d< BL_Qr_3d,      BL_Qr_3d     >  { static const intClass intFlag = OEL3d_STD; };   //! Qr - Qr
template<> struct intTableEL3d< BL_SK_3d,      BL_SK_3d     >  { static const intClass intFlag = OEL3d_SKT; };   //! Spectral - Spectral
template<> struct intTableEL3d< BL_RT0LT_3d,   BL_Pr_3d     >  { static const intClass intFlag = OEL3d_STD; };   //! RT0 - Pr
template<> struct intTableEL3d< BL_RT0LT_3d,   BL_RT0LT_3d  >  { static const intClass intFlag = OEL3d_STD; };   //! RT0 - RT0
template<> struct intTableEL3d< BL_RT0LOC_3d,  BL_Pr_3d     >  { static const intClass intFlag = OEL3d_STD; };   //! RT0LOC - Pr
template<> struct intTableEL3d< BL_RT0LOC_3d,  BL_RT0LOC_3d >  { static const intClass intFlag = OEL3d_STD; };   //! RT0LOC - RT0
template<> struct intTableEL3d< BL_OSwPr_3d,   BL_OSwPr_3d  >  { static const intClass intFlag = OEL3d_OSC; };   //! OscPr - OscPr
template<> struct intTableEL3d< BL_OSCRr_3d,   BL_OSCRr_3d  >  { static const intClass intFlag = OEL3d_OSC; };   //! OscCr - OscCr

#endif
