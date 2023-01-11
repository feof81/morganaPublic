/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef INTTABLEFEL3D_HPP
#define INTTABLEFEL3D_HPP

#include "fePr3d.hpp"
#include "feQr3d.hpp"
#include "feRt0Loc3d.hpp"
#include "feSpectralLH3d.hpp"
#include "morganaFields.hpp"
#include "morganaIntegrator.hpp"


template<FEBaseLabel> struct intTableFEL3d;  //! Generic class

template<> struct intTableFEL3d< BL_Pr_3d >     { static const intClass intFlag = FEL3d_STD; };   //! Pr
template<> struct intTableFEL3d< BL_Qr_3d >     { static const intClass intFlag = FEL3d_STD; };   //! Qr
template<> struct intTableFEL3d< BL_RT0LT_3d >  { static const intClass intFlag = FEL3d_STD; };   //! RT0
template<> struct intTableFEL3d< BL_RT0LOC_3d > { static const intClass intFlag = FEL3d_STD; };   //! RT0-LOC
template<> struct intTableFEL3d< BL_OSwPr_3d >  { static const intClass intFlag = FEL3d_OSC; };   //! Pr-OSC
template<> struct intTableFEL3d< BL_OSCRr_3d >  { static const intClass intFlag = FEL3d_OSC; };   //! Cr-OSC

#endif
