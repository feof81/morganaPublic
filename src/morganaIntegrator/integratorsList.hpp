/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTEGRATORSLIST_HPP
#define INTEGRATORSLIST_HPP

#include "morganaIntegrator.hpp"

#include "integratorOLA3d_STD.hpp"
#include "integratorOLA3d_STC.hpp"
#include "integratorOLA3d_OSC.hpp"
#include "integratorOHY3d_STD.hpp"
#include "integratorOHY3d_STC.hpp"
#include "integratorOHY3d_OSC.hpp"
#include "integratorOEL3d_STD.hpp"
#include "integratorOEL3d_STC.hpp"
#include "integratorOEL3d_SKT.hpp"
#include "integratorOEL3d_OSC.hpp"
#include "integratorOBB3d_STC.hpp"
#include "integratorOBB3d_STD.hpp"
#include "integratorOBB3d_OSC.hpp"

#include "integratorFLA3d_STD.hpp"
#include "integratorFLA3d_STC.hpp"
#include "integratorFLA3d_OSC.hpp"
#include "integratorFEL3d_STD.hpp"
#include "integratorFEL3d_STC.hpp"
#include "integratorFEL3d_OSC.hpp"
#include "integratorFBF3d_STD.hpp"
#include "integratorFBF3d_STC.hpp"
#include "integratorFBF3d_OSC.hpp"
#include "integratorFHY3d_STD.hpp"
#include "integratorFHY3d_STC.hpp"
#include "integratorFHY3d_OSC.hpp"

#include "integratorOLA2d_STD.hpp"
#include "integratorOHY2d_STD.hpp"
#include "integratorOEL2d_STD.hpp"
#include "integratorOEL2d_STC.hpp"
#include "integratorOEL2d_SKT.hpp"
#include "integratorOBB2d_STC.hpp"
#include "integratorOBB2d_STD.hpp"

#include "integratorFLA2d_STD.hpp"
#include "integratorFEL2d_STD.hpp"
#include "integratorFBF2d_STD.hpp"

#include "integratorOEL1d_STD.hpp"



template<intClass CLASS, typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList;



//_________________________________________________________________________________________________
// ELEMENTS 3D - OPERATOR
//-------------------------------------------------------------------------------------------------
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OLA3d_STD, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOLA3d_STD<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OLA3d_STC, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOLA3d_STC<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OLA3d_OSC, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOLA3d_OSC<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OHY3d_STD, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOHY3d_STD<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OHY3d_STC, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOHY3d_STC<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OHY3d_OSC, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOHY3d_OSC<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OEL3d_STD, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOEL3d_STD<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OEL3d_STC, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOEL3d_STC<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OEL3d_SKT, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOEL3d_SKT<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OEL3d_OSC, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOEL3d_OSC<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OBB3d_STD, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOBB3d_STD<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OBB3d_STC, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOBB3d_STC<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OBB3d_OSC, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOBB3d_OSC<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};



//_________________________________________________________________________________________________
// ELEMENTS 3D - FUNCTIONALS
//-------------------------------------------------------------------------------------------------
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FLA3d_STD, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFLA3d_STD<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FLA3d_STC, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFLA3d_STC<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FLA3d_OSC, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFLA3d_OSC<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FEL3d_STD, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFEL3d_STD<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FEL3d_STC, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFEL3d_STC<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FEL3d_OSC, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFEL3d_OSC<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FBF3d_STD, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFBF3d_STD<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FBF3d_STC, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFBF3d_STC<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FBF3d_OSC, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFBF3d_OSC<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FHY3d_STD, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFHY3d_STD<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FHY3d_STC, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFHY3d_STC<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FHY3d_OSC, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFHY3d_OSC<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};



//_________________________________________________________________________________________________
// ELEMENTS 2D - OPERATOR
//-------------------------------------------------------------------------------------------------
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OLA2d_STD, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOLA2d_STD<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OHY2d_STD, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOHY2d_STD<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OEL2d_STD, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOEL2d_STD<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OEL2d_STC, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOEL2d_STC<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OEL2d_SKT, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOEL2d_SKT<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OBB2d_STD, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOBB2d_STD<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};

template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OBB2d_STC, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOBB2d_STC<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};



//_________________________________________________________________________________________________
// ELEMENTS 2D - FUNCTIONALS
//-------------------------------------------------------------------------------------------------
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FLA2d_STD, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFLA2d_STD<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FEL2d_STD, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFEL2d_STD<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
struct integratorsList<FBF2d_STD, FUNCTIONAL, TYPE, PRECISION>
{
  typedef integratorFBF2d_STD<FUNCTIONAL,TYPE,PRECISION> INTEGRATOR;
};



//_________________________________________________________________________________________________
// ELEMENTS 1D - OPERATOR
//-------------------------------------------------------------------------------------------------
template<typename OPERATOR, intTypes TYPE, UInt PRECISION>
struct integratorsList<OEL1d_STD, OPERATOR, TYPE, PRECISION>
{
  typedef integratorOEL1d_STD<OPERATOR,TYPE,PRECISION> INTEGRATOR;
};



#endif
