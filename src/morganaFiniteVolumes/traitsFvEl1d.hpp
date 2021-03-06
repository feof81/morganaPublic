/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef TRAITSFVEL1D_HPP
#define TRAITSFVEL1D_HPP

#include "pMapItemShare.h"
#include "fePr1d.hpp"
#include "feStaticField1d.hpp"


/*! Empty forward declaration */
template<typename GEOSHAPE, typename DOFTYPE> class traitsFvEl1d;


/*! Linear Tetra specialization */
template<typename DOFTYPE> 
class traitsFvEl1d<linearLine,DOFTYPE>
{
  public:
    typedef pMapItemShare                                                        PMAPTYPE;
    typedef fePr1d<1,PMAPTYPE>                                                   FLUX_FETYPE;
    typedef feStaticField1d<FLUX_FETYPE,DOFTYPE,dms1d_vectMajor,dms1d_geoIdMode> FLUX_FIELD;
    
  public:
    static const UInt intPrec = 1;
};

#endif
