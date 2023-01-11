/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef TRAITSINTERPOLATOR_HPP
#define TRAITSINTERPOLATOR_HPP

#include "traitsBasic.h"

#include "mesh1d.hpp"
#include "mesh2d.hpp"
#include "mesh3d.hpp"

#include "feStaticField1d.hpp"
#include "feStaticField2d.hpp"
#include "feStaticField3d.hpp"

#include "feDynamicField1d.hpp"
#include "feDynamicField2d.hpp"
#include "feDynamicField3d.hpp"

#include "search1dA.hpp"
#include "search2dA.hpp"
#include "search3dA.hpp"

#include "geoMapSupport1d.hpp"
#include "geoMapSupport2d.hpp"
#include "geoMapSupport3d.hpp"

#include "feStaticField1dGlobalManip.hpp"
#include "feStaticField2dGlobalManip.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "feDynamicField1dGlobalManip.hpp"
#include "feDynamicField2dGlobalManip.hpp"
#include "feDynamicField3dGlobalManip.hpp"


/*! Traits for the simple interpolation */
template<typename SOURCEFIELD> class interpolatorTrait;



//One dymension------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
class interpolatorTrait<feStaticField1d<FETYPE,DOFTYPE,ORDER,MODE> >
{
  public:
    typedef typename feStaticField1d<FETYPE,DOFTYPE,ORDER,MODE>::PMAPTYPE  PMAPTYPE;
    typedef typename feStaticField1d<FETYPE,DOFTYPE,ORDER,MODE>::GEOSHAPE  GEOSHAPE;
    typedef mesh1d<GEOSHAPE,PMAPTYPE,PMAPTYPE>                             MESH;
    typedef connect1d<GEOSHAPE,PMAPTYPE,PMAPTYPE>                          CONNECT;
    typedef search1dA<GEOSHAPE,PMAPTYPE,PMAPTYPE>                          SEARCH;
    typedef geoMapSupport1d<GEOSHAPE>                                      SUPPORT;
    typedef feStaticField1dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>          MANIP;
};

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER, dmd1d_mode MODE>
class interpolatorTrait<feDynamicField1d<FETYPE,DOFTYPE,ORDER,MODE> >
{
  public:
    typedef typename feDynamicField1d<FETYPE,DOFTYPE,ORDER,MODE>::PMAPTYPE  PMAPTYPE;
    typedef typename feDynamicField1d<FETYPE,DOFTYPE,ORDER,MODE>::GEOSHAPE  GEOSHAPE;
    typedef mesh1d<GEOSHAPE,PMAPTYPE,PMAPTYPE>                              MESH;
    typedef connect1d<GEOSHAPE,PMAPTYPE,PMAPTYPE>                           CONNECT;
    typedef search1dA<GEOSHAPE,PMAPTYPE,PMAPTYPE>                           SEARCH;
    typedef geoMapSupport1d<GEOSHAPE>                                       SUPPORT;
    typedef feDynamicField1dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>          MANIP;
};


//Two dymensions-----------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
class interpolatorTrait<feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE> >
{
  public:
    typedef typename feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::PMAPTYPE  PMAPTYPE;
    typedef typename feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::GEOSHAPE  GEOSHAPE;
    typedef mesh2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>                             MESH;
    typedef connect2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>                          CONNECT;
    typedef search2dA<GEOSHAPE,PMAPTYPE,PMAPTYPE>                          SEARCH;
    typedef geoMapSupport2d<GEOSHAPE>                                      SUPPORT;
    typedef feStaticField2dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>          MANIP;
};

template<typename FETYPE, typename DOFTYPE, dmd2d_order ORDER, dmd2d_mode MODE>
class interpolatorTrait<feDynamicField2d<FETYPE,DOFTYPE,ORDER,MODE> >
{
  public:
    typedef typename feDynamicField2d<FETYPE,DOFTYPE,ORDER,MODE>::PMAPTYPE  PMAPTYPE;
    typedef typename feDynamicField2d<FETYPE,DOFTYPE,ORDER,MODE>::GEOSHAPE  GEOSHAPE;
    typedef mesh2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>                              MESH;
    typedef connect2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>                           CONNECT;
    typedef search2dA<GEOSHAPE,PMAPTYPE,PMAPTYPE>                           SEARCH;
    typedef geoMapSupport2d<GEOSHAPE>                                       SUPPORT;
    typedef feDynamicField2dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>          MANIP;
};


//Three dymensions---------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
class interpolatorTrait<feStaticField3d<FETYPE,DOFTYPE,ORDER,MODE> >
{
  public:
    typedef typename feStaticField3d<FETYPE,DOFTYPE,ORDER,MODE>::PMAPTYPE  PMAPTYPE;
    typedef typename feStaticField3d<FETYPE,DOFTYPE,ORDER,MODE>::GEOSHAPE  GEOSHAPE;
    typedef mesh3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>                             MESH;
    typedef connect3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>                          CONNECT;
    typedef search3dA<GEOSHAPE,PMAPTYPE,PMAPTYPE>                          SEARCH;
    typedef geoMapSupport3d<GEOSHAPE>                                      SUPPORT;
    typedef feStaticField3dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>          MANIP;
};

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER, dmd3d_mode MODE>
class interpolatorTrait<feDynamicField3d<FETYPE,DOFTYPE,ORDER,MODE> >
{
  public:
    typedef typename feDynamicField3d<FETYPE,DOFTYPE,ORDER,MODE>::PMAPTYPE  PMAPTYPE;
    typedef typename feDynamicField3d<FETYPE,DOFTYPE,ORDER,MODE>::GEOSHAPE  GEOSHAPE;
    typedef mesh3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>                              MESH;
    typedef connect3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>                           CONNECT;
    typedef search3dA<GEOSHAPE,PMAPTYPE,PMAPTYPE>                           SEARCH;
    typedef geoMapSupport3d<GEOSHAPE>                                       SUPPORT;
    typedef feDynamicField3dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>          MANIP;
};

#endif
