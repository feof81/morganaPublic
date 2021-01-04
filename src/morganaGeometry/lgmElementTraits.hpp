/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef LGMELEMENTTRAITS_HPP
#define LGMELEMENTTRAITS_HPP

#include "geoMapSupport1d.hpp"
#include "geoMapSupport2d.hpp"
#include "geoMapSupport3d.hpp"

#include "mesh3d.hpp"
#include "mesh2d.hpp"
#include "mesh1d.hpp"

#include "connect3d.hpp"
#include "connect2d.hpp"
#include "connect1d.hpp"


template<typename GEOSHAPE, typename ELMAP, typename NODEMAP, Int DIM>
class lgmElementTraits
{
};

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class lgmElementTraits<GEOSHAPE,ELMAP,NODEMAP,3>
{
  public:
    typedef geoMapSupport3d<GEOSHAPE>         GEOSUPPORT;
    typedef mesh3d<GEOSHAPE,ELMAP,NODEMAP>    MESH;
    typedef connect3d<GEOSHAPE,ELMAP,NODEMAP> CONNECT;
};

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class lgmElementTraits<GEOSHAPE,ELMAP,NODEMAP,2>
{
  public:
    typedef geoMapSupport2d<GEOSHAPE>         GEOSUPPORT;
    typedef mesh2d<GEOSHAPE,ELMAP,NODEMAP>    MESH;
    typedef connect2d<GEOSHAPE,ELMAP,NODEMAP> CONNECT;
};

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class lgmElementTraits<GEOSHAPE,ELMAP,NODEMAP,1>
{
  public:
    typedef geoMapSupport1d<GEOSHAPE>         GEOSUPPORT;
    typedef mesh1d<GEOSHAPE,ELMAP,NODEMAP>    MESH;
    typedef connect1d<GEOSHAPE,ELMAP,NODEMAP> CONNECT;
};

#endif
