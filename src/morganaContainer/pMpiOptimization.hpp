/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PMPIOPTIMIZATION_HPP
#define PMPIOPTIMIZATION_HPP

#include "pMapItem.h"
#include "pMapItemShare.h"
#include "pMapItemSendRecv.h"

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

/*! Maximum length in bytes of comm-data */
static const UInt mpiMaxSize = 1e7;

/*! Maximum number of channels */
static const UInt mpiLayer = 100;

/*! The optimizations of the pMapItems */
namespace boost
{
  namespace mpi
  {
    template<>
    struct is_mpi_datatype<pMapItem> : mpl::true_{ };
    
    template<>
    struct is_mpi_datatype<pMapItemShare> : mpl::true_{ };
    
    template<>
    struct is_mpi_datatype<pMapItemSendRecv> : mpl::true_{ };
  }
}

#endif
