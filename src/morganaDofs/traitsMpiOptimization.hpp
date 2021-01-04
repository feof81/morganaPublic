/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef TRAITSMPIOPTIMIZATION_HPP
#define TRAITSMPIOPTIMIZATION_HPP

#include "point2d.h"
#include "point3d.h"

#include "tensor2d.h"
#include "tensor3d.h"

#include "simpleFormats.hpp"
#include "staticVector.hpp"
#include "komplex.h"

#include <boost/serialization/access.hpp>
#include <boost/mpi/collectives/all_gather.hpp>
#include <boost/mpi.hpp>

using namespace std;
using namespace boost::mpi;


/*! Optimization of the types for parallel transmission */
namespace boost
{
  namespace mpi
  {
    template<>
    struct is_mpi_datatype<point3d> : mpl::true_{ };
    
    template<>
    struct is_mpi_datatype<point2d> : mpl::true_{ };
    
    template<>
    struct is_mpi_datatype<tensor3d> : mpl::true_{ };
    
    template<>
    struct is_mpi_datatype<tensor2d> : mpl::true_{ };
    
    template<>
    struct is_mpi_datatype<komplex> : mpl::true_{ };
    
    template<size_t N>
    struct is_mpi_datatype<staticVector<N> > : mpl::true_{ };
  }
}


/*! Boost Mpi collectives fix */
template<typename T>
void
morgana_all_gather(const communicator& comm, const T& in_value, std::vector<T>& out_values)
{
  //Allocations
  sVect<request> reqs(2 * comm.size() - 2);
  out_values.resize(comm.size());
  
  //Sending
  UInt j = 1;
  
  for(int k=0; k < comm.size(); ++k)
  {
    if(k != comm.rank())
    {
      reqs(j) = comm.isend(k, comm.rank(), in_value);
      ++j;
    }
  }
  
  //Reciving 
  for(int k=0; k < comm.size(); ++k)
  {
    if(k != comm.rank())
    {
      reqs(j) = comm.irecv(k, k, out_values[k]);
      ++j;
    }
    else
    {
      out_values[k] = in_value;
    }
  }
  
  //Waiting
  wait_all(reqs.begin(),reqs.end());
}

#endif

