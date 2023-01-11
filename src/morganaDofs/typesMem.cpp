/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "typesMem.h"


size_t dynamicSizeOf(const bool & A)
{ return(sizeof(A)); }

size_t dynamicSizeOf(const UInt & A)
{ return(sizeof(A)); }

size_t dynamicSizeOf(const Real & A)
{ return(sizeof(A)); }

size_t dynamicSizeOf(const point2d & A)
{ return(A.memSize()); }

size_t dynamicSizeOf(const point3d & A)
{ return(A.memSize()); }

size_t dynamicSizeOf(const tensor2d & A)
{ return(A.memSize()); }

size_t dynamicSizeOf(const tensor3d & A)
{ return(A.memSize()); }

size_t dynamicSizeOf(const stateVector & A)
{ return(A.memSize()); }

size_t dynamicSizeOf(const stateMatrix & A)
{ return(A.memSize()); }

template<size_t N>
size_t dynamicSizeOf(const staticVector<N> & A)
{ return(A.memSize()); }
