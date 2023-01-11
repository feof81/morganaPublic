/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef TYPESMEM_H
#define TYPESMEM_H

#include "point2d.h"
#include "point3d.h"
#include "tensor2d.h"
#include "tensor3d.h"
#include "stateVector.h"
#include "stateMatrix.h"
#include "staticVector.hpp"


size_t dynamicSizeOf(const bool & A);
size_t dynamicSizeOf(const UInt & A);
size_t dynamicSizeOf(const Real & A);
size_t dynamicSizeOf(const point2d & A);
size_t dynamicSizeOf(const point3d & A);
size_t dynamicSizeOf(const tensor2d & A);
size_t dynamicSizeOf(const tensor3d & A);
size_t dynamicSizeOf(const stateVector & A);
size_t dynamicSizeOf(const stateMatrix & A);

template<size_t N>
size_t dynamicSizeOf(const staticVector<N> & A);

#endif
