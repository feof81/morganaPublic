/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef LINEARMAP_DISTARRAY_HPP
#define LINEARMAP_DISTARRAY_HPP

#include "Epetra_Map.h"
#include "EpetraExt_DistArray.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Import.h"



/*! Redistribute the elements building a linear distribution map */
template<typename TYPE>
EpetraExt::DistArray<TYPE> linearMapDistArray(const EpetraExt::DistArray<TYPE> & distArray, const Epetra_MpiComm & epetraComm)
{
  const Epetra_BlockMap & OriginalMap = distArray.Map();
  assert(OriginalMap.UniqueGIDs());
  
  Epetra_Map    LinearMap(OriginalMap.NumGlobalElements(), OriginalMap.IndexBase(), epetraComm);
  Epetra_Import Importer(LinearMap, OriginalMap);
  
  EpetraExt::DistArray<TYPE> LinearX(LinearMap, distArray.RowSize());
  LinearX.Import(distArray, Importer, Add);
  
  return(LinearX);
}

#endif
