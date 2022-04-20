/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef EPETRABLOCKMAP_TO_EPETRAMAP_HPP
#define EPETRABLOCKMAP_TO_EPETRAMAP_HPP

#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_RCP.hpp"

#include "Epetra_ConfigDefs.h"
#include "Epetra_MpiComm.h"

#include "Epetra_DataAccess.h"
#include "Epetra_Map.h"

#include "typesInterface.hpp"


inline Epetra_Map epetraBlockMap_to_epetraMap(const Epetra_BlockMap & sourceMap)
{
  int       size = sourceMap.NumMyElements();
  int * elements = sourceMap.MyGlobalElements();
  
  Epetra_Map newMap(-1, size, elements, 0, sourceMap.Comm());
  
  return(newMap);
}


#endif
