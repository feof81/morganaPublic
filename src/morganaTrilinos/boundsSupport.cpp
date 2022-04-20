/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "boundsSupport.h"
#include "epetraBlockMap_to_epetraMap.hpp"

boundsSupport::
boundsSupport()
{
  pid = 0;
  commDevLoaded = false;
}

boundsSupport::
boundsSupport(const Epetra_MpiComm & CommDev)
{
  pid = 0;
  commDevLoaded = true;
  
  commDev = Teuchos::rcpFromRef(CommDev);
}

boundsSupport::
boundsSupport(const RCP_MPICOMM & CommDev)
{
  pid = 0;
  commDevLoaded = true;
  
  commDev = CommDev;
}

void
boundsSupport::
setCommDev(const Epetra_MpiComm & CommDev)
{
  pid = 0;
  commDevLoaded = true;
  
  commDev = Teuchos::rcpFromRef(CommDev);
}

void
boundsSupport::
setCommDev(const RCP_MPICOMM & CommDev)
{
  pid = 0;
  commDevLoaded = true;
  
  commDev = CommDev;
}

void
boundsSupport::
increment()
{
  assert(commDevLoaded);
  
  pid++;
  
  if(pid >= commDev->NumProc())
  { pid = 0; }
  
  commDev->Barrier();
}

int
boundsSupport::
getPid()
{
  return(pid);
}

Epetra_Map
boundsSupport::
getMapRef()
{
  assert(commDevLoaded);
  
  int numRow;
  int elements[1];
  
  if(pid == commDev->MyPID()) { numRow = 1; }
  else                        { numRow = 0; }
  
  elements[0] = 0;
  
  Epetra_Map map(-1,numRow,elements,0,*commDev);
  return(map);
}

boundsSupport::RCP_BLOCKMAP
boundsSupport::
getMapRcp()
{
  assert(commDevLoaded);
  
  int numRow;
  int elements[1];
  
  if(pid == commDev->MyPID()) { numRow = 1; }
  else                        { numRow = 0; }
  
  elements[0] = 0;
  
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(-1,numRow,elements,0,*commDev));
  return(map);
}

Epetra_MultiVector
boundsSupport::
getSourceRef(const Real & value)
{
  assert(commDevLoaded);
  
  int numRow;
  int elements[1];
  
  if(pid == commDev->MyPID()) { numRow = 1; }
  else                        { numRow = 0; }
    
  elements[0] = 0;
  
  Epetra_Map map(-1,numRow,elements,0,*commDev);
  Epetra_MultiVector vector(map,1);
  
  int newRow = 0;
  vector.SumIntoGlobalValue(newRow,0,double(value));
  return(vector);
}

boundsSupport::RCP_VECTOR
boundsSupport::
getSourceRcp(const Real & value)
{
  assert(commDevLoaded);
  
  int numRow;
  int elements[1];
  
  if(pid == commDev->MyPID()) { numRow = 1; }
  else                        { numRow = 0; }
  
  elements[0] = 0;  
  
  Epetra_Map map(-1,numRow,elements,0,*commDev);
  Teuchos::RCP<Epetra_MultiVector> vector = Teuchos::rcp(new Epetra_MultiVector(map,1));
  
  int newRow = 0;
  vector->SumIntoGlobalValue(newRow,0,double(value));
  return(vector);
}

Epetra_FECrsMatrix 
boundsSupport::
vector_to_matrix(const Epetra_MultiVector & vector)
{
  assert(commDevLoaded);
  
  //New maps
  Epetra_Map rowMap = epetraBlockMap_to_epetraMap(vector.Map());
  Epetra_Map colMap = getMapRef();
  
  //Extract data
  UInt rowSource = vector.MyLength();
  double values[rowSource];
  vector.ExtractCopy(values,0);
  int * indices = rowMap.MyGlobalElements();
  
  //New matrix
  Epetra_FECrsMatrix matrix(Copy,rowMap,1);
  
  int GlobalRow, NumEntries = 1;
  int I = 0;
  double V;
  
  for(UInt i=0; i < rowSource; ++i)
  {
    GlobalRow = indices[i];
    V         = values[i];
    
    if(V != 0.0)
    { matrix.InsertGlobalValues(GlobalRow,NumEntries,&V,&I); }
  }
  
  //Fill complete
  matrix.FillComplete(colMap,rowMap);
  return(matrix);
}

boundsSupport::RCP_MATRIX
boundsSupport::
vector_to_matrix(const RCP_VECTOR & vector)
{
  assert(commDevLoaded);
  
  //New maps
  Epetra_Map rowMap = epetraBlockMap_to_epetraMap(vector->Map());
  Epetra_Map colMap = getMapRef();
  
  //Extract data
  UInt rowSource = vector->MyLength();
  double values[rowSource];
  vector->ExtractCopy(values,0);
  int * indices = rowMap.MyGlobalElements();
  
  //New matrix
  Teuchos::RCP<Epetra_FECrsMatrix> matrix = Teuchos::rcp(new Epetra_FECrsMatrix(Copy,rowMap,1));
  
  int GlobalRow, NumEntries = 1;
  int I = 0;
  double V;
  
  for(UInt i=0; i < rowSource; ++i)
  {
    GlobalRow = indices[i];
    V         = values[i];
    
    if(V != 0.0)
    { matrix->InsertGlobalValues(GlobalRow,NumEntries,&V,&I); }
  }
  
  //Fill complete
  matrix->FillComplete(colMap,rowMap);
  return(matrix);
}
