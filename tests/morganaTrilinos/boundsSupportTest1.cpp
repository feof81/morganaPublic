/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <cmath>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_FEVector.h"

#include "typesInterface.hpp"
#include "pMapItem.h"
#include "pMapItemShare.h"
#include "pVect.hpp"
#include "pVectManip.hpp"
#include "pVectComm.hpp"
#include "pVectGlobalManip.hpp"

#include "epetraVector_to_pVect.h"
#include "epetraBlockMatrixAssemble.h"
#include "boundsSupport.h"



//! Run with two processors
int main(int argc, char *argv[])
{
  //Startup--------------------------------------------------------------------
  environment env(argc,argv);
  
  communicator world;
  Epetra_MpiComm epetraComm(MPI_COMM_WORLD);
  
  assert(world.size() == 2);
  
  
  //Maps-----------------------------------------------------------------------
  int numRowA, numColA;
  
  if(world.rank() == 0) { numRowA = 3; }
  else                  { numRowA = 2; }
  
  if(world.rank() == 0) { numColA = 3; }
  else                  { numColA = 3; }
  
  int elementsRowA[numRowA];
  int elementsColA[numColA];
  
  if(world.rank() == 0)
  {
    elementsRowA[0] = 2;
    elementsRowA[1] = 3;
    elementsRowA[2] = 4;
     
    elementsColA[0] = 0;
    elementsColA[1] = 1;
    elementsColA[2] = 2;
  }
  else
  {
    elementsRowA[0] = 0;
    elementsRowA[1] = 1;
    
    elementsColA[0] = 0;
    elementsColA[1] = 1;
    elementsColA[2] = 2;
  }

  Epetra_Map rowMapA(-1,numRowA,elementsRowA,0,epetraComm);
  Epetra_Map colMapA(-1,numColA,elementsColA,0,epetraComm);
  
  
  //Matrix matrixAA------------------------------------------------------------
  Epetra_CrsMatrix matrixAA(Copy, rowMapA, colMapA, 1);
  
  int length = 3;
  double valuesAA[length];
  int numEntriesAA, rowAA, indicesAA[length];
  
  numEntriesAA = 3;
  indicesAA[0] = 0; indicesAA[1] = 1; indicesAA[2] = 2;
  
  if(world.rank() == 0)
  {
    rowAA = 2;
    valuesAA[0] = 0; valuesAA[1] = 1; valuesAA[2] = 2;
    matrixAA.InsertGlobalValues(rowAA, numEntriesAA, valuesAA, indicesAA);
    
    rowAA = 3;
    valuesAA[0] = 3; valuesAA[1] = 2; valuesAA[2] = 1;
    matrixAA.InsertGlobalValues(rowAA, numEntriesAA, valuesAA, indicesAA);
    
    rowAA = 4;
    valuesAA[0] = 2; valuesAA[1] = 3; valuesAA[2] = 4;
    matrixAA.InsertGlobalValues(rowAA, numEntriesAA, valuesAA, indicesAA);
  }
  else
  {
    rowAA = 0;
    valuesAA[0] = 1; valuesAA[1] = 2; valuesAA[2] = 1;
    matrixAA.InsertGlobalValues(rowAA, numEntriesAA, valuesAA, indicesAA);
    
    rowAA = 1;
    valuesAA[0] = 2; valuesAA[1] = 1; valuesAA[2] = 2;
    matrixAA.InsertGlobalValues(rowAA, numEntriesAA, valuesAA, indicesAA);
  }
  
  matrixAA.FillComplete (colMapA, rowMapA);
  
 
  //Vector------------------------------------------------------------
  Epetra_MultiVector vectorA(rowMapA,1);
  
  int I;
  double tot;
  
  if(world.rank() == 0)
  {
    I = 2; tot = 3;
    vectorA.SumIntoGlobalValue(I,0,tot);
    
    I = 3; tot = 4;
    vectorA.SumIntoGlobalValue(I,0,tot);
    
    I = 4; tot = 5;
    vectorA.SumIntoGlobalValue(I,0,tot);
  }
  else
  {
    I = 0; tot = 1;
    vectorA.SumIntoGlobalValue(I,0,tot);
  
    I = 1; tot = 2;
    vectorA.SumIntoGlobalValue(I,0,tot);
  }
  
  
  boundsSupport bounds(epetraComm);
  cout << "Matrix" << endl;
  cout << bounds.vector_to_matrix(vectorA) << endl;
  
  cout << "Vector" << endl;
  cout << bounds.getSourceRef(4.0) << endl;
}
