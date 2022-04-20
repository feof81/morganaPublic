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
#include "epetraOperatorSum.h"


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
  
  if(world.rank() == 0) { numRowA = 2; }
  else                  { numRowA = 2; }
  
  if(world.rank() == 0) { numColA = 3; }
  else                  { numColA = 3; }
  
  int elementsRowA[numRowA];
  int elementsColA[numColA];
  
  if(world.rank() == 0)
  {
    elementsRowA[0] = 0;
    elementsRowA[1] = 1;
    
    elementsColA[0] = 0;
    elementsColA[1] = 1;
    elementsColA[2] = 2;
  }
  else
  {
    elementsRowA[0] = 2;
    elementsRowA[1] = 3;
    
    elementsColA[0] = 1;
    elementsColA[1] = 2;
    elementsColA[2] = 3;
  }  
  
  Epetra_Map rowMapA(-1,numRowA,elementsRowA,0,epetraComm);
  Epetra_Map colMapA(-1,numColA,elementsColA,0,epetraComm);
  
  
  //Matrix matrixA-------------------------------------------------------------
  Teuchos::RCP<Epetra_CrsMatrix> matrixA(new Epetra_CrsMatrix(Copy, rowMapA, colMapA, 1));
  Teuchos::RCP<Epetra_CrsMatrix> matrixB(new Epetra_CrsMatrix(Copy, rowMapA, colMapA, 1));
  
  int length = 4;
  double valuesA[length];
  int numEntriesA, rowA, indicesA[length];
  
  
  if(world.rank() == 0)
  {
    //MatrixA
    rowA = 0;
    numEntriesA = 3;
    indicesA[0] = 0; indicesA[1] = 1; indicesA[2] = 2;
    valuesA[0]  = 1; valuesA[1]  = 1; valuesA[2]  = 1;
    matrixA->InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
    
    rowA = 1;
    numEntriesA = 3;
    indicesA[0] = 0; indicesA[1] = 1; indicesA[2] = 2;
    valuesA[0]  = 0; valuesA[1]  = 1; valuesA[2]  = 1;
    matrixA->InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
    
    //MatrixB
    rowA = 0;
    indicesA[0] = 0; indicesA[1] = 1; indicesA[2] = 2;
    valuesA[0]  = 1; valuesA[1]  = 0; valuesA[2]  = 0;
    matrixB->InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
    
    rowA = 1;
    indicesA[0] = 0; indicesA[1] = 1; indicesA[2] = 2;
    valuesA[0]  = 0; valuesA[1]  = 1; valuesA[2]  = 0;
    matrixB->InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
  }
  else
  {
    //MatrixA
    rowA = 2;
    numEntriesA = 3;
    indicesA[0] =  1; indicesA[1] = 2; indicesA[2] = 3;
    valuesA[0]  = -1; valuesA[1]  = 1; valuesA[2]  = 1;
    matrixA->InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
    
    rowA = 3;
    numEntriesA = 3;
    indicesA[0] = 1; indicesA[1] = 2; indicesA[2] = 3;
    valuesA[0]  = 1; valuesA[1]  = 1; valuesA[2]  = 1;
    matrixA->InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
    
    //MatrixB
    rowA = 2;
    indicesA[0] = 1; indicesA[1] = 2; indicesA[2] = 3;
    valuesA[0]  = 0; valuesA[1]  = 1; valuesA[2]  = 0;
    matrixB->InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
    
    rowA = 3;
    indicesA[0] = 1; indicesA[1] = 2; indicesA[2] = 3;
    valuesA[0]  = 0; valuesA[1]  = 0; valuesA[2]  = 1;
    matrixB->InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
  }
  
  matrixA->FillComplete(rowMapA, rowMapA);
  matrixB->FillComplete(rowMapA, rowMapA);
  
  //Build the vectors----------------------------------------------------------
  Teuchos::RCP<Epetra_MultiVector> X(new Epetra_MultiVector(rowMapA,1));
  Teuchos::RCP<Epetra_MultiVector> Y(new Epetra_MultiVector(rowMapA,1));
  
  if(world.rank() == 0)
  {
    X->SumIntoMyValue(0,0,1.0);
    X->SumIntoMyValue(1,0,2.0);
  }
  else
  {
    X->SumIntoMyValue(0,0,3.0);
    X->SumIntoMyValue(1,0,4.0);
  }
  
  //Operator sum---------------------------------------------------------------
  epetraOperatorSum S(matrixA,epetraOperatorSum::op(matrixB,matrixB));
  S.Apply(*X,*Y);
  
  cout << *matrixB << endl;
  cout << "Y : " << *Y << endl;
}

