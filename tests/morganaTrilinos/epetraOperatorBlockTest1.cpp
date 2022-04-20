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
#include "epetraBlockVectorAssemble.h"
#include "epetraBlockMatrixAssemble.h"
#include "epetraOperatorBlock.h"



//! Run with two processors
int main(int argc, char *argv[])
{
  //Startup----------------------------------------------------------------------------------------
  environment env(argc,argv);
  
  communicator world;
  Epetra_MpiComm epetraComm(MPI_COMM_WORLD);
  
  assert(world.size() == 2);
  
  
  //Maps-------------------------------------------------------------------------------------------
  int numRowA, numRowB;
  int numColA, numColB;
  
  if(world.rank() == 0) { numRowA = 3; numRowB = 2; }
  else                  { numRowA = 2; numRowB = 2; }
  
  if(world.rank() == 0) { numColA = 3; numColB = 2; }
  else                  { numColA = 3; numColB = 2; }
  
  int elementsRowA[numRowA];
  int elementsRowB[numRowB];
  int elementsColA[numColA];
  int elementsColB[numColB];
  
  if(world.rank() == 0)
  {
    elementsRowA[0] = 2;
    elementsRowA[1] = 3;
    elementsRowA[2] = 4;
    
    elementsRowB[0] = 2;
    elementsRowB[1] = 3;
    
    elementsColA[0] = 0;
    elementsColA[1] = 1;
    elementsColA[2] = 2;
    
    elementsColB[0] = 0;
    elementsColB[1] = 1;
  }
  else
  {
    elementsRowA[0] = 0;
    elementsRowA[1] = 1;
    
    elementsRowB[0] = 0;
    elementsRowB[1] = 1;
    
    elementsColA[0] = 0;
    elementsColA[1] = 1;
    elementsColA[2] = 2;
    
    elementsColB[0] = 0;
    elementsColB[1] = 1;
  }

  Epetra_Map rowMapA(-1,numRowA,elementsRowA,0,epetraComm);
  Epetra_Map rowMapB(-1,numRowB,elementsRowB,0,epetraComm);
  
  Epetra_Map colMapA(-1,numColA,elementsColA,0,epetraComm);
  Epetra_Map colMapB(-1,numColB,elementsColB,0,epetraComm);
  
  
  //Matrix matrixAA--------------------------------------------------------------------------------
  Teuchos::RCP<Epetra_CrsMatrix> matrixAA(new Epetra_CrsMatrix(Copy, rowMapA, colMapA, 1));
  
  int length = 3;
  double valuesAA[length];
  int numEntriesAA, rowAA, indicesAA[length];
  
  numEntriesAA = 3;
  indicesAA[0] = 0; indicesAA[1] = 1; indicesAA[2] = 2;
  
  if(world.rank() == 0)
  {
    rowAA = 2;
    valuesAA[0] = 0; valuesAA[1] = 1; valuesAA[2] = 2;
    matrixAA->InsertGlobalValues(rowAA, numEntriesAA, valuesAA, indicesAA);
    
    rowAA = 3;
    valuesAA[0] = 3; valuesAA[1] = 2; valuesAA[2] = 1;
    matrixAA->InsertGlobalValues(rowAA, numEntriesAA, valuesAA, indicesAA);
    
    rowAA = 4;
    valuesAA[0] = 2; valuesAA[1] = 3; valuesAA[2] = 4;
    matrixAA->InsertGlobalValues(rowAA, numEntriesAA, valuesAA, indicesAA);
  }
  else
  {
    rowAA = 0;
    valuesAA[0] = 1; valuesAA[1] = 2; valuesAA[2] = 1;
    matrixAA->InsertGlobalValues(rowAA, numEntriesAA, valuesAA, indicesAA);
    
    rowAA = 1;
    valuesAA[0] = 2; valuesAA[1] = 1; valuesAA[2] = 2;
    matrixAA->InsertGlobalValues(rowAA, numEntriesAA, valuesAA, indicesAA);
  }
  
  matrixAA->FillComplete(rowMapA, rowMapA);
  
 
  
  //Matrix matrixAB--------------------------------------------------------------------------------
  Teuchos::RCP<Epetra_CrsMatrix> matrixAB(new Epetra_CrsMatrix(Copy, rowMapA, colMapB, 1));
  
  double valuesAB[length];
  int numEntriesAB, rowAB, indicesAB[length];
  
  numEntriesAB = 2;
  indicesAB[0] = 0; indicesAB[1] = 1; indicesAB[2] = 2;
  
  if(world.rank() == 0)
  {
    rowAB = 2;
    valuesAB[0] = 3; valuesAB[1] = 2;
    matrixAB->InsertGlobalValues(rowAB, numEntriesAB, valuesAB, indicesAB);
    
    rowAB = 3;
    valuesAB[0] = 2; valuesAB[1] = 1;
    matrixAB->InsertGlobalValues(rowAB, numEntriesAB, valuesAB, indicesAB);
    
    rowAB = 4;
    valuesAB[0] = 1; valuesAB[1] = 0;
    matrixAB->InsertGlobalValues(rowAB, numEntriesAB, valuesAB, indicesAB);
  }
  else
  {
    rowAB = 0;
    valuesAB[0] = 5; valuesAB[1] = 4;
    matrixAB->InsertGlobalValues(rowAB, numEntriesAB, valuesAB, indicesAB);
    
    rowAB = 1;
    valuesAB[0] = 4; valuesAB[1] = 3;
    matrixAB->InsertGlobalValues(rowAB, numEntriesAB, valuesAB, indicesAB);
  }
  
  matrixAB->FillComplete(rowMapB, rowMapA);
  

  
  //Matrix matrixBA--------------------------------------------------------------------------------
  Teuchos::RCP<Epetra_CrsMatrix> matrixBA(new Epetra_CrsMatrix(Copy, rowMapB, colMapA, 1));
  
  double valuesBA[length];
  int numEntriesBA, rowBA, indicesBA[length];
  
  numEntriesBA = 3;
  indicesBA[0] = 0; indicesBA[1] = 1; indicesBA[2] = 2;
  
  
  if(world.rank() == 0)
  {
    rowBA = 2;
    valuesBA[0] = 1; valuesBA[1] = 2; valuesBA[2] = 3;
    matrixBA->InsertGlobalValues(rowBA, numEntriesBA, valuesBA, indicesBA);
    
    rowBA = 3;
    valuesBA[0] = 4; valuesBA[1] = 5; valuesBA[2] = 6;
    matrixBA->InsertGlobalValues(rowBA, numEntriesBA, valuesBA, indicesBA);
  }
  else
  {
    rowBA = 0;
    valuesBA[0] = 9; valuesBA[1] = 8; valuesBA[2] = 1;
    matrixBA->InsertGlobalValues(rowBA, numEntriesBA, valuesBA, indicesBA);
    
    rowBA = 1;
    valuesBA[0] = 2; valuesBA[1] = 4; valuesBA[2] = 6;
    matrixBA->InsertGlobalValues(rowBA, numEntriesBA, valuesBA, indicesBA);
  }
  
  matrixBA->FillComplete (rowMapA, rowMapB);
  
  
  
  //Matrix matrixBB--------------------------------------------------------------------------------
  Teuchos::RCP<Epetra_CrsMatrix> matrixBB(new Epetra_CrsMatrix(Copy, rowMapB, colMapB, 1));
  
  double valuesBB[length];
  int numEntriesBB, rowBB, indicesBB[length];
  
  numEntriesBB = 2;
  indicesBB[0] = 0; indicesBB[1] = 1;
  
  if(world.rank() == 0)
  {
    rowBB = 2;
    valuesBB[0] = 1; valuesBB[1] = 6;
    matrixBB->InsertGlobalValues(rowBB, numEntriesBB, valuesBB, indicesBB);
    
    rowBB = 3;
    valuesBB[0] = 5; valuesBB[1] = 4;
    matrixBB->InsertGlobalValues(rowBB, numEntriesBB, valuesBB, indicesBB);
  }
  else
  {
    rowBB = 0;
    valuesBB[0] = 3; valuesBB[1] = 4;
    matrixBB->InsertGlobalValues(rowBB, numEntriesBB, valuesBB, indicesBB);
    
    rowBB = 1;
    valuesBB[0] = 4; valuesBB[1] = 2;
    matrixBB->InsertGlobalValues(rowBB, numEntriesBB, valuesBB, indicesBB);
  }
  
  matrixBB->FillComplete (rowMapB, rowMapB);
  
  
  //EpetraBlockMatrixAssemble startup--------------------------------------------------------------
  epetraOperatorBlock opBlock(2,2);
  opBlock.setEpetraComm(epetraComm);
  
  opBlock.setRowMap(1,rowMapA);
  opBlock.setRowMap(2,rowMapB);
  
  opBlock.setColMap(1,rowMapA);
  opBlock.setColMap(2,rowMapB);
  
  opBlock.startup();
  
  opBlock.setOperator(matrixAA,1,1);
  opBlock.setOperator(matrixAB,1,2);
  opBlock.setOperator(matrixBA,2,1);
  opBlock.setOperator(matrixBB,2,2);
  
  
  //Comparison matrix------------------------------------------------------------------------------
  epetraBlockMatrixAssemble assembler(2,2);
  assembler.setEpetraComm(epetraComm);
  
  assembler.setRowMap(1,rowMapA);
  assembler.setRowMap(2,rowMapB);
  
  assembler.setColMap(1,rowMapA);
  assembler.setColMap(2,rowMapB);
  
  assembler.startup();
  
  assembler.setMatrix(matrixAA,1,1);
  assembler.setMatrix(matrixAB,1,2);
  assembler.setMatrix(matrixBA,2,1);
  assembler.setMatrix(matrixBB,2,2);
  assembler.assemble();
  
  Teuchos::RCP<Epetra_FECrsMatrix> M = assembler.getMatrix();
  
  
  //Evaluation test--------------------------------------------------------------------------------
  Epetra_MultiVector X(opBlock.OperatorDomainMap(),1);
  Epetra_MultiVector Y(opBlock.OperatorRangeMap(), 1);
  
  if(world.rank() == 0)
  {
    X.SumIntoMyValue(0,0,1.0);
    X.SumIntoMyValue(1,0,0.5);
    X.SumIntoMyValue(2,0,2.5);
    X.SumIntoMyValue(4,0,5.5);
  }
  else
  {
    X.SumIntoMyValue(0,0,2.0);
    X.SumIntoMyValue(1,0,1.5);
  }
  
  Epetra_MultiVector XX(X);
  Epetra_MultiVector YY(Y);
  
  opBlock.Apply(X,Y);
  M->Multiply(false,XX,YY);
  
  cout << "X :  " << X << endl;
  cout << "Y :  " << Y << endl;
  
  cout << "XX : " << XX << endl;
  cout << "YY : " << YY << endl;
}
