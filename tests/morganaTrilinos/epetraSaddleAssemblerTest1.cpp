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
#include "epetraSaddleAssembler.h"



//! Run with two processors
int main(int argc, char *argv[])
{
  //Startup--------------------------------------------------------------------
  environment env(argc,argv);
  
  communicator world;
  Epetra_MpiComm epetraComm(MPI_COMM_WORLD);
  
  assert(world.size() == 2);
  
  
  //Maps-----------------------------------------------------------------------
  int numRowA, numRowB;
  int numColA, numColB;
  
  if(world.rank() == 0) { numRowA = 2; numRowB = 1; }
  else                  { numRowA = 2; numRowB = 1; }
  
  if(world.rank() == 0) { numColA = 3; numColB = 4; }
  else                  { numColA = 3; numColB = 4; }
  
  int elementsRowA[numRowA];
  int elementsRowB[numRowB];
  int elementsColA[numColA];
  int elementsColB[numColB];
  
  if(world.rank() == 0)
  {
    elementsRowA[0] = 0;
    elementsRowA[1] = 1;
    
    elementsRowB[0] = 0;
    
    elementsColA[0] = 0;
    elementsColA[1] = 1;
    elementsColA[2] = 2;
    
    elementsColB[0] = 0;
    elementsColB[1] = 1;
    elementsColB[2] = 2;
    elementsColB[3] = 3;
  }
  else
  {
    elementsRowA[0] = 2;
    elementsRowA[1] = 3;
    
    elementsRowB[0] = 1;
    
    elementsColA[0] = 1;
    elementsColA[1] = 2;
    elementsColA[2] = 3;
    
    elementsColB[0] = 0;
    elementsColB[1] = 1;
    elementsColB[2] = 2;
    elementsColB[3] = 3;
  }  
  
  Epetra_Map rowMapA(-1,numRowA,elementsRowA,0,epetraComm);
  Epetra_Map rowMapB(-1,numRowB,elementsRowB,0,epetraComm);
  
  Epetra_Map colMapA(-1,numColA,elementsColA,0,epetraComm);
  Epetra_Map colMapB(-1,numColB,elementsColB,0,epetraComm);
  
  
  //Matrix matrixA-------------------------------------------------------------
  Epetra_CrsMatrix matrixA(Copy, rowMapA, colMapA, 1);
  
  int length = 4;
  double valuesA[length];
  int numEntriesA, rowA, indicesA[length];
  
  
  if(world.rank() == 0)
  {
    rowA = 0;
    numEntriesA = 3;
    indicesA[0] = 0; indicesA[1] = 1; indicesA[2] = 2;
    valuesA[0]  = 6; valuesA[1]  = 4; valuesA[2]  = 6;
    matrixA.InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
    
    rowA = 1;
    numEntriesA = 3;
    indicesA[0] = 0; indicesA[1] = 1; indicesA[2] = 2;
    valuesA[0]  = 3; valuesA[1]  = 2; valuesA[2]  = 3;
    matrixA.InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
  }
  else
  {
    rowA = 2;
    numEntriesA = 3;
    indicesA[0] = 1; indicesA[1] = 2; indicesA[2] = 3;
    valuesA[0]  = 1; valuesA[1]  = 2; valuesA[2]  = 1;
    matrixA.InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
    
    rowA = 3;
    numEntriesA = 3;
    indicesA[0] = 1; indicesA[1] = 2; indicesA[2] = 3;
    valuesA[0]  = 2; valuesA[1]  = 7; valuesA[2]  = 2;
    matrixA.InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
  }
  
  matrixA.FillComplete(rowMapA, rowMapA);
  
  
  //Matrix matrixB-------------------------------------------------------------
  Epetra_CrsMatrix matrixB(Copy, rowMapB, colMapB, 1);
  
  double valuesB[length];
  int numEntriesB, rowB, indicesB[length];
  
  
  if(world.rank() == 0)
  {
    rowB = 0;
    numEntriesB = 4;
    indicesB[0] = 0; indicesB[1] = 1; indicesB[2] = 2; indicesB[3] = 3;
    valuesB[0]  = 1; valuesB[1]  = 2; valuesB[2]  = 3; valuesB[3]  = 4;
    matrixB.InsertGlobalValues(rowB, numEntriesB, valuesB, indicesB);
  }
  else
  {
    rowB = 1;
    numEntriesB = 4;
    indicesB[0] = 0; indicesB[1] = 1; indicesB[2] = 2; indicesB[3] = 3;
    valuesB[0]  = 5; valuesB[1]  = 6; valuesB[2]  = 7; valuesB[3]  = 8;
    matrixB.InsertGlobalValues(rowB, numEntriesB, valuesB, indicesB);
  }
  
  matrixB.FillComplete(rowMapA,rowMapB);
  
  
  
  // SaddleAssembler ----------------------------------------------------------
  epetraSaddleAssembler assembler(2);
  assembler.setEpetraComm(epetraComm);
  
  assembler.setMap(1,rowMapA);
  assembler.setMap(2,rowMapB);
  assembler.startup();
  
  assembler.setDiagMatrix(matrixA,1);
  assembler.setSymmMatrix(matrixB,2,1);
  assembler.assemble();
  
  cout << *assembler.getMap() << endl;
  cout << *assembler.getMatrix() << endl;
}

