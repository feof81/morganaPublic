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

#include "epetraMatrixManip.h"



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
  
  if(world.rank() == 0) { numRowA = 2; numRowB = 1; }
  else                  { numRowA = 2; numRowB = 1; }
  
  int elementsRowA[numRowA];
  int elementsRowB[numRowB];
  
  if(world.rank() == 0)
  {
    elementsRowA[0] = 0;
    elementsRowA[1] = 1;
    
    elementsRowB[0] = 0;
  }
  else
  {
    elementsRowA[0] = 2;
    elementsRowA[1] = 3;
    
    elementsRowB[0] = 1;
  }  
  
  Epetra_Map rowMapA(-1,numRowA,elementsRowA,0,epetraComm);
  Epetra_Map rowMapB(-1,numRowB,elementsRowB,0,epetraComm);
  
  
  //Matrix matrixA-------------------------------------------------------------
  Epetra_CrsMatrix matrixA(Copy, rowMapA, 1);
  
  int length = 4;
  double valuesA[length];
  int numEntriesA, rowA, indicesA[length];
  
  
  if(world.rank() == 0)
  {
    rowA = 0;
    numEntriesA = 2;
    indicesA[0] = 0; indicesA[1] = 2;
    valuesA[0]  = 4; valuesA[1]  = 1;
    matrixA.InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
    
    rowA = 1;
    numEntriesA = 2;
    indicesA[0] = 1; indicesA[1] = 3;
    valuesA[0]  = 2; valuesA[1]  = 0.5;
    matrixA.InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
  }
  else
  {
    rowA = 2;
    numEntriesA = 2;
    indicesA[0] = 0; indicesA[1] = 2;
    valuesA[0]  = 1; valuesA[1]  = 1;
    matrixA.InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
    
    rowA = 3;
    numEntriesA = 2;
    indicesA[0] = 1;   indicesA[1] = 3;
    valuesA[0]  = 0.5; valuesA[1]  = 4;
    matrixA.InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
  }
  
  matrixA.FillComplete(rowMapA, rowMapA);
  
  
  //Matrix matrixB-------------------------------------------------------------
  Epetra_CrsMatrix matrixB(Copy, rowMapB, 1);
  
  double valuesB[length];
  int numEntriesB, rowB, indicesB[length];
  
  
  if(world.rank() == 0)
  {
    rowB = 0;
    numEntriesB = 2;
    indicesB[0] = 0; indicesB[1] = 3;
    valuesB[0]  = 1; valuesB[1]  = -1;
    matrixB.InsertGlobalValues(rowB, numEntriesB, valuesB, indicesB);
  }
  else
  {
    rowB = 1;
    numEntriesB = 2;
    indicesB[0] = 0;  indicesB[1] = 2;
    valuesB[0]  = -1; valuesB[1]  = 1;
    matrixB.InsertGlobalValues(rowB, numEntriesB, valuesB, indicesB);
  }
  
  matrixB.FillComplete(rowMapA, rowMapB);
  

  //Manipulator
  epetraMatrixManip manip;
  //cout << manip.diag(matrixA) << endl;
  //cout << manip.diagInv(matrixA) << endl;
  //cout << manip.multiply(matrixB, matrixA, matrixA) << endl;
  cout << manip.shurDiag(matrixA, matrixB) << endl;
}

