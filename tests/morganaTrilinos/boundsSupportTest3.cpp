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
#include "boundsSupport.h"
#include "epetraBlockVectorAssemble.h"



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
  
  if(world.rank() == 0) { numRowA = 2; }
  else                  { numRowA = 2; }
  
  int elementsRowA[numRowA];
  
  if(world.rank() == 0)
  {
    elementsRowA[0] = 0;
    elementsRowA[1] = 1;
  }
  else
  {
    elementsRowA[0] = 2;
    elementsRowA[1] = 3;
  }  
  
  Epetra_Map rowMapA(-1,numRowA,elementsRowA,0,epetraComm);

  
  //Matrix matrixA-------------------------------------------------------------
  Epetra_CrsMatrix matrixA(Copy, rowMapA, 1);
  
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
  
  
  //Vectors----------------------------------------------------------
  Epetra_MultiVector sourceA(rowMapA,1);
  Epetra_MultiVector vectorV(rowMapA,1);
  
  int I;
  double tot;
  
  if(world.rank() == 0)
  {
    I = 0; tot = 4;
    vectorV.SumIntoGlobalValue(I,0,tot);
    
    I = 1; tot = 3;
    vectorV.SumIntoGlobalValue(I,0,tot);
    
    
    I = 0; tot = -1;
    sourceA.SumIntoGlobalValue(I,0,tot);
    
    I = 1; tot = -2;
    sourceA.SumIntoGlobalValue(I,0,tot);
  }
  else
  {
    I = 2; tot = 2;
    vectorV.SumIntoGlobalValue(I,0,tot);
  
    I = 3; tot = 1;
    vectorV.SumIntoGlobalValue(I,0,tot);
    
    
    I = 2; tot = -3;
    sourceA.SumIntoGlobalValue(I,0,tot);
  
    I = 3; tot = -4;
    sourceA.SumIntoGlobalValue(I,0,tot);
  }
  
  
  //Bounds-----------------------------------------------------------
  boundsSupport bounds(epetraComm);
  Epetra_Map         rowMapV = bounds.getMapRef();
  Epetra_CrsMatrix   matrixV = bounds.vector_to_matrix(vectorV);
  Epetra_MultiVector sourceV = bounds.getSourceRef(8.0);
  
  
  //SaddleAssembler -------------------------------------------------
  epetraSaddleAssembler assembler(2);
  assembler.setEpetraComm(epetraComm);
  
  assembler.setMap(1,rowMapA);
  assembler.setMap(2,rowMapV);
  assembler.startup();
 
  assembler.setDiagMatrix(matrixA,1);
  assembler.setSymmMatrix(matrixV,1,2);

  assembler.assemble();
  
  
  //Vector assembler-------------------------------------------------
  epetraBlockVectorAssemble vectass(2);
  vectass.setEpetraComm(epetraComm);
  
  vectass.setRowMap(1,rowMapA);
  vectass.setRowMap(2,rowMapV);
  vectass.startup();
  
  vectass.setVector(sourceA,1);
  vectass.setVector(sourceV,2);
  
  //cout << *assembler.getMap() << endl;
  //cout << *assembler.getMatrix() << endl;
  cout << *vectass.getVector() << endl;
}

