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
#include "Teuchos_Comm.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map_decl.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"

#include "pMapItem.h"
#include "pMapItemShare.h"
#include "pVect.hpp"
#include "pVectManip.hpp"
#include "pVectComm.hpp"
#include "pVectGlobalManip.hpp"

#include "tpetraSaddleAssembler.h"

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;

//! Run with two processors
int main(int argc, char *argv[])
{
  //Typedefs-------------------------------------------------------------------
  typedef tpetraBlockMatrixAssemble::MAP               MAP;
  typedef tpetraBlockMatrixAssemble::RCP_MAP       RCP_MAP;
  typedef tpetraBlockMatrixAssemble::MATRIX         MATRIX;
  typedef tpetraBlockMatrixAssemble::RCP_MATRIX RCP_MATRIX;
    
  typedef tpetraBlockMatrixAssemble::RCP_MAP                RCP_MAP;
  typedef tpetraBlockMatrixAssemble::TPETRA_GLOBAL_TYPE     TPETRA_GLOBAL_TYPE;
  typedef tpetraBlockMatrixAssemble::TPETRA_GLOBAL_ORDINAL  TPETRA_GLOBAL_ORDINAL;
  typedef tpetraBlockMatrixAssemble::TPETRA_SCALAR          TPETRA_SCALAR;  
    
  //Startup--------------------------------------------------------------------
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  
  Teuchos::RCP<const Comm<int> > teuchosComm = Teuchos::rcp(new MpiComm<int> (world));
  
  //Maps-----------------------------------------------------------------------
  TPETRA_GLOBAL_ORDINAL numRowA, numRowB;
  TPETRA_GLOBAL_ORDINAL numColA, numColB;
  
  if(world.rank() == 0) { numRowA = 2; numRowB = 1; }
  else                  { numRowA = 2; numRowB = 1; }
  
  if(world.rank() == 0) { numColA = 3; numColB = 4; }
  else                  { numColA = 3; numColB = 4; }
  
  TPETRA_GLOBAL_ORDINAL elementsRowA[numRowA];
  TPETRA_GLOBAL_ORDINAL elementsRowB[numRowB];
  TPETRA_GLOBAL_ORDINAL elementsColA[numColA];
  TPETRA_GLOBAL_ORDINAL elementsColB[numColB];
  
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
  
  TPETRA_GLOBAL_TYPE globSize = Teuchos::OrdinalTraits<TPETRA_GLOBAL_TYPE>::invalid();
  RCP_MAP rowMapA(new MAP(globSize,elementsRowA,numRowA,0,teuchosComm));
  RCP_MAP rowMapB(new MAP(globSize,elementsRowB,numRowB,0,teuchosComm));
  RCP_MAP colMapA(new MAP(globSize,elementsColA,numColA,0,teuchosComm));
  RCP_MAP colMapB(new MAP(globSize,elementsColB,numColB,0,teuchosComm));

  
  //Matrix matrixA-------------------------------------------------------------
  RCP_MATRIX matrixA = Teuchos::rcp(new MATRIX(rowMapA,3));
  
  TPETRA_GLOBAL_ORDINAL length = 4;
  TPETRA_SCALAR valuesA[length];
  TPETRA_GLOBAL_ORDINAL numEntriesA, rowA, indicesA[length];
  
  if(world.rank() == 0)
  {
    rowA = 0;
    numEntriesA = 3;
    indicesA[0] = 0; indicesA[1] = 1; indicesA[2] = 2;
    valuesA[0]  = 6; valuesA[1]  = 4; valuesA[2]  = 6;
    matrixA->insertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
    
    rowA = 1;
    numEntriesA = 3;
    indicesA[0] = 0; indicesA[1] = 1; indicesA[2] = 2;
    valuesA[0]  = 3; valuesA[1]  = 2; valuesA[2]  = 3;
    matrixA->insertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
  }
  else
  {
    rowA = 2;
    numEntriesA = 3;
    indicesA[0] = 1; indicesA[1] = 2; indicesA[2] = 3;
    valuesA[0]  = 1; valuesA[1]  = 2; valuesA[2]  = 1;
    matrixA->insertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
    
    rowA = 3;
    numEntriesA = 3;
    indicesA[0] = 1; indicesA[1] = 2; indicesA[2] = 3;
    valuesA[0]  = 2; valuesA[1]  = 7; valuesA[2]  = 2;
    matrixA->insertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
  }
  
  matrixA->fillComplete(rowMapA,rowMapA);
  
  
  //Matrix matrixB-------------------------------------------------------------
  RCP_MATRIX matrixB = Teuchos::rcp(new MATRIX(rowMapB,4));
  
  TPETRA_SCALAR valuesB[length];
  TPETRA_GLOBAL_ORDINAL numEntriesB, rowB, indicesB[length];
  
  if(world.rank() == 0)
  {
    rowB = 0;
    numEntriesB = 4;
    indicesB[0] = 0; indicesB[1] = 1; indicesB[2] = 2; indicesB[3] = 3;
    valuesB[0]  = 1; valuesB[1]  = 2; valuesB[2]  = 3; valuesB[3]  = 4;
    matrixB->insertGlobalValues(rowB, numEntriesB, valuesB, indicesB);
  }
  else
  {
    rowB = 1;
    numEntriesB = 4;
    indicesB[0] = 0; indicesB[1] = 1; indicesB[2] = 2; indicesB[3] = 3;
    valuesB[0]  = 5; valuesB[1]  = 6; valuesB[2]  = 7; valuesB[3]  = 8;
    matrixB->insertGlobalValues(rowB, numEntriesB, valuesB, indicesB);
  }
  
  matrixB->fillComplete(rowMapA,rowMapB);
  
  
  // SaddleAssembler ----------------------------------------------------------
  tpetraSaddleAssembler assembler(2);
  assembler.setTpetraComm(teuchosComm);
  
  assembler.setMap(1,rowMapA);
  assembler.setMap(2,rowMapB);
  assembler.startup();
  
  assembler.setDiagMatrix(matrixA,1);
  assembler.setSymmMatrix(matrixB,2,1);
  assembler.assemble();
  
  auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  
  RCP_MATRIX newMatrix = assembler.getMatrix();
  newMatrix->describe(*out,Teuchos::VERB_EXTREME);
  
  cout << *assembler.getMap() << endl;
}
