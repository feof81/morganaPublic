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
#include "factoryPreconditioner.hpp"

#include "epetraMatrixManip.h"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"



//! Run with two processors
int main(int argc, char *argv[])
{
  //Startup--------------------------------------------------------------------
  environment env(argc,argv);
  
  communicator world;
  Epetra_MpiComm epetraComm(MPI_COMM_WORLD);
  
  assert(world.size() == 2);
  
  
  //Maps-----------------------------------------------------------------------
  int numRowA;
  
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
  Teuchos::RCP<Epetra_CrsMatrix> matrixA(new Epetra_CrsMatrix(Copy, rowMapA, 1));
  
  int length = 4;
  double valuesA[length];
  int numEntriesA, rowA, indicesA[length];
  
  
  if(world.rank() == 0)
  {
    rowA = 0;
    numEntriesA = 2;
    indicesA[0] = 0; indicesA[1] = 2;
    valuesA[0]  = 4; valuesA[1]  = 1;
    matrixA->InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
    
    rowA = 1;
    numEntriesA = 2;
    indicesA[0] = 1; indicesA[1] = 3;
    valuesA[0]  = 2; valuesA[1]  = 0.5;
    matrixA->InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
  }
  else
  {
    rowA = 2;
    numEntriesA = 2;
    indicesA[0] = 0; indicesA[1] = 2;
    valuesA[0]  = 1; valuesA[1]  = 1;
    matrixA->InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
    
    rowA = 3;
    numEntriesA = 2;
    indicesA[0] = 1;   indicesA[1] = 3;
    valuesA[0]  = 0.5; valuesA[1]  = 4;
    matrixA->InsertGlobalValues(rowA, numEntriesA, valuesA, indicesA);
  }
  
  matrixA->FillComplete(rowMapA, rowMapA);
  
  
  //Preconditioner-------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> preconList(new Teuchos::ParameterList);
  preconList->set("precClass",precMl);
  preconList->sublist("mlList").set("max levels",2);
  preconList->sublist("mlList").set("increasing or decreasing","increasing");
  preconList->sublist("mlList").set("aggregation: type", "MIS");
  preconList->sublist("mlList").set("aggregation: nodes per aggregate", 2);
  preconList->sublist("mlList").set("smoother: pre or post", "both");
  preconList->sublist("mlList").set("coarse: type","Amesos-KLU");
  preconList->sublist("mlList").set("smoother: type", "IFPACK");
  preconList->sublist("mlList").set("smoother: ifpack level-of-fill", 2.0);
  preconList->sublist("mlList").set("smoother: ifpack overlap", 1);
  
  typedef factoryPreconditioner<double, Epetra_MultiVector, Epetra_RowMatrix, Belos::EpetraPrecOp> FACTORYPREC;
  typedef virtualPreconditioner<double, Epetra_MultiVector, Epetra_RowMatrix, Belos::EpetraPrecOp> PREC;
  
  FACTORYPREC factoryPrec;
  Teuchos::RCP<PREC> preconditioner = factoryPrec.create(preconList);
  preconditioner->setOperator(matrixA);
  preconditioner->initialize();
  preconditioner->compute();
  Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = preconditioner->getPreconditioner();
}

