/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "tpetraSaddleAssembler.h"
#include "Tpetra_RowMatrixTransposer_decl.hpp"
#include "TpetraExt_MatrixMatrix_def.hpp"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
tpetraSaddleAssembler::
tpetraSaddleAssembler() : tpetraBlockMatrixAssemble()
{
}

tpetraSaddleAssembler::
tpetraSaddleAssembler(const UInt & Dim) : tpetraBlockMatrixAssemble(Dim,Dim)
{
}

void
tpetraSaddleAssembler::
setDim(const UInt & Dim)
{
  tpetraBlockMatrixAssemble::setDim(Dim,Dim);
}

void
tpetraSaddleAssembler::
setEpetraComm(const RCP_MPICOMM & CommDev)
{
  tpetraBlockMatrixAssemble::setTpetraComm(CommDev);
}

void
tpetraSaddleAssembler::
setEpetraComm(const COMM & CommDev)
{
  tpetraBlockMatrixAssemble::setTpetraComm(CommDev);
}



//_________________________________________________________________________________________________
// TOPOLOGY
//-------------------------------------------------------------------------------------------------
void
tpetraSaddleAssembler::
setMap(const UInt & i, const MAP & map)
{
  tpetraBlockMatrixAssemble::setRowMap(i,map);
  tpetraBlockMatrixAssemble::setColMap(i,map);
}

void
tpetraSaddleAssembler::
setMap(const UInt & i, const RCP_MAP & map)
{
  tpetraBlockMatrixAssemble::setRowMap(i,map);
  tpetraBlockMatrixAssemble::setColMap(i,map);
}

void
tpetraSaddleAssembler::
startup()
{
  tpetraBlockMatrixAssemble::startup();
}



//_________________________________________________________________________________________________
// ASSEMBLE
//-------------------------------------------------------------------------------------------------
void
tpetraSaddleAssembler::
setDiagMatrix(const MATRIX & sourceMatrix, const UInt & i)
{
  assert(sourceMatrix.getGlobalNumRows() == sourceMatrix.getGlobalNumCols());
  
  tpetraBlockMatrixAssemble::setMatrix(sourceMatrix,i,i);
}

void
tpetraSaddleAssembler::
setDiagMatrix(const RCP_MATRIX & sourceMatrix, const UInt & i)
{
  assert(sourceMatrix->getGlobalNumRows() == sourceMatrix->getGlobalNumCols());
  
  tpetraBlockMatrixAssemble::setMatrix(sourceMatrix,i,i);
}

void
tpetraSaddleAssembler::
setMTM(const MATRIX & M, const UInt & i, const Real & scale)
{
  //Assert
  assert(M.getGlobalNumCols() == (rowMaps(i)->getMaxGlobalIndex() + 1));
  
  //Transpose matrix
  //RCP_MATRIX sourceTemp = M.clone(M.getNode());
  RCP_MATRIX sourceTemp = Teuchos::RCP<MATRIX>(new MATRIX(M));
  
  Tpetra::RowMatrixTransposer<> Transpositor(sourceTemp);
  RCP_MATRIX MT = Transpositor.createTranspose();
  
  //Product
  Teuchos::RCP<const MAP> constMap(new MAP(*MT->getRowMap()));
  //RCP_MATRIX MTM = Tpetra::createCrsMatrix<>(constMap);
  RCP_MATRIX MTM = Teuchos::RCP<MATRIX>(new MATRIX(constMap,1));
  Tpetra::MatrixMatrix::Multiply(*MT,false,M,false,*MTM);
  
  //Assemble
  MTM->scale(scale);
  tpetraBlockMatrixAssemble::setMatrix(MTM,i,i);
}

void
tpetraSaddleAssembler::
setEye(const UInt & i, const Real & scale)
{
  //Build matrix
  Teuchos::RCP<const MAP> constMap(new MAP(*rowMaps(i)));
  //RCP_MATRIX Eye = Tpetra::createCrsMatrix<double,int,int>(constMap);
  RCP_MATRIX Eye = Teuchos::RCP<MATRIX>(new MATRIX(constMap,1));
  
  //Alloc
  TPETRA_SCALAR values[1];
  TPETRA_GLOBAL_ORDINAL numEntries, newRow, indices[1];
  UInt numRows = Eye->getNodeNumRows();
  
  numEntries = 1;
  values[0]  = scale;

  //Fill loop
  for(UInt row = 0; row < numRows; ++row)
  {
    newRow     = rowMaps(i)->getGlobalElement(row);
    indices[0] = rowMaps(i)->getGlobalElement(row);
    
    Eye->insertGlobalValues(newRow,numEntries,values,indices);
  }
  
  Eye->fillComplete(constMap,constMap);
  
  //Assemble
  tpetraBlockMatrixAssemble::setMatrix(Eye,i,i);
}

void
tpetraSaddleAssembler::
setSymmMatrix(const MATRIX & sourceMatrix, const UInt & i, const UInt & j)
{
  //Transpose matrix
  //RCP_MATRIX sourceTemp = sourceMatrix.clone(sourceMatrix.getNode());
  RCP_MATRIX sourceTemp = Teuchos::RCP<MATRIX>(new MATRIX(sourceMatrix));
  Teuchos::RCP<const MAP> newMap(new MAP(*rowMaps(j)));
  
  Tpetra::RowMatrixTransposer<> Transpositor(sourceTemp);
  RCP_MATRIX transposeMatrix = Transpositor.createTranspose();
  
  //Assemble
  tpetraBlockMatrixAssemble::setMatrix(sourceMatrix,    i, j);  
  tpetraBlockMatrixAssemble::setMatrix(transposeMatrix, j, i);
}

void
tpetraSaddleAssembler::
setSymmMatrix(const RCP_MATRIX & sourceMatrix, const UInt & i, const UInt & j)
{
  setSymmMatrix(*sourceMatrix,i,j);
}

void
tpetraSaddleAssembler::
assemble()
{
  tpetraBlockMatrixAssemble::assemble();
}



//_________________________________________________________________________________________________
// RETURN
//-------------------------------------------------------------------------------------------------
tpetraSaddleAssembler::RCP_MAP
tpetraSaddleAssembler::
getMap() const
{
  return(tpetraBlockMatrixAssemble::getRowMap());
}

tpetraSaddleAssembler::RCP_MATRIX
tpetraSaddleAssembler::
getMatrix() const
{
  return(tpetraBlockMatrixAssemble::getMatrix());
}
