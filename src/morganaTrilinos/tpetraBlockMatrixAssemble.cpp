/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "tpetraBlockMatrixAssemble.h"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_RowMatrixTransposer_decl.hpp"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
tpetraBlockMatrixAssemble::
tpetraBlockMatrixAssemble()
{
  dimLoaded     = false;
  startupOk     = false;
  assembleOK    = false;
  commDevLoaded = false;
}

tpetraBlockMatrixAssemble::
tpetraBlockMatrixAssemble(const UInt & NumRows, const UInt & NumCols)
{
  dimLoaded     = true;
  startupOk     = false;
  assembleOK    = false;
  commDevLoaded = false;
  
  numRows = NumRows;
  numCols = NumCols;
  
  rowLoaded.resize(NumRows);
  colLoaded.resize(NumCols);
  
  rowMaps.resize(NumRows);
  colMaps.resize(NumCols);
  
  for(UInt i=1; i <= NumRows; ++i)
  { rowLoaded(i) = false; }
  
  for(UInt i=1; i <= NumCols; ++i)
  { colLoaded(i) = false; }
}

void
tpetraBlockMatrixAssemble::
setDim(const UInt & NumRows, const UInt & NumCols)
{
  dimLoaded  = true;
  startupOk  = false;
  assembleOK = false;
  
  numRows = NumRows;
  numCols = NumCols;
  
  rowLoaded.resize(NumRows);
  colLoaded.resize(NumCols);
  
  rowMaps.resize(NumRows);
  colMaps.resize(NumCols);
  
  for(UInt i=1; i <= NumRows; ++i)
  { rowLoaded(i) = false; }
  
  for(UInt i=1; i <= NumCols; ++i)
  { colLoaded(i) = false; }
}

void
tpetraBlockMatrixAssemble::
setTpetraComm(const RCP_MPICOMM & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

void
tpetraBlockMatrixAssemble::
setTpetraComm(const COMM & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcp(new COMM(CommDev));
}


//_________________________________________________________________________________________________
// TOPOLOGY
//-------------------------------------------------------------------------------------------------
void
tpetraBlockMatrixAssemble::
setRowMap(const UInt & i, const MAP & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numRows);
  
  rowLoaded(i) = true;
  rowMaps(i)   = Teuchos::rcp(new MAP(map));
  
  startupOk  = false;
  assembleOK = false;
}

void
tpetraBlockMatrixAssemble::
setRowMap(const UInt & i, const RCP_MAP & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numRows);
  
  rowLoaded(i) = true;
  rowMaps(i)   = map;
  
  startupOk  = false;
  assembleOK = false;
}

void
tpetraBlockMatrixAssemble::
setColMap(const UInt & i, const MAP & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numCols);
  
  colLoaded(i) = true;
  colMaps(i)   = Teuchos::rcp(new MAP(map));
  
  startupOk  = false;
  assembleOK = false;
}

void
tpetraBlockMatrixAssemble::
setColMap(const UInt & i, const RCP_MAP & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numCols);
  
  colLoaded(i) = true;
  colMaps(i)   = map;
  
  startupOk  = false;
  assembleOK = false;
}

void
tpetraBlockMatrixAssemble::
startup()
{
  //Checking---------------------------------------------------------------------------------------
  assert(commDevLoaded);  
  assert(dimLoaded);
  
  for(UInt i=1; i <= numRows; ++i)
  { assert(rowLoaded(i)); }
  
  for(UInt i=1; i <= numCols; ++i)
  { assert(colLoaded(i)); }
  
  //Max gids---------------------------------------------------------------------------------------
  offsetRow.resize(numRows);
  offsetCol.resize(numCols);
  
  if(numRows != 0)
  { offsetRow(1) = 0; }
  
  for(UInt i=2; i <= numRows; ++i)
  { offsetRow(i) = offsetRow(i-1) + (rowMaps(i-1)->getMaxAllGlobalIndex()  + 1); }
  
  if(numCols != 0)
  { offsetCol(1) = 0; }
  
  for(UInt i=2; i <= numCols; ++i)
  { offsetCol(i) = offsetCol(i-1) + (colMaps(i-1)->getMaxAllGlobalIndex()  + 1); }
  
  //Build row--------------------------------------------------------------------------------------
  UInt sizeRow = 0;
  
  for(UInt i=1; i <= numRows; ++i)
  { sizeRow += rowMaps(i)->getNodeNumElements() ; }
  
  Teuchos::ArrayView<const TPETRA_GLOBAL_ORDINAL> rowSource;
  sVect<TPETRA_GLOBAL_ORDINAL> rowElements(sizeRow);
  
  UInt prgr = 1;
  
  for(UInt i=1; i <= numRows; ++i)
  {
    rowSource = rowMaps(i)->getNodeElementList();
    
    for(UInt k=0; k < rowSource.size(); ++k)
    {
      rowElements(prgr) = rowSource[k] + offsetRow(i);
      assert(prgr <= sizeRow);
      prgr++;
    }
  }
  
  Teuchos::ArrayView<TPETRA_GLOBAL_ORDINAL> rowTeuchos(rowElements);
  TPETRA_GLOBAL_TYPE globSize = Teuchos::OrdinalTraits<TPETRA_GLOBAL_TYPE>::invalid();
  outRowMap = Teuchos::rcp(new Tpetra::Map<>(globSize,rowTeuchos,0,commDev));
  
  //Build col--------------------------------------------------------------------------------------
  UInt sizeCol = 0;
  
  for(UInt i=1; i <= numCols; ++i)
  { sizeCol += colMaps(i)->getNodeNumElements() ; }
  
  Teuchos::ArrayView<const TPETRA_GLOBAL_ORDINAL> colSource;
  sVect<TPETRA_GLOBAL_ORDINAL> colElements(sizeCol);
  
  prgr = 1;
  
  for(UInt i=1; i <= numCols; ++i)
  {
    colSource = colMaps(i)->getNodeElementList();
    
    for(UInt k=0; k < colSource.size(); ++k)
    {
      colElements(prgr) = colSource[k] + offsetCol(i);
      assert(prgr <= sizeCol);
      prgr++;
    }
  }
  
  Teuchos::ArrayView<TPETRA_GLOBAL_ORDINAL> colTeuchos(colElements);
  outColMap = Teuchos::rcp(new Tpetra::Map<>(globSize,colTeuchos,0,commDev));
  
  //Alloc the final matrix-------------------------------------------------------------------------
  //Teuchos::RCP<const MAP> constMap(new MAP(*outRowMap));
  //outMatrix = Teuchos::RCP<MATRIX>(new MATRIX(constMap,10));
  
  //Resize the maximum number----------------------------------------------------------------------
  maxElements.resize(numCols);
  
  for(UInt i=1; i <= maxElements.size(); ++i)
  { maxElements(i) = 0; }
  
  matrixList.clear();
  Ilist.clear();
  Jlist.clear();
  
  //Flag-------------------------------------------------------------------------------------------
  startupOk  = true;
  assembleOK = false;
}


//_________________________________________________________________________________________________
// ASSEMBLE
//-------------------------------------------------------------------------------------------------
void
tpetraBlockMatrixAssemble::
setMatrix(const MATRIX & sourceMatrix, const UInt & i, const UInt & j)
{
  RCP_MATRIX rcpMatrix(new MATRIX(sourceMatrix));
  setMatrix(rcpMatrix,i,j);
    
  /*
  //Asserts
  assert(startupOk);
  assert(i >= 1); assert(i <= numRows);
  assert(i >= 1); assert(j <= numCols);
  assert(rowMaps(i)->getNodeNumElements() == sourceMatrix.getNodeNumRows());
  assert(rowMaps(i)->isSameAs(*sourceMatrix.getRowMap()));
  assert(rowMaps(i)->isSameAs(*sourceMatrix.getRangeMap()));
  assert(colMaps(j)->isSameAs(*sourceMatrix.getDomainMap()));
  
  //Estract data
  UInt    length = sourceMatrix.getGlobalMaxNumRowEntries() + 1;
  UInt rowSource = sourceMatrix.getNodeNumRows();
  
  //Alloc
  std::vector<TPETRA_SCALAR>                valuesInit(length);
  std::vector<TPETRA_GLOBAL_ORDINAL>        indicesInit(length);
  Teuchos::ArrayView<TPETRA_SCALAR>         values(valuesInit);
  Teuchos::ArrayView<TPETRA_GLOBAL_ORDINAL> indices(indicesInit);
  
  TPETRA_SCALAR         * valuesNew  = new TPETRA_SCALAR[length];
  TPETRA_GLOBAL_ORDINAL * indicesNew = new TPETRA_GLOBAL_ORDINAL[length];
  
  std::size_t numEntries;
  TPETRA_GLOBAL_ORDINAL tpetraRow, newRow;
  
  
  //Copy cycle
  for(UInt row = 0; row < rowSource; ++row)
  {
    //Estrazione
    tpetraRow = rowMaps(i)->getGlobalElement(row);
    sourceMatrix.getGlobalRowCopy(tpetraRow,indices,values,numEntries);
    
    //Aggiornamento indici
    newRow = rowMaps(i)->getGlobalElement(row) + offsetRow(i);
    
    for(int k=0; k < numEntries; ++k)
    {
      assert(indices[k] < int(sourceMatrix.getGlobalNumCols()));
      indicesNew[k] = indices[k] + offsetCol(j);
      valuesNew[k]  = values[k];
    }
    
    //Memorizzazione 
    outMatrix->insertGlobalValues(newRow,
                                  TPETRA_GLOBAL_ORDINAL(numEntries),
                                  valuesNew,
                                  indicesNew);
  }
  
  delete[] valuesNew;
  delete[] indicesNew;
  assembleOK = false;*/
}

void
tpetraBlockMatrixAssemble::
setMatrix(const RCP_MATRIX & sourceMatrix, const UInt & i, const UInt & j)
{
  matrixList.push_back(sourceMatrix);
  Ilist.push_back(i);
  Jlist.push_back(j);
  maxElements(j) = std::max(maxElements(j),
                            UInt(sourceMatrix->getGlobalMaxNumRowEntries()));    
    
  //setMatrix(*sourceMatrix,i,j);
}

void
tpetraBlockMatrixAssemble::
setMatrixTransp(const MATRIX & sourceMatrix, const UInt & i, const UInt & j)
{
  //RCP_MATRIX matrixRcp = sourceMatrix.clone(sourceMatrix.getNode());  
  RCP_MATRIX matrixRcp = Teuchos::RCP<MATRIX>(new MATRIX(sourceMatrix));
  setMatrixTransp(matrixRcp,i,j);
}

void
tpetraBlockMatrixAssemble::
setMatrixTransp(const RCP_MATRIX & sourceMatrix, const UInt & i, const UInt & j)
{
  Tpetra::RowMatrixTransposer<> transposer(sourceMatrix);
  RCP_MATRIX tempMat = transposer.createTranspose();
  setMatrix(tempMat,i,j);
}

void
tpetraBlockMatrixAssemble::
addMatrix(const RCP_MATRIX & sourceMatrix, const UInt & i, const UInt & j)
{
  //Asserts----------------------------------------------------------------------------------------
  assert(startupOk);
  assert(i >= 1); assert(i <= numRows);
  assert(i >= 1); assert(j <= numCols);
  assert(rowMaps(i)->getNodeNumElements() == sourceMatrix->getNodeNumRows());
  assert(rowMaps(i)->isSameAs(*sourceMatrix->getRowMap()));
  assert(rowMaps(i)->isSameAs(*sourceMatrix->getRangeMap()));
  assert(colMaps(j)->isSameAs(*sourceMatrix->getDomainMap()));
  
  //Estract data-----------------------------------------------------------------------------------
  UInt    length = sourceMatrix->getGlobalMaxNumRowEntries() + 1;
  UInt rowSource = sourceMatrix->getNodeNumRows();
  
  //Alloc------------------------------------------------------------------------------------------
  std::vector<TPETRA_SCALAR>                valuesInit(length);
  std::vector<TPETRA_GLOBAL_ORDINAL>        indicesInit(length);
  Teuchos::ArrayView<TPETRA_SCALAR>         values(valuesInit);
  Teuchos::ArrayView<TPETRA_GLOBAL_ORDINAL> indices(indicesInit);
  
  TPETRA_SCALAR         * valuesNew  = new TPETRA_SCALAR[length];
  TPETRA_GLOBAL_ORDINAL * indicesNew = new TPETRA_GLOBAL_ORDINAL[length];
  
  std::size_t numEntries;
  TPETRA_GLOBAL_ORDINAL tpetraRow, newRow;
  
  //Copy cycle-------------------------------------------------------------------------------------
  for(UInt row = 0; row < rowSource; ++row)
  {
    //Estrazione
    tpetraRow = rowMaps(i)->getGlobalElement(row);
    sourceMatrix->getGlobalRowCopy(tpetraRow,indices,values,numEntries);
    
    //Aggiornamento indici
    newRow = rowMaps(i)->getGlobalElement(row) + offsetRow(i);
    
    for(int k=0; k < numEntries; ++k)
    {
      assert(indices[k] < int(sourceMatrix->getGlobalNumCols()));
      indicesNew[k] = indices[k] + offsetCol(j);
      valuesNew[k]  = values[k];
    }
    
    //Memorizzazione 
    outMatrix->insertGlobalValues(newRow,
                                  TPETRA_GLOBAL_ORDINAL(numEntries),
                                  valuesNew,
                                  indicesNew);
  }
  
  delete[] valuesNew;
  delete[] indicesNew;
  assembleOK = false;
}

void
tpetraBlockMatrixAssemble::
assemble()
{
  //Alloc the matrix-------------------------------------------------------------------------------
  UInt maxEl = 0;
  for(UInt i=1; i <= maxElements.size(); ++i)
  { maxEl += maxElements(i); }
  
  Teuchos::RCP<const MAP> constMap(new MAP(*outRowMap));
  outMatrix = Teuchos::RCP<MATRIX>(new MATRIX(constMap,maxEl));
    
  //Copy the matrices------------------------------------------------------------------------------
  assert(matrixList.size() == Ilist.size());
  assert(matrixList.size() == Jlist.size());
  
  for(UInt k=1; k <= matrixList.size(); ++k)
  {
    addMatrix(matrixList(k),
              Ilist(k), 
              Jlist(k));
  }    
    
  //Final assemble---------------------------------------------------------------------------------
  outMatrix->fillComplete(outColMap,outRowMap);
  assembleOK = true;
}

tpetraBlockMatrixAssemble::RCP_MAP
tpetraBlockMatrixAssemble::
getRowMap() const
{
  assert(startupOk);
  return(outRowMap);
}

tpetraBlockMatrixAssemble::RCP_MAP
tpetraBlockMatrixAssemble::
getColMap() const
{
  assert(startupOk);
  return(outColMap);
}

tpetraBlockMatrixAssemble::RCP_MATRIX
tpetraBlockMatrixAssemble::
getMatrix() const
{
  assert(assembleOK);
  return(outMatrix);
}

void
tpetraBlockMatrixAssemble::
printToFile(const string & name) const
{
  Tpetra::MatrixMarket::Writer<MATRIX>::writeSparseFile(name,outMatrix);
  assert(assembleOK);
}
