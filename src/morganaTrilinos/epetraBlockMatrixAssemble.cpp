#include "epetraBlockMatrixAssemble.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
epetraBlockMatrixAssemble::
epetraBlockMatrixAssemble()
{
  dimLoaded     = false;
  startupOk     = false;
  assembleOK    = false;
  commDevLoaded = false;
}

epetraBlockMatrixAssemble::
epetraBlockMatrixAssemble(const UInt & NumRows, const UInt & NumCols)
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
epetraBlockMatrixAssemble::
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
epetraBlockMatrixAssemble::
setEpetraComm(const RCP_MPICOMM & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

void
epetraBlockMatrixAssemble::
setEpetraComm(const Epetra_MpiComm & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcp(new Epetra_MpiComm(CommDev));
}

    

//_________________________________________________________________________________________________
// TOPOLOGY
//-------------------------------------------------------------------------------------------------
void
epetraBlockMatrixAssemble::
setRowMap(const UInt & i, const Epetra_BlockMap & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numRows);
  assert(map.UniqueGIDs());
  
  rowLoaded(i) = true;
  rowMaps(i)   = Teuchos::rcp(new Epetra_BlockMap(map));
  
  startupOk  = false;
  assembleOK = false;
}

void
epetraBlockMatrixAssemble::
setRowMap(const UInt & i, const RCP_BLOCKMAP & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numRows);
  assert(map->UniqueGIDs());
  
  rowLoaded(i) = true;
  rowMaps(i)   = map;
  
  startupOk  = false;
  assembleOK = false;
}

void
epetraBlockMatrixAssemble::
setColMap(const UInt & i, const Epetra_BlockMap & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numCols);
  
  colLoaded(i) = true;
  colMaps(i)   = Teuchos::rcp(new Epetra_BlockMap(map));
  
  startupOk  = false;
  assembleOK = false;
}

void
epetraBlockMatrixAssemble::
setColMap(const UInt & i, const RCP_BLOCKMAP & map)
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
epetraBlockMatrixAssemble::
startup()
{
  //Typedefs
  typedef Teuchos::RCP<int>  RCP_INT;
  
  //Checking
  assert(commDevLoaded);  
  assert(dimLoaded);
  
  for(UInt i=1; i <= numRows; ++i)
  { assert(rowLoaded(i)); }
  
  for(UInt i=1; i <= numCols; ++i)
  { assert(colLoaded(i)); }
  
  
  //Max gids
  offsetRow.resize(numRows);
  offsetCol.resize(numCols);
  
  if(numRows != 0)
  { offsetRow(1) = 0; }
  
  for(UInt i=2; i <= numRows; ++i)
  { offsetRow(i) = offsetRow(i-1) + (rowMaps(i-1)->MaxAllGID() + 1); }
  
  if(numCols != 0)
  { offsetCol(1) = 0; }
  
  for(UInt i=2; i <= numCols; ++i)
  { offsetCol(i) = offsetCol(i-1) + (colMaps(i-1)->MaxAllGID() + 1); }
  
  
  //Build row
  UInt sizeRow = 0;
  
  for(UInt i=1; i <= numRows; ++i)
  { sizeRow += rowMaps(i)->NumMyElements(); }
  
  int   sizeSourceRow;
  int * rowSource;
  sVect<int> rowElements(sizeRow);
  
  UInt prgr = 1;
  
  for(UInt i=1; i <= numRows; ++i)
  {
    sizeSourceRow = rowMaps(i)->NumMyElements();
    rowSource     = rowMaps(i)->MyGlobalElements();
    
    for(UInt k=0; k < UInt(sizeSourceRow); ++k)
    {
      rowElements(prgr) = rowSource[k] + offsetRow(i);
      assert(prgr <= sizeRow);
      prgr++;
    }
  }
  
  outRowMap = Teuchos::rcp(new Epetra_Map(-1, int(sizeRow), &rowElements[0], 0, *commDev) );
  assert(outRowMap->UniqueGIDs());
  
  
  //Build col
  UInt sizeCol = 0;
  
  for(UInt i=1; i <= numCols; ++i)
  { sizeCol += colMaps(i)->NumMyElements(); }
  
  int   sizeSourceCol;
  int * colSource;
  sVect<int> colElements(sizeCol);
  
  prgr = 1;
  
  for(UInt i=1; i <= numCols; ++i)
  {
    sizeSourceCol = colMaps(i)->NumMyElements();
    colSource     = colMaps(i)->MyGlobalElements();
    
    for(UInt k=0; k < UInt(sizeSourceCol); ++k)
    {
      colElements(prgr) = colSource[k] + offsetCol(i);
      assert(prgr <= sizeCol);
      prgr++;
    }
  }
  
  outColMap = Teuchos::rcp(new Epetra_Map(-1, int(sizeCol), &colElements[0], 0, *commDev) );
  
  //Alloc the final matrix
  outMatrix = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, *outRowMap, 1) );
  
  //Flag
  startupOk  = true;
  assembleOK = false;
}


//_________________________________________________________________________________________________
// ASSEMBLE
//-------------------------------------------------------------------------------------------------
void
epetraBlockMatrixAssemble::
setMatrix(const Epetra_CrsMatrix & sourceMatrix, const UInt & i, const UInt & j)
{
  //Asserts
  assert(startupOk);
  assert(i >= 1); assert(i <= numRows);
  assert(i >= 1); assert(j <= numCols);
  assert(rowMaps(i)->NumMyElements() == sourceMatrix.NumMyRows());
  assert(rowMaps(i)->SameAs(sourceMatrix.RowMap()));
  assert(rowMaps(i)->SameAs(sourceMatrix.RangeMap()));
  assert(colMaps(j)->SameAs(sourceMatrix.DomainMap()));
  
  //Estract data
  UInt    length = sourceMatrix.MaxNumEntries() + 1;
  UInt rowSource = sourceMatrix.NumMyRows();
  
  //Alloc
  double * values  = new double[length];
  int    * indices = new int[length]; 
  int numEntries, newRow;
  int ret;
  
  //Copy cycle
  for(UInt row = 0; row < rowSource; ++row)
  {
    //Estrazione
    ret = sourceMatrix.ExtractGlobalRowCopy(rowMaps(i)->GID(row), length, numEntries, values, indices);
    assert(ret == 0);
    
    //Aggiornamento indici
    newRow = rowMaps(i)->GID(row) + offsetRow(i);    
    
    for(int k=0; k < numEntries; ++k)
    {
      assert(indices[k] < int(sourceMatrix.NumGlobalCols()));
      indices[k] += offsetCol(j);
    }
    
    //Memorizzazione 
    outMatrix->InsertGlobalValues(newRow, numEntries, values, indices);
  }
  
  delete[] values;
  delete[] indices;
  assembleOK = false;
}

void
epetraBlockMatrixAssemble::
setMatrix(const RCP_MATRIX & sourceMatrix, const UInt & i, const UInt & j)
{
  //Asserts
  assert(startupOk);
  assert(i >= 1); assert(i <= numRows);
  assert(i >= 1); assert(j <= numCols);
  assert(rowMaps(i)->NumMyElements() == sourceMatrix->NumMyRows());
  assert(rowMaps(i)->SameAs(sourceMatrix->RowMap()));
  assert(rowMaps(i)->SameAs(sourceMatrix->RangeMap()));
  assert(colMaps(j)->SameAs(sourceMatrix->DomainMap()));
  
  //Estract data
  UInt    length = sourceMatrix->MaxNumEntries() + 1;
  UInt rowSource = sourceMatrix->NumMyRows();
  
  //Alloc
  double * values  = new double[length];
  int    * indices = new int[length]; 
  int numEntries, newRow;
  int ret;
  
  //Copy cycle
  for(UInt row = 0; row < rowSource; ++row)
  {
    //Estrazione
    ret = sourceMatrix->ExtractGlobalRowCopy(rowMaps(i)->GID(row), length, numEntries, values, indices);
    assert(ret == 0);
    
    //Aggiornamento indici
    newRow = rowMaps(i)->GID(row) + offsetRow(i);
    
    for(int k=0; k < numEntries; ++k)
    {
      assert(indices[k] < int(sourceMatrix->NumGlobalCols()));
      indices[k] += offsetCol(j);
    }
    
    //Memorizzazione 
    outMatrix->InsertGlobalValues(newRow, numEntries, values, indices);
  }
  
  delete[] values;
  delete[] indices;
  assembleOK = false;
}


void
epetraBlockMatrixAssemble::
setMatrixTransp(const Epetra_CrsMatrix & sourceMatrix, const UInt & i, const UInt & j)
{
  Epetra_BlockMap rowMap = *rowMaps(i);
  Epetra_Map * rowMapPtr = (Epetra_Map*) &rowMap;
  Epetra_CrsMatrix mat(sourceMatrix);
  Epetra_CrsMatrix * matT;
  Epetra_RowMatrixTransposer Traspositore(&mat);
  Traspositore.CreateTranspose(false,matT,rowMapPtr);
  
  setMatrix(*matT,i,j);  
  delete matT;
}

void
epetraBlockMatrixAssemble::
setMatrixTransp(const RCP_MATRIX & sourceMatrix, const UInt & i, const UInt & j)
{
  Epetra_BlockMap rowMap = *rowMaps(i);
  Epetra_Map * rowMapPtr = (Epetra_Map*) &rowMap;
  Epetra_CrsMatrix mat(*sourceMatrix);
  Epetra_CrsMatrix * matT;
  Epetra_RowMatrixTransposer Traspositore(&mat);
  Traspositore.CreateTranspose(false,matT,rowMapPtr);
  
  setMatrix(*matT,i,j);  
  delete matT;
}

void
epetraBlockMatrixAssemble::
assemble()
{
  outMatrix->FillComplete(*outColMap,*outRowMap);
  assembleOK = true;
}

epetraBlockMatrixAssemble::RCP_MAP
epetraBlockMatrixAssemble::
getRowMap() const
{
  assert(startupOk);
  return(outRowMap);
}

epetraBlockMatrixAssemble::RCP_MAP
epetraBlockMatrixAssemble::
getColMap() const
{
  assert(startupOk);
  return(outColMap);
}

epetraBlockMatrixAssemble::RCP_FEMATRIX
epetraBlockMatrixAssemble::
getMatrix() const
{
  assert(assembleOK);
  return(outMatrix);
}

void
epetraBlockMatrixAssemble::
printToFile(const string & name) const
{
  assert(assembleOK);
  EpetraExt::RowMatrixToMatlabFile(name.c_str(),*outMatrix);
}
