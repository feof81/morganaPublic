#include "epetraBlockVectorAssemble.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
epetraBlockVectorAssemble::
epetraBlockVectorAssemble()
{
  dimLoaded     = false;
  startupOk     = false;
  commDevLoaded = false;
}

epetraBlockVectorAssemble::
epetraBlockVectorAssemble(const UInt & NumRows)
{
  dimLoaded     = true;
  startupOk     = false;
  commDevLoaded = false;
  
  numRows = NumRows;
  
  rowLoaded.resize(NumRows);
  rowMaps.resize(NumRows);
  
  for(UInt i=1; i <= NumRows; ++i)
  { rowLoaded(i) = false; }
}

void
epetraBlockVectorAssemble::
setDim(const UInt & NumRows)
{
  dimLoaded = true;
  startupOk = false;
  
  numRows = NumRows;
  
  rowLoaded.resize(NumRows);  
  rowMaps.resize(NumRows);
  
  for(UInt i=1; i <= NumRows; ++i)
  { rowLoaded(i) = false; }
}

void
epetraBlockVectorAssemble::
setEpetraComm(const RCP_MPICOMM & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

void
epetraBlockVectorAssemble::
setEpetraComm(const Epetra_MpiComm & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcp(new Epetra_MpiComm(CommDev));
}



//_________________________________________________________________________________________________
// TOPOLOGY
//-------------------------------------------------------------------------------------------------
void
epetraBlockVectorAssemble::
setRowMap(const UInt & i, const Epetra_BlockMap & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numRows);
  assert(map.UniqueGIDs());
  
  rowLoaded(i) = true;
  rowMaps(i)   = Teuchos::rcp(new Epetra_BlockMap(map));
  
  startupOk = false;
}

void
epetraBlockVectorAssemble::
setRowMap(const UInt & i, const RCP_BLOCKMAP & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numRows);
  assert(map->UniqueGIDs());
  
  rowLoaded(i) = true;
  rowMaps(i)   = map;
  
  startupOk = false;
}

void
epetraBlockVectorAssemble::
startup()
{
  //Typedefs
  typedef Teuchos::RCP<int>  RCP_INT;
  
  //Checking
  assert(commDevLoaded);  
  assert(dimLoaded);
  
  for(UInt i=1; i <= numRows; ++i)
  { assert(rowLoaded(i)); }
  
  
  //Max gids
  offsetRow.resize(numRows);
  
  if(numRows != 0)
  { offsetRow(1) = 0; }
  
  for(UInt i=2; i <= numRows; ++i)
  { offsetRow(i) = offsetRow(i-1) + (rowMaps(i-1)->MaxAllGID() + 1); }
  
  
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
  
  
  //Alloc the final matrix
  outVector = Teuchos::rcp(new Epetra_FEVector(*outRowMap,1));
  
  //Flag
  startupOk  = true;
}



//_________________________________________________________________________________________________
// ASSEMBLE
//-------------------------------------------------------------------------------------------------
void
epetraBlockVectorAssemble::
setVector(const Epetra_MultiVector & sourceVector, const UInt & i)
{
  //Asserts
  assert(startupOk);
  assert(i >= 1); assert(i <= numRows);
  assert(rowMaps(i)->NumMyElements() == sourceVector.MyLength());
  
  //Estract data
  UInt rowSource = sourceVector.MyLength();
  
  //Alloc
  int newRow;
  double val;
  
  //Data extraction
  double Temp[rowSource];
  sourceVector.ExtractCopy(Temp,0);
  
  //Copy cycle
  for(UInt row = 0; row < rowSource; ++row)
  {
    val    = Temp[row];
    newRow = rowMaps(i)->GID(row) + offsetRow(i);
    
    outVector->SumIntoGlobalValue(newRow,0,val);
  }
}

void
epetraBlockVectorAssemble::
setVector(const RCP_VECTOR & sourceVector, const UInt & i)
{
  //Asserts
  assert(startupOk);
  assert(i >= 1); assert(i <= numRows);
  assert(rowMaps(i)->NumMyElements() == sourceVector->MyLength());
  
  //Estract data
  UInt rowSource = sourceVector->MyLength();
  
  //Alloc
  int newRow;
  double val;
  
  //Data extraction
  double Temp[rowSource];
  sourceVector->ExtractCopy(Temp,0);
  
  //Copy cycle
  for(UInt row = 0; row < rowSource; ++row)
  {
    val    = Temp[row];
    newRow = rowMaps(i)->GID(row) + offsetRow(i);
    
    outVector->SumIntoGlobalValue(newRow,0,val);
  }
}

epetraBlockVectorAssemble::RCP_MAP
epetraBlockVectorAssemble::
getRowMap() const
{
  assert(startupOk);
  return(outRowMap);
}

epetraBlockVectorAssemble::RCP_FEVECTOR
epetraBlockVectorAssemble::
getVector() const
{
  assert(startupOk);
  return(outVector);
}
