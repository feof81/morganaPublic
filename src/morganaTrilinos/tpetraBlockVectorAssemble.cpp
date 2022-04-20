#include "tpetraBlockVectorAssemble.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
tpetraBlockVectorAssemble::
tpetraBlockVectorAssemble()
{
  dimLoaded     = false;
  startupOk     = false;
  commDevLoaded = false;
}

tpetraBlockVectorAssemble::
tpetraBlockVectorAssemble(const UInt & NumRows)
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
tpetraBlockVectorAssemble::
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
tpetraBlockVectorAssemble::
setTpetraComm(const RCP_MPICOMM & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

void
tpetraBlockVectorAssemble::
setTpetraComm(const COMM & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcp(new COMM(CommDev));
}


//_________________________________________________________________________________________________
// TOPOLOGY
//-------------------------------------------------------------------------------------------------
void
tpetraBlockVectorAssemble::
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
tpetraBlockVectorAssemble::
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
tpetraBlockVectorAssemble::
startup()
{
  //Checking---------------------------------------------------------
  assert(commDevLoaded);  
  assert(dimLoaded);
  
  for(UInt i=1; i <= numRows; ++i)
  { assert(rowLoaded(i)); }
  
  //Max gids---------------------------------------------------------
  offsetRow.resize(numRows);
  
  if(numRows != 0)
  { offsetRow(1) = 0; }
  
  for(UInt i=2; i <= numRows; ++i)
  { offsetRow(i) = offsetRow(i-1) + (rowMaps(i-1)->getMaxAllGlobalIndex()  + 1); }
  
  //Build row--------------------------------------------------------
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
  
  //Alloc the final matrix-------------------------------------------
  outVector = Teuchos::rcp(new VECTOR(outRowMap,1));
  
  //Flag-------------------------------------------------------------
  startupOk  = true;
}


//_________________________________________________________________________________________________
// ASSEMBLE
//-------------------------------------------------------------------------------------------------
void
tpetraBlockVectorAssemble::
setVector(const VECTOR & sourceVector, const UInt & i)
{
  //Asserts
  assert(startupOk);
  assert(i >= 1); assert(i <= numRows);
  assert(rowMaps(i)->getNodeNumElements() == sourceVector.getLocalLength());
  
  //Estract data
  UInt rowSource = sourceVector.getLocalLength();
  
  //Alloc
  int newRow;
  Teuchos::ArrayRCP<const TPETRA_SCALAR> val;
  
  //Copy cycle
  for(UInt row = 0; row < rowSource; ++row)
  {
    val    = sourceVector.getData(0);
    newRow = rowMaps(i)->getGlobalElement(row) + offsetRow(i);
    
    outVector->sumIntoGlobalValue(newRow,0,val[row]);
  }
}

void
tpetraBlockVectorAssemble::
setVector(const RCP_VECTOR & sourceVector, const UInt & i)
{
  setVector(*sourceVector,i);
}

tpetraBlockVectorAssemble::RCP_MAP
tpetraBlockVectorAssemble::
getRowMap() const
{
  assert(startupOk);
  return(outRowMap);
}

tpetraBlockVectorAssemble::RCP_VECTOR
tpetraBlockVectorAssemble::
getVector() const
{
  assert(startupOk);
  return(outVector);
}
