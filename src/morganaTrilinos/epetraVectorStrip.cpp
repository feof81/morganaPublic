#include "epetraVectorStrip.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
epetraVectorStrip::
epetraVectorStrip()
{
  dimLoaded     = false;
  startupOk     = false;
  commDevLoaded = false;
}

epetraVectorStrip::
epetraVectorStrip(const UInt & NumRows)
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
epetraVectorStrip::
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
epetraVectorStrip::
setEpetraComm(const RCP_MPICOMM & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

void
epetraVectorStrip::
setEpetraComm(const Epetra_MpiComm & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcp(new Epetra_MpiComm(CommDev));
}



//_________________________________________________________________________________________________
// TOPOLOGY
//-------------------------------------------------------------------------------------------------
void
epetraVectorStrip::
setRowMap(const UInt & i, const Epetra_BlockMap & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numRows);
  
  rowLoaded(i) = true;
  rowMaps(i)   = Teuchos::rcp(new Epetra_BlockMap(map));
  
  startupOk = false;
}

void
epetraVectorStrip::
setRowMap(const UInt & i, const RCP_BLOCKMAP & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numRows);
  
  rowLoaded(i) = true;
  rowMaps(i)   = map;
  
  startupOk = false;
}

void
epetraVectorStrip::
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
  
  //Flag
  startupOk  = true;
}



//_________________________________________________________________________________________________
// ASSEMBLE
//-------------------------------------------------------------------------------------------------
Epetra_MultiVector
epetraVectorStrip::
getVectorRef(const Epetra_MultiVector & inVector, const UInt & i)
{  
  assert(startupOk);
  
  //Extract data
  Epetra_BlockMap  inMap(inVector.Map());
  const Epetra_MultiVector & vectRef = inVector;
  double * vectPtr = vectRef[0];
  
  //Alloc
  Epetra_MultiVector outVector(*rowMaps(i),1);
  int gid, lid;
  double val;
  
  for(int k=0; k < rowMaps(i)->NumMyElements(); ++k)
  {
    gid = rowMaps(i)->GID(k) + offsetRow(i);
    lid = inMap.LID(gid);
    val = vectPtr[lid];
    
    outVector.SumIntoMyValue(k,0,val);
  }
  
  return(outVector);
}

epetraVectorStrip::RCP_VECTOR
epetraVectorStrip::
getVectorRcp(const RCP_VECTOR & inVector, const UInt & i)
{
  assert(startupOk);
  
  //Extract data
  Epetra_BlockMap  inMap(inVector->Map());
  const Epetra_MultiVector & vectRef = *inVector;
  double * vectPtr = vectRef[0];
  
  //Alloc
  RCP_VECTOR outVector = Teuchos::rcp(new Epetra_MultiVector(*rowMaps(i),1));
  
  int gid, lid;
  double val;
  
  for(int k=0; k < rowMaps(i)->NumMyElements(); ++k)
  {
    gid = rowMaps(i)->GID(k) + offsetRow(i);
    lid = inMap.LID(gid);
    val = vectPtr[lid];
    
    outVector->SumIntoMyValue(lid,0,val);
  }
  
  return(outVector);
}


