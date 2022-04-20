#include "tpetraVectorStrip.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
tpetraVectorStrip::
tpetraVectorStrip()
{
  dimLoaded     = false;
  startupOk     = false;
  commDevLoaded = false;
}

tpetraVectorStrip::
tpetraVectorStrip(const UInt & NumRows)
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
tpetraVectorStrip::
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
tpetraVectorStrip::
setEpetraComm(const RCP_MPICOMM & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

void
tpetraVectorStrip::
setEpetraComm(const COMM & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcp(new COMM(CommDev));
}


//_________________________________________________________________________________________________
// TOPOLOGY
//-------------------------------------------------------------------------------------------------
void
tpetraVectorStrip::
setRowMap(const UInt & i, const MAP & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numRows);
  
  rowLoaded(i) = true;
  rowMaps(i)   = Teuchos::rcp(new MAP(map));
  
  startupOk = false;
}

void
tpetraVectorStrip::
setRowMap(const UInt & i, const RCP_MAP & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numRows);
  
  rowLoaded(i) = true;
  rowMaps(i)   = map;
  
  startupOk = false;
}

void
tpetraVectorStrip::
startup()
{  
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
  { offsetRow(i) = offsetRow(i-1) + (rowMaps(i-1)->getMaxAllGlobalIndex() + 1); }
  
  //Flag
  startupOk  = true;
}



//_________________________________________________________________________________________________
// ASSEMBLE
//-------------------------------------------------------------------------------------------------
tpetraVectorStrip::VECTOR
tpetraVectorStrip::
getVectorRef(const VECTOR & inVector, const UInt & i)
{  
  assert(startupOk);
  
  //Extract data
  RCP_MAP inMap = inVector.getMap();
  Teuchos::ArrayRCP<const double> val = inVector.getData(0);
  
  //Alloc
  VECTOR outVector(rowMaps(i),1);
  int gid, lid;
  
  for(int k=0; k < rowMaps(i)->getNodeNumElements(); ++k)
  {
    gid = rowMaps(i)->getGlobalElement(k) + offsetRow(i);
    lid = inMap->getLocalElement(gid);

    outVector.sumIntoLocalValue(k,0,val[lid]);
  }
  
  return(outVector);
}

tpetraVectorStrip::RCP_VECTOR
tpetraVectorStrip::
getVectorRcp(const RCP_VECTOR & inVector, const UInt & i)
{
  RCP_VECTOR outVector(new VECTOR(getVectorRef(*inVector,i)));
  return(outVector);
}

