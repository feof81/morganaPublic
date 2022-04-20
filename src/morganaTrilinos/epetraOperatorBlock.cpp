/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "epetraOperatorBlock.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
epetraOperatorBlock::
epetraOperatorBlock()
{
  dimLoaded     = false;
  startupOk     = false;
  commDevLoaded = false;
}

epetraOperatorBlock::
epetraOperatorBlock(const UInt & NumRows, const UInt & NumCols)
{
  dimLoaded     = true;
  startupOk     = false;
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
epetraOperatorBlock::
setDim(const UInt & NumRows, const UInt & NumCols)
{
  dimLoaded  = true;
  startupOk  = false;
  
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
epetraOperatorBlock::
setEpetraComm(const RCP_MPICOMM & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

void
epetraOperatorBlock::
setEpetraComm(const Epetra_MpiComm & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcp(new Epetra_MpiComm(CommDev));
}


//_________________________________________________________________________________________________
// TOPOLOGY
//-------------------------------------------------------------------------------------------------
void
epetraOperatorBlock::
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
epetraOperatorBlock::
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
epetraOperatorBlock::
setColMap(const UInt & i, const Epetra_BlockMap & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numCols);
  assert(map.UniqueGIDs());
  
  colLoaded(i) = true;
  colMaps(i)   = Teuchos::rcp(new Epetra_BlockMap(map));
  
  startupOk = false;
}

void
epetraOperatorBlock::
setColMap(const UInt & i, const RCP_BLOCKMAP & map)
{
  assert(dimLoaded);
  assert(i >= 1);
  assert(i <= numCols);
  assert(map->UniqueGIDs());
  
  colLoaded(i) = true;
  colMaps(i)   = map;
  
  startupOk = false;
}

void
epetraOperatorBlock::
startup()
{
  //Typedefs---------------------------------------------------------------------------------------
  typedef Teuchos::RCP<int>  RCP_INT;
  
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
  { offsetRow(i) = offsetRow(i-1) + (rowMaps(i-1)->MaxAllGID() + 1); }
  
  if(numCols != 0)
  { offsetCol(1) = 0; }
  
  for(UInt i=2; i <= numCols; ++i)
  { offsetCol(i) = offsetCol(i-1) + (colMaps(i-1)->MaxAllGID() + 1); }
  
  //Build row--------------------------------------------------------------------------------------
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
  
  //Build col--------------------------------------------------------------------------------------
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
  assert(outColMap->UniqueGIDs());
  
  //Vector strip - Cols----------------------------------------------------------------------------
  vectStrip.setDim(numCols);
  vectStrip.setEpetraComm(commDev);
  
  for(UInt i=1; i <= colMaps.size(); ++i)
  { vectStrip.setRowMap(i,colMaps(i)); }
  
  vectStrip.startup();
 
  //Vector assemble startup - Rows-----------------------------------------------------------------
  blockVect.setDim(numRows);
  blockVect.setEpetraComm(commDev);
  
  for(UInt i=1; i <= rowMaps.size(); ++i)
  { blockVect.setRowMap(i,rowMaps(i)); }
  
  blockVect.startup();
  
  //Alloc matrix-----------------------------------------------------------------------------------
  operatorMatrix.resize(numRows);
  boolMatrix.resize(numRows);
  
  for(UInt i=1; i <= numRows; ++i)
  {
    operatorMatrix(i).resize(numCols);
    boolMatrix(i).resize(numCols);
    
    for(UInt j=1; j <= numCols; ++j)
    { boolMatrix(i)(j) = false; }
  }
  
  //Create rows subvectors-------------------------------------------------------------------------
  subRowVect.resize(numRows);
  
  for(UInt i=1; i <= numRows; ++i)
  { subRowVect(i) = Teuchos::rcp(new Epetra_MultiVector(*rowMaps(i),1)); }
  
  //Create cols subvectors-------------------------------------------------------------------------
  subColVect.resize(numCols);
  
  for(UInt i=1; i <= numCols; ++i)
  { subColVect(i) = Teuchos::rcp(new Epetra_MultiVector(*colMaps(i),1)); }

  //Flag-------------------------------------------------------------------------------------------
  startupOk = true;
}

void
epetraOperatorBlock::
setOperator(const RCP_OPERATOR & sourceOperator, const UInt & i, const UInt & j)
{
  //Asserts----------------------------------------------------------------------------------------
  assert(startupOk);
  assert(i >= 1); assert(i <= numRows);
  assert(i >= 1); assert(j <= numCols);
  assert(rowMaps(i)->SameAs(sourceOperator->OperatorRangeMap()));
  assert(colMaps(j)->SameAs(sourceOperator->OperatorDomainMap()));
  
  //Alloc------------------------------------------------------------------------------------------
  operatorMatrix(i)(j) = sourceOperator;
  boolMatrix(i)(j)     = true;
}


//_________________________________________________________________________________________________
// INHERITED FUNCTIONS
//-------------------------------------------------------------------------------------------------
int
epetraOperatorBlock::
SetUseTranspose(bool UseTranspose)
{
  assert(startupOk);
  
  return(false);
}

int
epetraOperatorBlock::
Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  //Assert-----------------------------------------------------------------------------------------
  assert(startupOk);
  
  /*for(UInt i=1; i <= numRows; ++i)
  {
    for(UInt j=1; j <= numCols; ++j)
    { assert(boolMatrix(i)(j)); }
  }*/
  
  //Vector stripping-------------------------------------------------------------------------------
  for(UInt i=1; i <= numCols; ++i)
  { *subColVect(i) = vectStrip.getVectorRef(X,i); }
  
  //Operator computing-----------------------------------------------------------------------------
  for(UInt i=1; i <= numRows; ++i)
  {
    subRowVect(i)->Scale(0.0);
    
    Epetra_MultiVector tempRow(*subRowVect(i));
    
    for(UInt j=1; j <= numCols; ++j)
    {
      if(boolMatrix(i)(j))
      {
        operatorMatrix(i)(j)->Apply(*subColVect(j),tempRow);
        subRowVect(i)->Update(1.0,tempRow,1.0);
      }
    }
  }
  
  //Build final vector-----------------------------------------------------------------------------
  blockVect.startup();
  
  for(UInt i=1; i <= numRows; ++i)
  { blockVect.setVector(subRowVect(i),i); }
  
  Y = *blockVect.getVector();
  
  return(0);
}

int
epetraOperatorBlock::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  assert(startupOk);
  assert(false);
  
  return(1);
}

double
epetraOperatorBlock::
NormInf() const
{
  //Assert-----------------------------------------------------------------------------------------
  assert(startupOk);
  
  /*for(UInt i=1; i <= numRows; ++i)
  {
    for(UInt j=1; j <= numCols; ++j)
    { assert(boolMatrix(i)(j)); }
  }*/
  
  //Compute----------------------------------------------------------------------------------------
  Real maxNorm = 1.0;
  
  for(UInt i=1; i <= numRows; ++i)
  {
    for(UInt j=1; j <= numCols; ++j)
    {
      if(boolMatrix(i)(j))
      { maxNorm = std::max(maxNorm, operatorMatrix(i)(j)->NormInf()); }
    }
  }
  
  return(maxNorm);
}


const char *
epetraOperatorBlock::
Label() const
{
  assert(startupOk);
  
  return("operatorBlok");
}

bool
epetraOperatorBlock::
UseTranspose() const
{
  assert(startupOk);
  
  return(false);
}

bool
epetraOperatorBlock::
HasNormInf() const
{
  assert(startupOk);
  
  return(true);
}

const Epetra_Comm &
epetraOperatorBlock::
Comm() const
{
  assert(startupOk);
  assert(boolMatrix(1)(1));
  
  return(operatorMatrix(1)(1)->Comm());
}

const Epetra_Map &
epetraOperatorBlock::
OperatorDomainMap() const
{
  assert(startupOk);
  
  return(*outColMap);
}

const Epetra_Map &
epetraOperatorBlock::
OperatorRangeMap() const
{
  assert(startupOk);
  
  return(*outRowMap);
}
