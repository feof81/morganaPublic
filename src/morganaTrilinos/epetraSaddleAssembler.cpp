#include "epetraSaddleAssembler.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
epetraSaddleAssembler::
epetraSaddleAssembler() : epetraBlockMatrixAssemble()
{
}

epetraSaddleAssembler::
epetraSaddleAssembler(const UInt & Dim) : epetraBlockMatrixAssemble(Dim,Dim)
{
}

void
epetraSaddleAssembler::
setDim(const UInt & Dim)
{
  epetraBlockMatrixAssemble::setDim(Dim,Dim);
}

void
epetraSaddleAssembler::
setEpetraComm(const RCP_MPICOMM & CommDev)
{
  epetraBlockMatrixAssemble::setEpetraComm(CommDev);
}

void
epetraSaddleAssembler::
setEpetraComm(const Epetra_MpiComm & CommDev)
{
  epetraBlockMatrixAssemble::setEpetraComm(CommDev);
}



//_________________________________________________________________________________________________
// TOPOLOGY
//-------------------------------------------------------------------------------------------------
void
epetraSaddleAssembler::
setMap(const UInt & i, const Epetra_BlockMap & map)
{
  epetraBlockMatrixAssemble::setRowMap(i,map);
  epetraBlockMatrixAssemble::setColMap(i,map);
}

void
epetraSaddleAssembler::
setMap(const UInt & i, const RCP_BLOCKMAP & map)
{
  epetraBlockMatrixAssemble::setRowMap(i,map);
  epetraBlockMatrixAssemble::setColMap(i,map);
}

void
epetraSaddleAssembler::
startup()
{
  epetraBlockMatrixAssemble::startup();
}



//_________________________________________________________________________________________________
// ASSEMBLE
//-------------------------------------------------------------------------------------------------
void
epetraSaddleAssembler::
setDiagMatrix(const Epetra_CrsMatrix & sourceMatrix, const UInt & i)
{
  assert(sourceMatrix.NumGlobalRows() == sourceMatrix.NumGlobalCols());
  
  epetraBlockMatrixAssemble::setMatrix(sourceMatrix, i, i);
}

void
epetraSaddleAssembler::
setDiagMatrix(const RCP_MATRIX & sourceMatrix, const UInt & i)
{
  assert(sourceMatrix->NumGlobalRows() == sourceMatrix->NumGlobalCols());
  
  epetraBlockMatrixAssemble::setMatrix(sourceMatrix, i, i);
}

void
epetraSaddleAssembler::
setMTM(const Epetra_CrsMatrix & M, const UInt & i, const Real & scale)
{
  //Assert
  assert(M.NumGlobalCols() == (rowMaps(i)->MaxAllGID() + 1));
  
  //Transpose matrix
  Epetra_CrsMatrix * MT;
  Epetra_CrsMatrix sourceTemp(M);
  Epetra_Map * newMap = static_cast<Epetra_Map*>(rowMaps(i).get());  
    
  Epetra_RowMatrixTransposer Transpositor(&sourceTemp);
  Transpositor.CreateTranspose(true,MT,newMap);
  
  //Product
  Epetra_FECrsMatrix MTM(Copy,*newMap,1);
  
  EpetraExt::MatrixMatrix multiplier;
  multiplier.Multiply(*MT,false, M,false, MTM);
  
  //Assemble
  MTM.Scale(scale);
  epetraBlockMatrixAssemble::setMatrix(MTM, i, i);
}

void
epetraSaddleAssembler::
setMTM(const RCP_MATRIX & M, const UInt & i, const Real & scale)
{
  //Assert
  assert(M->NumGlobalCols() == (rowMaps(i)->MaxAllGID() + 1));
  
  //Transpose matrix
  Epetra_CrsMatrix * MT;
  Epetra_CrsMatrix sourceTemp(*M);
  Epetra_Map * newMap = static_cast<Epetra_Map*>(rowMaps(i).get());  
    
  Epetra_RowMatrixTransposer Transpositor(&sourceTemp);
  Transpositor.CreateTranspose(true,MT,newMap);
  
  //Product
  Epetra_FECrsMatrix MTM(Copy,*newMap,1);
  
  EpetraExt::MatrixMatrix multiplier;
  multiplier.Multiply(*MT,false, *M,false, MTM);
  
  //Assemble
  MTM.Scale(scale);
  epetraBlockMatrixAssemble::setMatrix(MTM, i, i);
}

void
epetraSaddleAssembler::
setMMT(const Epetra_CrsMatrix & M, const UInt & i, const Real & scale)
{
  //Assert
  assert(M.NumGlobalRows() == (rowMaps(i)->MaxAllGID() + 1));
  
  //Product
  Epetra_FECrsMatrix MTM(Copy,M.RowMap(),1);
  
  //Multiply
  EpetraExt::MatrixMatrix multiplier;
  multiplier.Multiply(M,false, M,true, MTM);
  
  //Assemble
  MTM.Scale(scale);
  epetraBlockMatrixAssemble::setMatrix(MTM, i, i);
}

void
epetraSaddleAssembler::
setMMT(const RCP_MATRIX & M, const UInt & i, const Real & scale)
{
  //Assert
  assert(M->NumGlobalRows() == (rowMaps(i)->MaxAllGID() + 1));
  
  //Product
  Epetra_FECrsMatrix MTM(Copy,M->RowMap(),1);
  
  //Multiply
  EpetraExt::MatrixMatrix multiplier;
  multiplier.Multiply(*M,false, *M,true, MTM);
  
  //Assemble
  MTM.Scale(scale);
  epetraBlockMatrixAssemble::setMatrix(MTM, i, i);
}

void
epetraSaddleAssembler::
setEye(const UInt & i, const Real & scale)
{
  //Build matrix
  Epetra_Map * newMap = static_cast<Epetra_Map*>(rowMaps(i).get());
  Epetra_FECrsMatrix Eye(Copy,*newMap,1);
  
  //Alloc
  double values[1];
  int numEntries, newRow, indices[1];
  UInt numRows = Eye.NumMyRows();
  
  numEntries = 1;
  values[0]  = scale;
  
  //Fill loop
  for(UInt row = 0; row < numRows; ++row)
  {
    newRow     = rowMaps(i)->GID(row);
    indices[0] = rowMaps(i)->GID(row);
    
    Eye.InsertGlobalValues(newRow, numEntries, values, indices);
  }
  
  Eye.FillComplete(*newMap,*newMap);
  
  //Assemble
  epetraBlockMatrixAssemble::setMatrix(Eye, i, i);
}

void
epetraSaddleAssembler::
setSymmMatrix(const Epetra_CrsMatrix & sourceMatrix, const UInt & i, const UInt & j)
{
  //Transpose matrix
  Epetra_CrsMatrix * transposeMatrix;
  Epetra_CrsMatrix sourceTemp(sourceMatrix);
  Epetra_Map * newMap = static_cast<Epetra_Map*>(rowMaps(j).get());  
    
  Epetra_RowMatrixTransposer Transpositor(&sourceTemp);
  Transpositor.CreateTranspose(true,transposeMatrix,newMap);
  
  //Assemble
  epetraBlockMatrixAssemble::setMatrix(sourceMatrix,     i, j);  
  epetraBlockMatrixAssemble::setMatrix(*transposeMatrix, j, i);
  
  //Deleting
  delete transposeMatrix;
}

void
epetraSaddleAssembler::
setSymmMatrix(const RCP_MATRIX & sourceMatrix, const UInt & i, const UInt & j)
{
  //Transpose matrix
  Epetra_CrsMatrix * transposeMatrix;
  Epetra_CrsMatrix sourceTemp(*sourceMatrix);
  Epetra_Map * newMap = static_cast<Epetra_Map*>(rowMaps(j).get());
    
  Epetra_RowMatrixTransposer Transpositor(&sourceTemp);
  Transpositor.CreateTranspose(true,transposeMatrix,newMap);
  
  //Assemble
  epetraBlockMatrixAssemble::setMatrix(*sourceMatrix,    i, j);
  epetraBlockMatrixAssemble::setMatrix(*transposeMatrix, j, i);
  
  //Deleting
  delete transposeMatrix;
}

void
epetraSaddleAssembler::
assemble()
{
  epetraBlockMatrixAssemble::assemble();
}



//_________________________________________________________________________________________________
// RETURN
//-------------------------------------------------------------------------------------------------
epetraSaddleAssembler::RCP_MAP
epetraSaddleAssembler::
getMap() const
{
  return(epetraBlockMatrixAssemble::getRowMap());
}

epetraSaddleAssembler::RCP_FEMATRIX
epetraSaddleAssembler::
getMatrix() const
{
  return(epetraBlockMatrixAssemble::getMatrix());
}

