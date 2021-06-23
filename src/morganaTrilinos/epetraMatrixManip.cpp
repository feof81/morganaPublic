/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <set>
#include <complex>

#include "epetraMatrixManip.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseMatrix.h"


epetraMatrixManip::
epetraMatrixManip()
{
}

Epetra_CrsMatrix
epetraMatrixManip::
transpose(const Epetra_CrsMatrix & A) const
{
  Epetra_Map rowMapA = A.DomainMap();
  
  Epetra_CrsMatrix Acopy(A);
  Epetra_CrsMatrix * AT;
  Epetra_RowMatrixTransposer Traspositore(&Acopy);
  Traspositore.CreateTranspose(false, AT, &rowMapA);
  
  Epetra_CrsMatrix ATcopy(*AT);
  delete(AT);
  
  return(ATcopy);
}
    
epetraMatrixManip::RCP
epetraMatrixManip::
transpose(const RCP_CONST & A) const
{
  Epetra_CrsMatrix M = transpose(*A);
  return(Teuchos::rcp(new Epetra_CrsMatrix(M)));
}

Epetra_CrsMatrix
epetraMatrixManip::
diag(const Epetra_CrsMatrix & A) const
{
  assert(A.NumGlobalRows() == A.NumGlobalCols());
  
  Epetra_CrsMatrix out(Copy,A.RowMap(),1);
  Epetra_Vector    diag(A.RowMap());
  Epetra_Map       rowMap(A.RowMap());
  
  A.ExtractDiagonalCopy(diag);
  Epetra_Vector& Dref = diag;
  
  int indexI, indexJ[1];
  double val[1];
  
  for(int i=0; i < A.NumMyRows(); ++i)
  {
    val[0]    = Dref[i];
    indexI    = rowMap.GID(i);
    indexJ[0] = rowMap.GID(i);
    
    out.InsertGlobalValues(indexI,1,val,indexJ);
  }
  
  out.FillComplete(A.DomainMap(),A.RangeMap());
  
  return(out);
}

epetraMatrixManip::RCP
epetraMatrixManip::
diag(const RCP_CONST & A) const
{
  Epetra_CrsMatrix M = diag(*A);
  return(Teuchos::rcp(new Epetra_CrsMatrix(M)));
}

Epetra_CrsMatrix 
epetraMatrixManip::
getEyeRef(const Real       & scale,
          const Epetra_Map & Map) const
{
  Epetra_FECrsMatrix Eye(Copy,Map,1);
  
  //Alloc
  double values[1];
  int numEntries, newRow, indices[1];
  UInt numRows = Eye.NumMyRows();
  
  numEntries = 1;
  values[0]  = scale;
  
  //Fill loop
  for(UInt row = 0; row < numRows; ++row)
  {
    newRow     = Map.GID(row);
    indices[0] = Map.GID(row);
    
    Eye.InsertGlobalValues(newRow, numEntries, values, indices);
  }
  
  Eye.FillComplete(Map,Map);
  
  return(Eye);
}
    
epetraMatrixManip::RCP
epetraMatrixManip::
getEyeRcp(const Real       & scale,
          const Epetra_Map & Map) const
{
  RCP Eye = Teuchos::rcp(new Epetra_FECrsMatrix(Copy,Map,1));
  
  //Alloc
  double values[1];
  int numEntries, newRow, indices[1];
  UInt numRows = Eye->NumMyRows();
  
  numEntries = 1;
  values[0]  = scale;
  
  //Fill loop
  for(UInt row = 0; row < numRows; ++row)
  {
    newRow     = Map.GID(row);
    indices[0] = Map.GID(row);
    
    Eye->InsertGlobalValues(newRow, numEntries, values, indices);
  }
  
  Eye->FillComplete(Map,Map);
  
  return(Eye);
}

Epetra_CrsMatrix
epetraMatrixManip::
diagInv(const Epetra_CrsMatrix & A) const
{
  assert(A.NumGlobalRows() == A.NumGlobalCols());
  
  Epetra_CrsMatrix out(Copy,A.RowMap(),1);
  Epetra_Vector    diag(A.RowMap());
  Epetra_Map       rowMap(A.RowMap());
  
  A.ExtractDiagonalCopy(diag);
  Epetra_Vector& Dref = diag;
  
  int indexI, indexJ[1];
  double val[1];
  
  for(int i=0; i < A.NumMyRows(); ++i)
  {
    val[0]    = 1.0 / Dref[i];
    indexI    = rowMap.GID(i);
    indexJ[0] = rowMap.GID(i);
    
    if(abs(Dref[i]) == 0) {cout << "Warning! divide by zero, diagInv" << endl;}
    
    out.InsertGlobalValues(indexI,1,val,indexJ);
  }
  
  out.FillComplete(A.DomainMap(),A.RangeMap());
  
  return(out);
}

epetraMatrixManip::RCP
epetraMatrixManip::
diagInv(const RCP_CONST & A) const
{
  Epetra_CrsMatrix M = diagInv(*A);
  return(Teuchos::rcp(new Epetra_CrsMatrix(M)));
}

Epetra_CrsMatrix
epetraMatrixManip::
multiply(const Epetra_CrsMatrix & A,
         const Epetra_CrsMatrix & B) const
{
  Epetra_Map rowMapA = A.RowMap();
  Epetra_CrsMatrix AB(Copy,rowMapA,1);  
  EpetraExt::MatrixMatrix multiplier;
  
  multiplier.Multiply(A,false, B,false, AB);
  
  return(AB);
}

epetraMatrixManip::RCP
epetraMatrixManip::
multiply(const RCP_CONST & A,
         const RCP_CONST & B) const
{
  Epetra_CrsMatrix M = multiply(*A,*B);
  return(Teuchos::rcp(new Epetra_CrsMatrix(M)));
}

Epetra_CrsMatrix
epetraMatrixManip::
multiply(const Epetra_CrsMatrix & A,
         const Epetra_CrsMatrix & B,
         const Epetra_CrsMatrix & C) const
{
  Epetra_Map rowMapA = A.RowMap();
  Epetra_CrsMatrix AB(Copy,rowMapA,1);  
  Epetra_CrsMatrix ABC(Copy,rowMapA,1);  
  EpetraExt::MatrixMatrix multiplier;
  
  multiplier.Multiply(A,false,  B,false, AB);
  multiplier.Multiply(AB,false, C,false, ABC);
  
  return(ABC);
}

epetraMatrixManip::RCP
epetraMatrixManip::
multiply(const RCP_CONST & A,
         const RCP_CONST & B,
         const RCP_CONST & C) const
{
  Epetra_CrsMatrix M = multiply(*A,*B,*C);
  return(Teuchos::rcp(new Epetra_CrsMatrix(M)));
}

Epetra_CrsMatrix
epetraMatrixManip::
shurDiag(const Epetra_CrsMatrix & A,
         const Epetra_CrsMatrix & B) const
{
  //Maps
  Epetra_Map rowMapA = A.RowMap();
  Epetra_Map rowMapB = B.RowMap();
  
  //Transpose 
  Epetra_CrsMatrix Bcopy(B);
  Epetra_CrsMatrix * BT;
  Epetra_RowMatrixTransposer Traspositore(&Bcopy);
  Traspositore.CreateTranspose(false, BT, &rowMapA);
  
  //Diag
  Epetra_CrsMatrix D = diagInv(A);
  
  //Assemble
  Epetra_CrsMatrix O = multiply(B,D,*BT);
  delete BT;
  
  return(O);
}

epetraMatrixManip::RCP
epetraMatrixManip::
shurDiag(const RCP_CONST & A,
         const RCP_CONST & B) const
{
  Epetra_CrsMatrix M = shurDiag(*A,*B);
  return(Teuchos::rcp(new Epetra_CrsMatrix(M)));
}

Epetra_CrsMatrix
epetraMatrixManip::
shurDiag(const Epetra_CrsMatrix & A,
         const Epetra_CrsMatrix & B,
         const Epetra_CrsMatrix & C) const
{
  //Maps
  Epetra_Map rowMapA = A.RowMap();
  Epetra_Map rowMapB = B.RowMap();
  
  //Transpose 
  Epetra_CrsMatrix Ccopy(C);
  Epetra_CrsMatrix * CT;
  Epetra_RowMatrixTransposer Traspositore(&Ccopy);
  Traspositore.CreateTranspose(false, CT, &rowMapA);
  
  //Diag
  Epetra_CrsMatrix D = diagInv(A);
  
  //Assemble
  Epetra_CrsMatrix O = multiply(B,D,*CT);
  delete CT;
  
  return(O);
}

epetraMatrixManip::RCP
epetraMatrixManip::
shurDiag(const RCP_CONST & A,
         const RCP_CONST & B,
         const RCP_CONST & C) const
{
  Epetra_CrsMatrix M = shurDiag(*A,*B,*C);
  return(Teuchos::rcp(new Epetra_CrsMatrix(M)));
}

Epetra_CrsMatrix
epetraMatrixManip::
shur(const Epetra_CrsMatrix & A,
     const Epetra_CrsMatrix & invB,
     const Epetra_CrsMatrix & C) const
{
  //Maps
  Epetra_Map rowMap = invB.RowMap();
  
  //Transpose 
  Epetra_CrsMatrix Ccopy(C);
  Epetra_CrsMatrix * CT;
  Epetra_RowMatrixTransposer Traspositore(&Ccopy);
  Traspositore.CreateTranspose(false, CT, &rowMap);
  
  //Assemble
  Epetra_CrsMatrix O = multiply(A,invB,*CT);
  delete CT;
  
  return(O);
}
    
epetraMatrixManip::RCP
epetraMatrixManip::
shur(const RCP_CONST & A,
     const RCP_CONST & invB,
     const RCP_CONST & C) const
{
  Epetra_CrsMatrix M = shur(*A,*invB,*C);
  return(Teuchos::rcp(new Epetra_CrsMatrix(M)));
}

Epetra_CrsMatrix
epetraMatrixManip::
BCT(const Epetra_CrsMatrix & B, 
    const Epetra_CrsMatrix & C) const
{
  //Maps
  Epetra_Map rowMapB = B.RowMap();
  Epetra_Map colMapB = B.DomainMap();
  
  //Transpose 
  Epetra_CrsMatrix Ccopy(C);
  Epetra_CrsMatrix * CT;
  Epetra_RowMatrixTransposer Traspositore(&Ccopy);
  Traspositore.CreateTranspose(false, CT, &colMapB);
  
  //Multiply
  Epetra_CrsMatrix BCt(Copy,rowMapB,1);
  
  EpetraExt::MatrixMatrix multiplier;
  multiplier.Multiply(B,false, *CT,false, BCt);
  
  return(BCt);
}
    
epetraMatrixManip::RCP
epetraMatrixManip::
BCT(const RCP_CONST & B,
    const RCP_CONST & C) const
{
  Epetra_CrsMatrix M = BCT(*B,*C);
  return(Teuchos::rcp(new Epetra_CrsMatrix(M)));
}

Epetra_CrsMatrix
epetraMatrixManip::
inv4x4(const Epetra_CrsMatrix & A) const
{
  //Checks
  assert(A.NumGlobalRows() == A.NumGlobalCols());
  assert((A.NumGlobalRows() % 4) == 0);
  
  //Extract data
  int length    = A.GlobalMaxNumEntries();
  int rowSource = A.NumGlobalRows();
  int numBlocks = rowSource / 4;
  
  //Alloc inv matrix
  Epetra_CrsMatrix invA(Copy,A.RowMap(),4);
  
  //Alloc
  Epetra_SerialDenseMatrix M(4,4);
  Epetra_SerialDenseSolver solver;
  
  int row0, row1, row2, row3;
  int numEntries, indices[length];
  double values[length];
  
  
  //Main loop
  for(int i=0; i < numBlocks; ++i)
  {
    if(A.RowMap().MyGID(i))
    {
      row0 = i;
      row1 = numBlocks + i;
      row2 = (2 * numBlocks) + i;
      row3 = (3 * numBlocks) + i;
      
      //Extract local row 1----------------------------------
      A.ExtractGlobalRowCopy(row0,length,numEntries,values,indices);
    
      assert(numEntries == 4);
    
      M(0,0) = values[0] * (indices[0] == row0)  +  values[1] * (indices[1] == row0)  +  values[2] * (indices[2] == row0)  +  values[3] * (indices[3] == row0);
      M(0,1) = values[0] * (indices[0] == row1)  +  values[1] * (indices[1] == row1)  +  values[2] * (indices[2] == row1)  +  values[3] * (indices[3] == row1);
      M(0,2) = values[0] * (indices[0] == row2)  +  values[1] * (indices[1] == row2)  +  values[2] * (indices[2] == row2)  +  values[3] * (indices[3] == row2);
      M(0,3) = values[0] * (indices[0] == row3)  +  values[1] * (indices[1] == row3)  +  values[2] * (indices[2] == row3)  +  values[3] * (indices[3] == row3);
      
      //Extract local row 2----------------------------------
      A.ExtractGlobalRowCopy(row1,length,numEntries,values,indices);
    
      assert(numEntries == 4);
    
      M(1,0) = values[0] * (indices[0] == row0)  +  values[1] * (indices[1] == row0)  +  values[2] * (indices[2] == row0)  +  values[3] * (indices[3] == row0);
      M(1,1) = values[0] * (indices[0] == row1)  +  values[1] * (indices[1] == row1)  +  values[2] * (indices[2] == row1)  +  values[3] * (indices[3] == row1);
      M(1,2) = values[0] * (indices[0] == row2)  +  values[1] * (indices[1] == row2)  +  values[2] * (indices[2] == row2)  +  values[3] * (indices[3] == row2);
      M(1,3) = values[0] * (indices[0] == row3)  +  values[1] * (indices[1] == row3)  +  values[2] * (indices[2] == row3)  +  values[3] * (indices[3] == row3);
    
      //Extract local row 3----------------------------------
      A.ExtractGlobalRowCopy(row2,length,numEntries,values,indices);
    
      assert(numEntries == 4);
    
      M(2,0) = values[0] * (indices[0] == row0)  +  values[1] * (indices[1] == row0)  +  values[2] * (indices[2] == row0)  +  values[3] * (indices[3] == row0);
      M(2,1) = values[0] * (indices[0] == row1)  +  values[1] * (indices[1] == row1)  +  values[2] * (indices[2] == row1)  +  values[3] * (indices[3] == row1);
      M(2,2) = values[0] * (indices[0] == row2)  +  values[1] * (indices[1] == row2)  +  values[2] * (indices[2] == row2)  +  values[3] * (indices[3] == row2);
      M(2,3) = values[0] * (indices[0] == row3)  +  values[1] * (indices[1] == row3)  +  values[2] * (indices[2] == row3)  +  values[3] * (indices[3] == row3);
    
      //Extract local row 4----------------------------------
      A.ExtractGlobalRowCopy(row3,length,numEntries,values,indices);
    
      assert(numEntries == 4);
    
      M(3,0) = values[0] * (indices[0] == row0)  +  values[1] * (indices[1] == row0)  +  values[2] * (indices[2] == row0)  +  values[3] * (indices[3] == row0);
      M(3,1) = values[0] * (indices[0] == row1)  +  values[1] * (indices[1] == row1)  +  values[2] * (indices[2] == row1)  +  values[3] * (indices[3] == row1);
      M(3,2) = values[0] * (indices[0] == row2)  +  values[1] * (indices[1] == row2)  +  values[2] * (indices[2] == row2)  +  values[3] * (indices[3] == row2);
      M(3,3) = values[0] * (indices[0] == row3)  +  values[1] * (indices[1] == row3)  +  values[2] * (indices[2] == row3)  +  values[3] * (indices[3] == row3);
      
      //Matrix inversion-------------------------------------
      solver.SetMatrix(M);
      solver.Invert();
      
      //Columns set------------------------------------------
      numEntries = 4;
      indices[0] = row0;
      indices[1] = row1;
      indices[2] = row2;
      indices[3] = row3;
    
      //Insert row 1-----------------------------------------
      values[0] = M(0,0);
      values[1] = M(0,1);
      values[2] = M(0,2);
      values[3] = M(0,3);
    
      invA.InsertGlobalValues(row0,numEntries,values,indices);
    
      //Insert row 2-----------------------------------------
      values[0] = M(1,0);
      values[1] = M(1,1);
      values[2] = M(1,2);
      values[3] = M(1,3);
    
      invA.InsertGlobalValues(row1,numEntries,values,indices);
    
      //Insert row 3-----------------------------------------
      values[0] = M(2,0);
      values[1] = M(2,1);
      values[2] = M(2,2);
      values[3] = M(2,3);
    
      invA.InsertGlobalValues(row2,numEntries,values,indices);
    
      //Insert row 3-----------------------------------------
      values[0] = M(3,0);
      values[1] = M(3,1);
      values[2] = M(3,2);
      values[3] = M(3,3);
    
      invA.InsertGlobalValues(row3,numEntries,values,indices);
    }
  }
  
  //Fill complete
  invA.FillComplete(A.DomainMap(),A.RangeMap());

  return(invA);
}

epetraMatrixManip::RCP 
epetraMatrixManip::
inv4x4(const RCP_CONST & A) const
{
  Epetra_CrsMatrix M = inv4x4(*A);
  return(Teuchos::rcp(new Epetra_CrsMatrix(M)));
}


Epetra_CrsMatrix
epetraMatrixManip::
sum(const Epetra_CrsMatrix & A,
    const Epetra_CrsMatrix & B)
{
  Epetra_CrsMatrix * Cptr = NULL;
  
  EpetraExt::MatrixMatrix::Add(A,false,1.0,
                               B,false,1.0,
                               Cptr);
  
  Epetra_Map  rangeMap = A.RangeMap();
  Epetra_Map domainMap = B.DomainMap();
  Cptr->FillComplete(domainMap,rangeMap);
  
  Epetra_CrsMatrix C(*Cptr);
  delete Cptr;
  
  return(C);
}
    
epetraMatrixManip::RCP
epetraMatrixManip::
sum(const RCP_CONST & A,
    const RCP_CONST & B) const
{  
  Epetra_CrsMatrix C = sum(*A,*B);
  return(Teuchos::rcp(new Epetra_CrsMatrix(C)));
  
}
    
Epetra_CrsMatrix
epetraMatrixManip::
sum(const Epetra_CrsMatrix & A,
    const Epetra_CrsMatrix & B,
    const Epetra_CrsMatrix & C)
{
  Epetra_CrsMatrix * D1ptr = NULL;
  Epetra_CrsMatrix * D2ptr = NULL;
  
  EpetraExt::MatrixMatrix::Add(A,false,1.0,
                               B,false,1.0,
                               D1ptr);
  
  Epetra_Map  rangeMap = A.RangeMap();
  Epetra_Map domainMap = B.DomainMap();
  D1ptr->FillComplete(domainMap,rangeMap);
  
  EpetraExt::MatrixMatrix::Add(*D1ptr,false,1.0,
                                C,    false,1.0,
                                D2ptr);
  
  rangeMap = D1ptr->RangeMap();
  domainMap = C.DomainMap();
  D2ptr->FillComplete(domainMap,rangeMap);
  
  Epetra_CrsMatrix D(*D2ptr);
  delete D2ptr;
  
  return(D);
}
    
epetraMatrixManip::RCP
epetraMatrixManip::
sum(const RCP_CONST & A,
    const RCP_CONST & B,
    const RCP_CONST & C) const
{
  Epetra_CrsMatrix D = sum(*A,*B,*C);
  return(Teuchos::rcp(new Epetra_CrsMatrix(D)));
}

Epetra_CrsMatrix
epetraMatrixManip::
matrixLV(const Epetra_CrsMatrix & L,
         const Real & toll)
{
  //Typedefs---------------------------------------------------------
  typedef std::set<int>::iterator SETINT_ITERATOR;
  
  //Maps-------------------------------------------------------------
  Epetra_Map mapL = L.RowMap();
  Epetra_Map mapQ = L.DomainMap();
  
  //Identify the active local Q-rows---------------------------------
  Real maxVal = L.NormInf() * toll;
  int maxSize = L.MaxNumEntries();
  
  int      numEntries, rowG;
  int    * indices;
  double * values;
  
  values  = new double[maxSize];
  indices = new int[maxSize];
  
  std::set<int> localInactiveQ;
  
  for(int rowL = 0; rowL < mapL.NumMyElements(); rowL++)
  { 
    rowG = mapL.GID(rowL);
    L.ExtractGlobalRowCopy(rowG,maxSize,numEntries,values,indices);

    for(UInt i=0; i < numEntries; ++i)
    {        
      if(abs(values[i]) >= maxVal)
      { localInactiveQ.insert(indices[i]); }
    }
  }
  
  delete values;
  delete indices;
  
  assert(mapL.NumMyElements() == localInactiveQ.size());
  
  //Compute global---------------------------------------------------
  sVect<int> inMerge, outMerge;
  int locNumInactiveQ, globNumInactiveQ;
  
  for(SETINT_ITERATOR iter = localInactiveQ.begin(); iter != localInactiveQ.end(); ++iter)
  { inMerge.push_back(*iter); }
  
  locNumInactiveQ = inMerge.size();
  mapQ.Comm().MaxAll(&locNumInactiveQ, &globNumInactiveQ, 1);
  
  inMerge.resize(globNumInactiveQ);
  for(UInt i=(locNumInactiveQ+1); i <= inMerge.size(); ++i)
  { inMerge(i) = -1; }
  
  locNumInactiveQ   = inMerge.size();
  globNumInactiveQ *= mapQ.Comm().NumProc();
  outMerge.resize(globNumInactiveQ);
  
  mapQ.Comm().GatherAll(&inMerge[0],
                        &outMerge[0],
                        locNumInactiveQ);
  
  for(UInt i=1; i <= outMerge.size(); ++i)
  {
    if(outMerge(i) != -1)
    { localInactiveQ.insert(outMerge(i)); }
  }
  
  //Build V matrix---------------------------------------------------
  Epetra_Map mapV(-1, mapQ.NumMyElements() - mapL.NumMyElements(), 0, mapQ.Comm());
  Epetra_CrsMatrix V(Copy,mapQ,1);
  
  int col = 0;
  int indexI, indexJ[1];
  double val[1];
  
  indexJ[0] = 0;
  val[0]    = 1.0;
  
  for(int rowL=0; rowL < V.NumMyRows(); ++rowL)
  {
    rowG = mapQ.GID(rowL);
      
    if(localInactiveQ.count(rowG) == 0)
    {
      indexI    = rowG;
      indexJ[0] = mapV.GID(col);
      
      V.InsertGlobalValues(indexI,1,val,indexJ);
      col++;
    }
  }
  
  V.FillComplete(mapV,mapQ);
  
  return(V);
}

epetraMatrixManip::RCP
epetraMatrixManip::
matrixLV(const RCP_CONST & LV,
         const Real & toll) const
{
  Epetra_CrsMatrix D = matrixLV(*LV,toll);
  return(Teuchos::rcp(new Epetra_CrsMatrix(D)));
}

Epetra_CrsMatrix
epetraMatrixManip::
taylorInv(const Epetra_CrsMatrix & A,
          const UInt & it)
{
  Epetra_CrsMatrix PMP  = diagInv(A);
  Epetra_CrsMatrix invA = diagInv(A);
  
  for(UInt k=1; k <= it; ++k)
  {    
    PMP = multiply(invA,A,invA);
    PMP.Scale(-1.0);
    
    invA.Scale(2.0);
    invA = sum(invA,PMP);
  }
  
  return(invA);
}
    
epetraMatrixManip::RCP
epetraMatrixManip::
taylorInv(const RCP_CONST & A,
          const UInt & it)
{
  Epetra_CrsMatrix invA = taylorInv(*A,it);
  return(Teuchos::rcp(new Epetra_CrsMatrix(invA)));
}
