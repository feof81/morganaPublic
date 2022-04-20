/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "epetraOperatorPrecSaddleHD.h"
#include "epetraMatrixManip.h"

epetraOperatorPrecSaddleHD::
epetraOperatorPrecSaddleHD(const communicator & CommDev)
{
  matrixOk = false;
  precOk   = false;
  
  commDev = Teuchos::rcp(new communicator(CommDev));
}

epetraOperatorPrecSaddleHD::
epetraOperatorPrecSaddleHD(const Teuchos::RCP<communicator> & CommDev)

{
  matrixOk = false;
  precOk   = false;
  
  commDev = CommDev;
}

epetraOperatorPrecSaddleHD::
epetraOperatorPrecSaddleHD(const communicator     & CommDev,
                          const Epetra_CrsMatrix & AA,
                          const Epetra_CrsMatrix & LL)
{
  matrixOk = true;
  precOk   = false;
  
  commDev = Teuchos::rcp(new communicator(CommDev));
  
  A = Teuchos::rcp(new Epetra_CrsMatrix(AA));
  L = Teuchos::rcp(new Epetra_CrsMatrix(LL));
}

epetraOperatorPrecSaddleHD::
epetraOperatorPrecSaddleHD(const Teuchos::RCP<communicator>     & CommDev,
                          const Teuchos::RCP<Epetra_CrsMatrix> & AA,
                          const Teuchos::RCP<Epetra_CrsMatrix> & LL)
{
  matrixOk = true;
  precOk   = false;
  
  commDev = CommDev;
  A = AA;
  L = LL;
}

void
epetraOperatorPrecSaddleHD::
setMatrix(const Epetra_CrsMatrix & AA,
          const Epetra_CrsMatrix & LL)
{
  matrixOk = true;
  precOk   = false;
  
  A = Teuchos::rcp(new Epetra_CrsMatrix(AA));
  L = Teuchos::rcp(new Epetra_CrsMatrix(LL));
}
    
void
epetraOperatorPrecSaddleHD::
setMatrix(const Teuchos::RCP<Epetra_CrsMatrix> & AA,
          const Teuchos::RCP<Epetra_CrsMatrix> & LL)
{
  matrixOk = true;
  precOk   = false;
  
  A = AA;
  L = LL;
}
    
void
epetraOperatorPrecSaddleHD::
buildPreconditioner()
{ 
  //Controls---------------------------------------------------------
  assert(matrixOk);
  precOk = true;
  
  //Download maps----------------------------------------------------
  Epetra_MpiComm epetraComm(*commDev);
  rowMapA = Teuchos::rcp(new Epetra_Map(A->RangeMap()));
  rowMapL = Teuchos::rcp(new Epetra_Map(L->RangeMap()));
  
  //Block assembler--------------------------------------------------
  assembler = Teuchos::rcp(new epetraBlockVectorAssemble(2));
  assembler->setEpetraComm(epetraComm);
  assembler->setRowMap(1,*rowMapA);
  assembler->setRowMap(2,*rowMapL);
  assembler->startup();
  
  map = assembler->getRowMap();
  
  //Stripper---------------------------------------------------------
  splitter = Teuchos::rcp(new epetraVectorStrip(2));
  splitter->setEpetraComm(epetraComm);
  splitter->setRowMap(1,*rowMapA);
  splitter->setRowMap(2,*rowMapL);
  splitter->startup();
  
  //Vectors----------------------------------------------------------
  Xa = Teuchos::rcp(new Epetra_MultiVector(*rowMapA,1));
  Xl = Teuchos::rcp(new Epetra_MultiVector(*rowMapL,1));
  Ta = Teuchos::rcp(new Epetra_MultiVector(*rowMapA,1));
  Tl = Teuchos::rcp(new Epetra_MultiVector(*rowMapL,1));
  
  //Matrices---------------------------------------------------------
  epetraMatrixManip matManip;
  LT  = Teuchos::rcp(new Epetra_CrsMatrix(matManip.transpose(*L)));
  LLT = Teuchos::rcp(new Epetra_CrsMatrix(matManip.multiply(*L,*LT)));
  
  //Direct solvers---------------------------------------------------
  problemA = Teuchos::rcp(new Epetra_LinearProblem(A.get(),   Xa.get(), Ta.get()));
  problemL = Teuchos::rcp(new Epetra_LinearProblem(LLT.get(), Xl.get(), Tl.get()));
  
  Amesos Factory;
  char* SolverType = "Amesos_Mumps";
  
  solverA = Teuchos::rcp(Factory.Create(SolverType,*problemA));
  solverL = Teuchos::rcp(Factory.Create(SolverType,*problemL));
  
  /*AMESOS_CHK_ERR(solverA->SymbolicFactorization());
  AMESOS_CHK_ERR(solverL->SymbolicFactorization());
  
  AMESOS_CHK_ERR(solverA->NumericFactorization());
  AMESOS_CHK_ERR(solverL->NumericFactorization());*/
  
  solverA->SymbolicFactorization();
  solverL->SymbolicFactorization();
  
  solverA->NumericFactorization();
  solverL->NumericFactorization();
}
    
int
epetraOperatorPrecSaddleHD::
SetUseTranspose(bool UseTranspose)
{
  return(0);
}

int
epetraOperatorPrecSaddleHD::
Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  //Divide the two sub-vectors---------------------------------------
  Teuchos::RCP<Epetra_MultiVector> Xt(new Epetra_MultiVector(splitter->getVectorRef(X,1)));
  
  //First part-------------------------------------------------------
  *Ta = splitter->getVectorRef(X,1);
  AMESOS_CHK_ERR(solverA->Solve()); //Result in Xa
  
  //Second part------------------------------------------------------
  *Tl = splitter->getVectorRef(X,2);
  AMESOS_CHK_ERR(solverL->Solve());
  
  LT->Multiply(false,*Tl,*Xt);
  A->Multiply(false,*Xt,*Xt);
  L->Multiply(false,*Xt,*Tl);
  
  AMESOS_CHK_ERR(solverL->Solve()); //Result in Xl
  
  //Rebuild the output vector----------------------------------------
  assembler->startup();
  assembler->setVector(*Xa,1);
  assembler->setVector(*Xl,2);
  Teuchos::RCP<Epetra_MultiVector> T = assembler->getVector();
  
  Y = *T;
  
  return(0);
}

int
epetraOperatorPrecSaddleHD::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  return(1);
}

double
epetraOperatorPrecSaddleHD::
NormInf() const
{
  assert(matrixOk);
  return(A->NormInf());
}

const char *
epetraOperatorPrecSaddleHD::
Label() const
{
  return("operatorPrecSaddleHD");
}

bool
epetraOperatorPrecSaddleHD::
UseTranspose() const
{
  return(false);
}

bool
epetraOperatorPrecSaddleHD::
HasNormInf() const
{
  return(true);
}

const Epetra_Comm &
epetraOperatorPrecSaddleHD::
Comm() const
{
  assert(matrixOk);
  return(L->Comm());
}

const Epetra_Map &
epetraOperatorPrecSaddleHD::
OperatorDomainMap() const
{
  assert(precOk);
  return(*map);
}

const Epetra_Map &
epetraOperatorPrecSaddleHD::
OperatorRangeMap() const
{
  assert(precOk);
  return(*map);
}
