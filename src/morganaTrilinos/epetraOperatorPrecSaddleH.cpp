/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "epetraOperatorPrecSaddleH.h"
#include "epetraMatrixManip.h"

epetraOperatorPrecSaddleH::
epetraOperatorPrecSaddleH(const communicator & CommDev) : sweepCoeff(-1.0)
{
  matrixOk = false;
  precOk   = false;
  
  commDev = Teuchos::rcp(new communicator(CommDev));
}

epetraOperatorPrecSaddleH::
epetraOperatorPrecSaddleH(const Teuchos::RCP<communicator> & CommDev) : sweepCoeff(-1.0)

{
  matrixOk = false;
  precOk   = false;
  
  commDev = CommDev;
}

epetraOperatorPrecSaddleH::
epetraOperatorPrecSaddleH(const communicator     & CommDev,
                          const Epetra_CrsMatrix & AA,
                          const Epetra_CrsMatrix & LL,
                          const Epetra_CrsMatrix & KK,
                          const Epetra_CrsMatrix & MM) : sweepCoeff(-1.0)
{
  matrixOk = true;
  precOk   = false;
  
  commDev = Teuchos::rcp(new communicator(CommDev));
  
  A = Teuchos::rcp(new Epetra_CrsMatrix(AA));
  L = Teuchos::rcp(new Epetra_CrsMatrix(LL));
  K = Teuchos::rcp(new Epetra_CrsMatrix(KK));
  M = Teuchos::rcp(new Epetra_CrsMatrix(MM));
}

epetraOperatorPrecSaddleH::
epetraOperatorPrecSaddleH(const Teuchos::RCP<communicator>     & CommDev,
                          const Teuchos::RCP<Epetra_CrsMatrix> & AA,
                          const Teuchos::RCP<Epetra_CrsMatrix> & LL,
                          const Teuchos::RCP<Epetra_CrsMatrix> & KK,
                          const Teuchos::RCP<Epetra_CrsMatrix> & MM) : sweepCoeff(-1.0)
{
  matrixOk = true;
  precOk   = false;
  
  commDev = CommDev;
  A = AA;
  L = LL;
  K = KK;
  M = MM;
}

void
epetraOperatorPrecSaddleH::
setMatrix(const Epetra_CrsMatrix & AA,
          const Epetra_CrsMatrix & LL,
          const Epetra_CrsMatrix & KK,
          const Epetra_CrsMatrix & MM)
{
  matrixOk = true;
  precOk   = false;
  
  A = Teuchos::rcp(new Epetra_CrsMatrix(AA));
  L = Teuchos::rcp(new Epetra_CrsMatrix(LL));
  K = Teuchos::rcp(new Epetra_CrsMatrix(KK));
  M = Teuchos::rcp(new Epetra_CrsMatrix(MM));
}
    
void
epetraOperatorPrecSaddleH::
setMatrix(const Teuchos::RCP<Epetra_CrsMatrix> & AA,
          const Teuchos::RCP<Epetra_CrsMatrix> & LL,
          const Teuchos::RCP<Epetra_CrsMatrix> & KK,
          const Teuchos::RCP<Epetra_CrsMatrix> & MM)
{
  matrixOk = true;
  precOk   = false;
  
  A = AA;
  L = LL;
  K = KK;
  M = MM;
}

void
epetraOperatorPrecSaddleH::
setSweepCoeff(const Real & Sweep)
{
  sweepCoeff = Sweep;
}

void
epetraOperatorPrecSaddleH::
buildPreconditioner(const Teuchos::ParameterList & PreconditionerList,
                    const Teuchos::ParameterList & SolverList)
{
  Teuchos::RCP<Teuchos::ParameterList> PrecList(new Teuchos::ParameterList(PreconditionerList));
  Teuchos::RCP<Teuchos::ParameterList>  SolList(new Teuchos::ParameterList(SolverList));
  
  buildPreconditioner(PrecList,SolList);
}

void
epetraOperatorPrecSaddleH::
buildPreconditioner(const Teuchos::RCP<Teuchos::ParameterList> & PreconditionerList,
                    const Teuchos::RCP<Teuchos::ParameterList> & SolverList)
{
  //Copy solverList--------------------------------------------------
  solverList = SolverList;
  
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
  
  //Matrices---------------------------------------------------------
  epetraMatrixManip matManip;
  LT  = Teuchos::rcp(new Epetra_CrsMatrix(matManip.transpose(*L)));
  LLT = Teuchos::rcp(new Epetra_CrsMatrix(matManip.multiply(*L,*LT)));
  
  //Build preconditioner P
  P = Teuchos::rcp(new Epetra_CrsMatrix(*K));
  EpetraExt::MatrixMatrix::Add(*M,false,sweepCoeff,*P,1.0);
  
  //Preconditioners--------------------------------------------------
  FACTORYPREC factoryPrec;
  
  Teuchos::RCP<PREC> preconditionerA = factoryPrec.create(PreconditionerList);
  preconditionerA->setOperator(P);
  preconditionerA->initialize();
  preconditionerA->compute();
  belosPrecA = preconditionerA->getPreconditioner();
  
  Teuchos::RCP<PREC> preconditionerLLT = factoryPrec.create(PreconditionerList);
  preconditionerLLT->setOperator(LLT);
  preconditionerLLT->initialize();
  preconditionerLLT->compute();
  belosPrecLLT = preconditionerLLT->getPreconditioner();
}


//_________________________________________________________________________________________________
// INHERITED FUNCTIONS
//-------------------------------------------------------------------------------------------------
int
epetraOperatorPrecSaddleH::
SetUseTranspose(bool UseTranspose)
{
  return(0);
}

int
epetraOperatorPrecSaddleH::
Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  //Divide the two sub-vectors---------------------------------------
  Teuchos::RCP<Epetra_MultiVector> xi(new Epetra_MultiVector(splitter->getVectorRef(X,1)));
  Teuchos::RCP<Epetra_MultiVector> yi(new Epetra_MultiVector(splitter->getVectorRef(X,2)));

  Teuchos::RCP<Epetra_MultiVector> xo(new Epetra_MultiVector(*xi));
  Teuchos::RCP<Epetra_MultiVector> xt(new Epetra_MultiVector(*xi));
  Teuchos::RCP<Epetra_MultiVector> yo(new Epetra_MultiVector(*yi));
  
  //First part-------------------------------------------------------
  typedef virtualLinearSolver<double, Epetra_MultiVector, Epetra_Operator> SOLVER;
  
  FACTORYSOLVER factorySolver;
  Teuchos::RCP<SOLVER> linearSolver = factorySolver.create(solverList);
  linearSolver->setProblem(A,xo,xi);
  linearSolver->setRightPrec(belosPrecA);
  bool success = linearSolver->solve();
  
  //Second part------------------------------------------------------
  linearSolver->setProblem(LLT,yo,yi);
  linearSolver->setRightPrec(belosPrecLLT);
  success = success && linearSolver->solve();
  
  LT->Multiply(false,*yo,*xi);
  A->Multiply(false,*xi,*xt);
  L->Multiply(false,*xt,*yi);
  
  linearSolver->setProblem(LLT,yo,yi);
  linearSolver->setRightPrec(belosPrecLLT);
  success = success && linearSolver->solve();
  
  //Rebuild the output vector----------------------------------------
  assembler->startup();
  assembler->setVector(*xo,1);
  assembler->setVector(*yo,2);
  Teuchos::RCP<Epetra_MultiVector> T = assembler->getVector();
  
  Y = *T; 
  
  return(success == false);
}

int
epetraOperatorPrecSaddleH::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{  
  return(1);
}

double
epetraOperatorPrecSaddleH::
NormInf() const
{
  assert(matrixOk);
  return(A->NormInf());
}

const char *
epetraOperatorPrecSaddleH::
Label() const
{
  return("operatorPrecSaddleA");
}

bool
epetraOperatorPrecSaddleH::
UseTranspose() const
{
  return(false);
}

bool
epetraOperatorPrecSaddleH::
HasNormInf() const
{
  return(true);
}

const Epetra_Comm &
epetraOperatorPrecSaddleH::
Comm() const
{
  assert(matrixOk);
  return(L->Comm());
}

const Epetra_Map &
epetraOperatorPrecSaddleH::
OperatorDomainMap() const
{
  assert(precOk);
  return(*map);
}

const Epetra_Map &
epetraOperatorPrecSaddleH::
OperatorRangeMap() const
{
  assert(precOk);
  return(*map);
}
