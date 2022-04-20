/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "epetraOperatorPrecSaddleB.h"
#include "epetraMatrixManip.h"

epetraOperatorPrecSaddleB::
epetraOperatorPrecSaddleB(const communicator & CommDev)
{
  matrixOk = false;
  precOk   = false;
  
  commDev = Teuchos::rcp(new communicator(CommDev));
}

epetraOperatorPrecSaddleB::
epetraOperatorPrecSaddleB(const Teuchos::RCP<communicator> & CommDev)

{
  matrixOk = false;
  precOk   = false;
  
  commDev = CommDev;
}

epetraOperatorPrecSaddleB::
epetraOperatorPrecSaddleB(const communicator     & CommDev,
                          const Epetra_CrsMatrix & AA,
                          const Epetra_CrsMatrix & LL)
{
  matrixOk = true;
  precOk   = false;
  
  commDev = Teuchos::rcp(new communicator(CommDev));
  A       = Teuchos::rcp(new Epetra_CrsMatrix(AA));
  L       = Teuchos::rcp(new Epetra_CrsMatrix(LL));
}

epetraOperatorPrecSaddleB::
epetraOperatorPrecSaddleB(const Teuchos::RCP<communicator>     & CommDev,
                          const Teuchos::RCP<Epetra_CrsMatrix> & AA,
                          const Teuchos::RCP<Epetra_CrsMatrix> & LL)
{
  matrixOk = true;
  precOk   = false;
  
  commDev = CommDev;
  A       = AA;
  L       = LL;
}

void
epetraOperatorPrecSaddleB::
setMatrix(const Epetra_CrsMatrix & AA,
          const Epetra_CrsMatrix & LL)
{
  matrixOk = true;
  precOk   = false;
  
  A = Teuchos::rcp(new Epetra_CrsMatrix(AA));
  L = Teuchos::rcp(new Epetra_CrsMatrix(LL));
}
    
void
epetraOperatorPrecSaddleB::
setMatrix(const Teuchos::RCP<Epetra_CrsMatrix> & AA,
          const Teuchos::RCP<Epetra_CrsMatrix> & LL)
{
  matrixOk = true;
  precOk   = false;
  
  A = AA;
  L = LL;
}

void
epetraOperatorPrecSaddleB::
buildPreconditioner(const Teuchos::ParameterList & solverList,
                    const Teuchos::ParameterList & preconditionerList)
{
  Teuchos::RCP<Teuchos::ParameterList> SolverList        (new Teuchos::ParameterList(solverList));
  Teuchos::RCP<Teuchos::ParameterList> PreconditionerList(new Teuchos::ParameterList(preconditionerList));
  
  buildPreconditioner(SolverList,PreconditionerList);
}

void
epetraOperatorPrecSaddleB::
buildPreconditioner(const Teuchos::RCP<Teuchos::ParameterList> & solverList,
                    const Teuchos::RCP<Teuchos::ParameterList> & preconditionerList)
{
  //Update list
  solList = solverList;
  
  //Controls
  assert(matrixOk);
  precOk = true;
  
  //Download maps
  Epetra_MpiComm epetraComm(*commDev);
  rowMapA = Teuchos::rcp(new Epetra_Map(A->RangeMap()));
  rowMapL = Teuchos::rcp(new Epetra_Map(L->RangeMap()));
  
  //Block assembler
  assembler = Teuchos::rcp(new epetraBlockVectorAssemble(2));
  assembler->setEpetraComm(epetraComm);
  assembler->setRowMap(1,*rowMapA);
  assembler->setRowMap(2,*rowMapL);
  assembler->startup();
  
  map = assembler->getRowMap();
  
  //Stripper
  splitter = Teuchos::rcp(new epetraVectorStrip(2));
  splitter->setEpetraComm(epetraComm);
  splitter->setRowMap(1,*rowMapA);
  splitter->setRowMap(2,*rowMapL);
  splitter->startup();
  
  //Matrices
  epetraMatrixManip matManip;
  LT  = Teuchos::rcp(new Epetra_CrsMatrix(matManip.transpose(*L)));
  LLT = Teuchos::rcp(new Epetra_CrsMatrix(matManip.multiply(*L,*LT)));
  
  //Preconditioners
  FACTORYPREC factoryPrec;
  
  Teuchos::RCP<PREC> preconditionerA = factoryPrec.create(preconditionerList);
  preconditionerA->setOperator(A);
  preconditionerA->initialize();
  preconditionerA->compute();
  belosPrecA = preconditionerA->getPreconditioner();
  
  Teuchos::RCP<PREC> preconditionerLLT = factoryPrec.create(preconditionerList);
  preconditionerLLT->setOperator(LLT);
  preconditionerLLT->initialize();
  preconditionerLLT->compute();
  belosPrecLLT = preconditionerLLT->getPreconditioner();
}


//_________________________________________________________________________________________________
// INHERITED FUNCTIONS
//-------------------------------------------------------------------------------------------------
int
epetraOperatorPrecSaddleB::
SetUseTranspose(bool UseTranspose)
{
  return(0);
}

int
epetraOperatorPrecSaddleB::
Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  //Typedefs---------------------------------------------------------------------------------------
  typedef factoryLinearSolver<double, Epetra_MultiVector, Epetra_Operator> FACTORYSOLVER;
  typedef virtualLinearSolver<double, Epetra_MultiVector, Epetra_Operator> SOLVER;
  
  FACTORYSOLVER factorySolver;
  Teuchos::RCP<SOLVER> linearSolver = factorySolver.create(solList);
  
  //Divide the two sub-vectors---------------------------------------------------------------------
  Teuchos::RCP<Epetra_MultiVector> xi(new Epetra_MultiVector(splitter->getVectorRef(X,1)));
  Teuchos::RCP<Epetra_MultiVector> yi(new Epetra_MultiVector(splitter->getVectorRef(X,2)));
  
  Teuchos::RCP<Epetra_MultiVector> xo(new Epetra_MultiVector(*xi));
  Teuchos::RCP<Epetra_MultiVector> xt(new Epetra_MultiVector(*xi));
  Teuchos::RCP<Epetra_MultiVector> yo(new Epetra_MultiVector(*yi));
  
  //First part-------------------------------------------------------------------------------------  
  linearSolver->setProblem(A,xo,xi);
  linearSolver->setRightPrec(belosPrecA);
  bool success = linearSolver->solve();

  //Second part------------------------------------------------------------------------------------
  linearSolver->setProblem(LLT,yo,yi);
  linearSolver->setRightPrec(belosPrecLLT);
  success = success && linearSolver->solve();
  
  LT->Multiply(false,*yo,*xi);
  A->Multiply(false,*xi,*xt);
  L->Multiply(false,*xt,*yi);
  
  linearSolver->setProblem(LLT,yo,yi);
  linearSolver->setRightPrec(belosPrecLLT);
  success = success && linearSolver->solve();
  
  //Rebuild the output vector----------------------------------------------------------------------
  assembler->startup();
  assembler->setVector(*xo,1);
  assembler->setVector(*yo,2);
  Teuchos::RCP<Epetra_MultiVector> T = assembler->getVector();
  
  Y = *T;
  
  return(int(!success));
}

int
epetraOperatorPrecSaddleB::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  return(Apply(X,Y));
}

double
epetraOperatorPrecSaddleB::
NormInf() const
{
  assert(matrixOk);
  return(A->NormInf());
}

const char *
epetraOperatorPrecSaddleB::
Label() const
{
  return("operatorPrecSaddleA");
}

bool
epetraOperatorPrecSaddleB::
UseTranspose() const
{
  return(false);
}

bool
epetraOperatorPrecSaddleB::
HasNormInf() const
{
  return(true);
}

const Epetra_Comm &
epetraOperatorPrecSaddleB::
Comm() const
{
  assert(matrixOk);
  return(L->Comm());
}

const Epetra_Map &
epetraOperatorPrecSaddleB::
OperatorDomainMap() const
{
  assert(precOk);
  return(*map);
}

const Epetra_Map &
epetraOperatorPrecSaddleB::
OperatorRangeMap() const
{
  assert(precOk);
  return(*map);
}
