/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "epetraOperatorPrecSaddleA.h"
#include "epetraMatrixManip.h"

epetraOperatorPrecSaddleA::
epetraOperatorPrecSaddleA(const communicator & CommDev)
{
  matrixOk = false;
  precOk   = false;
  
  commDev = Teuchos::rcp(new communicator(CommDev));
}

epetraOperatorPrecSaddleA::
epetraOperatorPrecSaddleA(const Teuchos::RCP<communicator> & CommDev)

{
  matrixOk = false;
  precOk   = false;
  
  commDev = CommDev;
}

epetraOperatorPrecSaddleA::
epetraOperatorPrecSaddleA(const communicator     & CommDev,
                          const Epetra_CrsMatrix & AA,
                          const Epetra_CrsMatrix & LL)
{
  matrixOk = true;
  precOk   = false;
  
  commDev = Teuchos::rcp(new communicator(CommDev));
  A       = Teuchos::rcp(new Epetra_CrsMatrix(AA));
  L       = Teuchos::rcp(new Epetra_CrsMatrix(LL));
}

epetraOperatorPrecSaddleA::
epetraOperatorPrecSaddleA(const Teuchos::RCP<communicator>     & CommDev,
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
epetraOperatorPrecSaddleA::
setMatrix(const Epetra_CrsMatrix & AA,
          const Epetra_CrsMatrix & LL)
{
  matrixOk = true;
  precOk   = false;
  
  A = Teuchos::rcp(new Epetra_CrsMatrix(AA));
  L = Teuchos::rcp(new Epetra_CrsMatrix(LL));
}
    
void
epetraOperatorPrecSaddleA::
setMatrix(const Teuchos::RCP<Epetra_CrsMatrix> & AA,
          const Teuchos::RCP<Epetra_CrsMatrix> & LL)
{
  matrixOk = true;
  precOk   = false;
  
  A = AA;
  L = LL;
}

void
epetraOperatorPrecSaddleA::
buildPreconditioner(const Teuchos::ParameterList & preconditionerList)
{
  Teuchos::RCP<Teuchos::ParameterList> PreconditionerList(new Teuchos::ParameterList(preconditionerList));
  buildPreconditioner(PreconditionerList);
}

void
epetraOperatorPrecSaddleA::
buildPreconditioner(const Teuchos::RCP<Teuchos::ParameterList> & preconditionerList)
{
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
epetraOperatorPrecSaddleA::
SetUseTranspose(bool UseTranspose)
{
  return(0);
}

int
epetraOperatorPrecSaddleA::
Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  //Divide the two sub-vectors
  Epetra_MultiVector xi = splitter->getVectorRef(X,1);
  Epetra_MultiVector yi = splitter->getVectorRef(X,2);
  
  Epetra_MultiVector xo(xi), xt(xi);
  Epetra_MultiVector yo(yi);
  
  //First part
  belosPrecA->Apply(xi,xo);

  //Second part
  belosPrecLLT->Apply(yi,yo);
  LT->Multiply(false,yo,xi);
  A->Multiply(false,xi,xt);
  L->Multiply(false,xt,yi);  
  belosPrecLLT->Apply(yi,yo);
  
  //Rebuild the output vector
  assembler->startup();
  assembler->setVector(xo,1);
  assembler->setVector(yo,2);
  Teuchos::RCP<Epetra_MultiVector> T = assembler->getVector();
  
  Y = *T; 
  
  return(0);
}

int
epetraOperatorPrecSaddleA::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  //Divide the two sub-vectors
  Epetra_MultiVector xi = splitter->getVectorRef(X,1);
  Epetra_MultiVector yi = splitter->getVectorRef(X,2);
  
  Epetra_MultiVector xo(xi), xt(xi);
  Epetra_MultiVector yo(yi);
  
  //First part
  belosPrecA->ApplyInverse(xi,xo);

  //Second part
  belosPrecLLT->ApplyInverse(yi,yo);
  LT->Multiply(false,yo,xi);
  A->Multiply(false,xi,xt);
  L->Multiply(false,xt,yi);  
  belosPrecLLT->ApplyInverse(yi,yo);
  
  //Rebuild the output vector
  assembler->startup();
  assembler->setVector(xo,1);
  assembler->setVector(yo,2);
  Teuchos::RCP<Epetra_MultiVector> T = assembler->getVector();
  
  Y = *T; 
  
  return(0);
}

double
epetraOperatorPrecSaddleA::
NormInf() const
{
  assert(matrixOk);
  return(A->NormInf());
}

const char *
epetraOperatorPrecSaddleA::
Label() const
{
  return("operatorPrecSaddleA");
}

bool
epetraOperatorPrecSaddleA::
UseTranspose() const
{
  return(false);
}

bool
epetraOperatorPrecSaddleA::
HasNormInf() const
{
  return(true);
}

const Epetra_Comm &
epetraOperatorPrecSaddleA::
Comm() const
{
  assert(matrixOk);
  return(L->Comm());
}

const Epetra_Map &
epetraOperatorPrecSaddleA::
OperatorDomainMap() const
{
  assert(precOk);
  return(*map);
}

const Epetra_Map &
epetraOperatorPrecSaddleA::
OperatorRangeMap() const
{
  assert(precOk);
  return(*map);
}
