/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "epetraOperatorPrj.h"
#include "epetraMatrixManip.h"

//_________________________________________________________________________________________________
// CONSTRUCTOR AND FUNCTIONS
//-------------------------------------------------------------------------------------------------
epetraOperatorPrj::
epetraOperatorPrj()
{
  matrixOk = false;
  precOk   = false;
}

epetraOperatorPrj::
epetraOperatorPrj(const Epetra_CrsMatrix & LL)
{
  matrixOk = true;
  precOk   = false;
  
  L = Teuchos::rcp(new Epetra_CrsMatrix(LL));
}

epetraOperatorPrj::
epetraOperatorPrj(const Teuchos::RCP<Epetra_CrsMatrix> & LL)
{
  matrixOk = true;
  precOk   = false;
  
  L = Teuchos::rcp(new Epetra_CrsMatrix(*LL));
}

void
epetraOperatorPrj::
setMatrix(const Epetra_CrsMatrix & LL)
{
  matrixOk = true;
  precOk   = false;
  
  L = Teuchos::rcp(new Epetra_CrsMatrix(LL));
}

void
epetraOperatorPrj::
setMatrix(const Teuchos::RCP<Epetra_CrsMatrix> & LL)
{
  matrixOk = true;
  precOk   = false;
  
  L = LL;
}

void
epetraOperatorPrj::
buildPreconditioner(const Teuchos::ParameterList & solverList,
                    const Teuchos::ParameterList & preconditionerList)
{
  Teuchos::RCP<Teuchos::ParameterList> SolverList(new Teuchos::ParameterList(solverList));
  Teuchos::RCP<Teuchos::ParameterList> PreconditionerList(new Teuchos::ParameterList(preconditionerList));
  
  buildPreconditioner(SolverList,PreconditionerList);
}

void
epetraOperatorPrj::
buildPreconditioner(const Teuchos::RCP<Teuchos::ParameterList> & solverList,
                    const Teuchos::RCP<Teuchos::ParameterList> & preconditionerList)
{
  assert(matrixOk);
  precOk = true;
  
  epetraMatrixManip matManip;
  LT = Teuchos::rcp(new Epetra_CrsMatrix(matManip.transpose(*L)));
  A  = Teuchos::rcp(new Epetra_CrsMatrix(matManip.multiply(*L,*LT)));
  
  FACTORYPREC factoryPrec;
  Teuchos::RCP<PREC> preconditioner = factoryPrec.create(preconditionerList);
  preconditioner->setOperator(A);
  preconditioner->initialize();
  preconditioner->compute();
  belosPrec = preconditioner->getPreconditioner();
  
  FACTORYSOLVER factorySolver;
  linearSolver = factorySolver.create(solverList);
  linearSolver->setRightPrec(belosPrec);
}


//_________________________________________________________________________________________________
// INHERITED FUNCTIONS
//-------------------------------------------------------------------------------------------------
int
epetraOperatorPrj::
SetUseTranspose(bool UseTranspose)
{
  return(0);
}

int
epetraOperatorPrj::
Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  assert(matrixOk);
  assert(precOk);
  
  Epetra_Map mapSmall = L->RangeMap();
  Teuchos::RCP<Epetra_Vector> W(new Epetra_Vector(mapSmall)), Z(new Epetra_Vector(mapSmall));
  
  L->Multiply(false,X,*Z);
  
  linearSolver->setProblem(A,W,Z);
  bool success = linearSolver->solve();
  
  LT->Multiply(false,*W,Y);
  Y.Update(1.0,X,-1.0);
  
  return(0);
}

int
epetraOperatorPrj::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  return(1);
}

double
epetraOperatorPrj::
NormInf() const
{
  return(L->NormInf());
}

const char *
epetraOperatorPrj::
Label() const
{
  return("operatorPrj");
}

bool
epetraOperatorPrj::
UseTranspose() const
{
  return(false);
}

bool
epetraOperatorPrj::
HasNormInf() const
{
  return(true);
}

const Epetra_Comm &
epetraOperatorPrj::
Comm() const
{
  assert(matrixOk);
  return(L->Comm());
}

const Epetra_Map &
epetraOperatorPrj::
OperatorDomainMap() const
{
  assert(matrixOk);
  return(L->OperatorDomainMap());
}

const Epetra_Map &
epetraOperatorPrj::
OperatorRangeMap() const
{
  assert(matrixOk);
  return(L->OperatorDomainMap());
}

