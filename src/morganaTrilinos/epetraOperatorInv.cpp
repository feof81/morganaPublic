/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "epetraOperatorInv.h"


//_________________________________________________________________________________________________
// CONSTRUCTOR AND FUNCTIONS
//-------------------------------------------------------------------------------------------------
epetraOperatorInv::
epetraOperatorInv()
{
  //Control logic
  matrixOk = false;
  solverOk = false;
}

epetraOperatorInv::
epetraOperatorInv(const Teuchos::RCP<Epetra_CrsMatrix> & AA)
{
  //Control logic
  matrixOk = true;
  solverOk = false;
  
  //Set matrix
  A = AA;
}
    
epetraOperatorInv::
epetraOperatorInv(const Teuchos::RCP<Epetra_CrsMatrix> & AA,
                  const Teuchos::ParameterList & SolverList,
                  const Teuchos::ParameterList & PreconditionerList)
{
  //Control logic
  matrixOk = true;
  solverOk = true;
  
  //Set matrix
  A = AA;
  
  //Set parameters
  solverList = SolverList;
  preconList = PreconditionerList;
  
  //Build preconditioner
  FACTORYPREC factoryPrec;
  preconditioner = factoryPrec.create(Teuchos::rcp(new Teuchos::ParameterList(preconList)));
  preconditioner->setOperator(Teuchos::rcp(new Epetra_CrsMatrix(*A)));
  preconditioner->initialize();
  preconditioner->compute();
  belosPrec = preconditioner->getPreconditioner();
  
  //Build solver
  FACTORYSOLVER factorySolver;
  linearSolver = factorySolver.create(Teuchos::rcp(new Teuchos::ParameterList(solverList)));
  linearSolver->setRightPrec(belosPrec);
}
    
epetraOperatorInv::
epetraOperatorInv(const epetraOperatorInv & op)
{
  //Control logic
  matrixOk = true;
  solverOk = true;
  
  //Set matrix
  A = op.A;
  
  //Set parameters
  solverList = op.solverList;
  preconList = op.preconList;
  
  //Build preconditioner
  FACTORYPREC factoryPrec;
  Teuchos::RCP<PREC> preconditioner = factoryPrec.create(Teuchos::rcp(new Teuchos::ParameterList(preconList)));
  preconditioner->setOperator(Teuchos::rcp(new Epetra_CrsMatrix(*A)));
  preconditioner->initialize();
  preconditioner->compute();
  belosPrec = preconditioner->getPreconditioner();
  
  //Build solver
  FACTORYSOLVER factorySolver;
  linearSolver = factorySolver.create(Teuchos::rcp(new Teuchos::ParameterList(solverList)));
  linearSolver->setRightPrec(belosPrec);
}    

epetraOperatorInv
epetraOperatorInv::
operator=(const epetraOperatorInv & op)
{
  //Control logic
  matrixOk = true;
  solverOk = true;
  
  //Set matrix
  A = op.A;
  
  //Set parameters
  solverList = op.solverList;
  preconList = op.preconList;
  
  //Build preconditioner
  FACTORYPREC factoryPrec;
  Teuchos::RCP<PREC> preconditioner = factoryPrec.create(Teuchos::rcp(new Teuchos::ParameterList(preconList)));
  preconditioner->setOperator(Teuchos::rcp(new Epetra_CrsMatrix(*A)));
  preconditioner->initialize();
  preconditioner->compute();
  belosPrec = preconditioner->getPreconditioner();
  
  //Build solver
  FACTORYSOLVER factorySolver;
  linearSolver = factorySolver.create(Teuchos::rcp(new Teuchos::ParameterList(solverList)));
  linearSolver->setRightPrec(belosPrec);
  
  return(*this);
}    

void
epetraOperatorInv::
setOperator(const Teuchos::RCP<Epetra_CrsMatrix> & AA)
{
  //Control logic
  matrixOk = true;
  
  //Set matrix
  A = AA;
}

void
epetraOperatorInv::
setSolver(const Teuchos::ParameterList & SolverList,
          const Teuchos::ParameterList & PreconditionerList)
{
  //Assert logic
  assert(matrixOk);
  
  //Control logic
  solverOk = true;
  
  //Set parameters
  solverList = SolverList;
  preconList = PreconditionerList;
  
  //Build preconditioner
  FACTORYPREC factoryPrec;
  preconditioner = factoryPrec.create(Teuchos::rcp(new Teuchos::ParameterList(preconList)));
  preconditioner->setOperator(Teuchos::rcp(new Epetra_CrsMatrix(*A)));
  preconditioner->initialize();
  preconditioner->compute();
  belosPrec = preconditioner->getPreconditioner();
  
  //Build solver
  FACTORYSOLVER factorySolver;
  linearSolver = factorySolver.create(Teuchos::rcp(new Teuchos::ParameterList(solverList)));
  linearSolver->setRightPrec(belosPrec);
}

Teuchos::RCP<epetraOperatorInv>
epetraOperatorInv::
op(const Teuchos::RCP<Epetra_CrsMatrix> & AA)
{
  return(Teuchos::rcp(new epetraOperatorInv(AA)));
}



//_________________________________________________________________________________________________
// INHERITED FUNCTIONS
//-------------------------------------------------------------------------------------------------
int
epetraOperatorInv::
SetUseTranspose(bool UseTranspose)
{
  //Assert logic
  assert(matrixOk);
  assert(solverOk);
  
  return(0);
}

int
epetraOperatorInv::
Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  //Assert logic
  assert(matrixOk);
  assert(solverOk);
  
  assert( X.Map().SameAs(A->OperatorRangeMap())  );
  assert( Y.Map().SameAs(A->OperatorDomainMap()) );
  
  //Linear system
  Teuchos::RCP<Epetra_MultiVector> U(new Epetra_MultiVector(Y));
  Teuchos::RCP<Epetra_MultiVector> T(new Epetra_MultiVector(X));
  
  linearSolver->setProblem(Teuchos::rcp(new Epetra_CrsMatrix(*A)), U, T);
  bool success = linearSolver->solve();
  
  Y = *U;
  
  return(int(!success));
}

int
epetraOperatorInv::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  //Assert logic
  assert(matrixOk);
  assert(solverOk);
  
  return(0);
}

double
epetraOperatorInv::
NormInf() const
{
  return(A->NormInf());
}

const char *
epetraOperatorInv::
Label() const
{
  return("operatorInv");
}

bool
epetraOperatorInv::
UseTranspose() const
{
  return(false);
}

bool
epetraOperatorInv::
HasNormInf() const
{
  return(false);
}

const Epetra_Comm &
epetraOperatorInv::
Comm() const
{
  return(A->Comm());
}

const Epetra_Map &
epetraOperatorInv::
OperatorDomainMap() const
{
  return(A->OperatorRangeMap());
}

const Epetra_Map &
epetraOperatorInv::
OperatorRangeMap() const
{
  return(A->OperatorDomainMap());
}
