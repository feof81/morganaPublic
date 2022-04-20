/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "epetraOperatorSum.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS AND FUNCTIONS
//-------------------------------------------------------------------------------------------------

epetraOperatorSum::
epetraOperatorSum()
{
  matrixOk = false;
}

epetraOperatorSum::
epetraOperatorSum(const Teuchos::RCP<Epetra_Operator> & AA,
                  const Teuchos::RCP<Epetra_Operator> & BB)
{
  //Assert
  assert(AA->OperatorDomainMap().SameAs(BB->OperatorDomainMap()));
  assert(AA->OperatorRangeMap().SameAs(BB->OperatorRangeMap()));
  
  //Control logic
  matrixOk = true;
  
  //Copy data
  A = AA;
  B = BB;
}

epetraOperatorSum::
epetraOperatorSum(const epetraOperatorSum & op)
{
  //Assert
  assert(op.A->OperatorDomainMap().SameAs(op.B->OperatorDomainMap()));
  assert(op.A->OperatorRangeMap().SameAs(op.B->OperatorRangeMap()));
  
  //Control logic
  matrixOk = op.matrixOk;
  
  //Copy data
  A = op.A;
  B = op.B;
}    

epetraOperatorSum
epetraOperatorSum::
operator=(const epetraOperatorSum & op)
{
  //Assert
  assert(op.A->OperatorDomainMap().SameAs(op.B->OperatorDomainMap()));
  assert(op.A->OperatorRangeMap().SameAs(op.B->OperatorRangeMap()));
  
  //Control logic
  matrixOk = op.matrixOk;
  
  //Copy data
  A = op.A;
  B = op.B;
  
  //Return
  return(*this);
}

void
epetraOperatorSum::
setOperators(const Teuchos::RCP<Epetra_Operator> & AA,
             const Teuchos::RCP<Epetra_Operator> & BB)
{
  //Assert
  assert(AA->OperatorDomainMap().SameAs(BB->OperatorDomainMap()));
  assert(AA->OperatorRangeMap().SameAs(BB->OperatorRangeMap()));
  
  //Control logic
  matrixOk = true;
  
  //Copy data
  A = AA;
  B = BB;
}

Teuchos::RCP<epetraOperatorSum>
epetraOperatorSum::
op(const Teuchos::RCP<Epetra_Operator> & AA,
   const Teuchos::RCP<Epetra_Operator> & BB)
{
  return(Teuchos::rcp(new epetraOperatorSum(AA,BB)));
}


//_________________________________________________________________________________________________
// INHERITED FUNCTIONS
//-------------------------------------------------------------------------------------------------
int
epetraOperatorSum::
SetUseTranspose(bool UseTranspose)
{
  assert(matrixOk);
  
  return(false);
}

int
epetraOperatorSum::
Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  //Control logic
  assert(matrixOk);
  assert( X.Map().SameAs(B->OperatorDomainMap()));
  assert( Y.Map().SameAs(B->OperatorRangeMap()));
 
  //Sum 
  Epetra_MultiVector Yb(Y);
  
  A->Apply(X,Y);
  B->Apply(X,Yb);
  
  Y.Update(1.0,Yb,1.0);
  
  return(0);
}

int
epetraOperatorSum::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  assert(matrixOk);
  assert(false);
  return(1);
}

double
epetraOperatorSum::
NormInf() const
{
  assert(matrixOk);
  
  return(A->NormInf() + B->NormInf());
}

const char *
epetraOperatorSum::
Label() const
{
  assert(matrixOk);
  
  return("operatorSum");
}

bool
epetraOperatorSum::
UseTranspose() const
{
  assert(matrixOk);
  
  return(false);
}

bool
epetraOperatorSum::
HasNormInf() const
{
  assert(matrixOk);
  
  return(A->HasNormInf() && B->HasNormInf());
}

const Epetra_Comm &
epetraOperatorSum::
Comm() const
{
  assert(matrixOk);
  
  return(A->Comm());
}

const Epetra_Map  &
epetraOperatorSum::
OperatorDomainMap() const
{
  assert(matrixOk);
  
  return(A->OperatorDomainMap());
}

const Epetra_Map  &
epetraOperatorSum::
OperatorRangeMap() const
{
  assert(matrixOk);
  
  return(A->OperatorRangeMap());
}
