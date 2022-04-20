/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "epetraOperatorLcomb.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS AND FUNCTIONS
//-------------------------------------------------------------------------------------------------
epetraOperatorLcomb::
epetraOperatorLcomb()
{
  matrixOk = false;
}

epetraOperatorLcomb::
epetraOperatorLcomb(const Real                          & aa,
                    const Teuchos::RCP<Epetra_Operator> & AA,
                    const Real                          & bb,
                    const Teuchos::RCP<Epetra_Operator> & BB)
{
  //Assert
  assert(AA->OperatorDomainMap().SameAs(BB->OperatorDomainMap()));
  assert(AA->OperatorRangeMap().SameAs(BB->OperatorRangeMap()));
  
  //Control logic
  matrixOk = true;
  
  //Copy data
  a = aa;
  b = bb;
  
  A = AA;
  B = BB;
}

epetraOperatorLcomb::
epetraOperatorLcomb(const epetraOperatorLcomb & op)
{
  //Assert
  assert(op.A->OperatorDomainMap().SameAs(op.B->OperatorDomainMap()));
  assert(op.A->OperatorRangeMap().SameAs(op.B->OperatorRangeMap()));
  
  //Control logic
  matrixOk = op.matrixOk;
  
  //Copy data
  a = op.a;
  b = op.b;
  
  A = op.A;
  B = op.B;
}    

epetraOperatorLcomb
epetraOperatorLcomb::
operator=(const epetraOperatorLcomb & op)
{
  //Assert
  assert(op.A->OperatorDomainMap().SameAs(op.B->OperatorDomainMap()));
  assert(op.A->OperatorRangeMap().SameAs(op.B->OperatorRangeMap()));
  
  //Control logic
  matrixOk = op.matrixOk;
  
  //Copy data
  a = op.a;
  b = op.b;
  
  A = op.A;
  B = op.B;
  
  //Return
  return(*this);
}

void
epetraOperatorLcomb::
setOperators(const Real                          & aa,
             const Teuchos::RCP<Epetra_Operator> & AA,
             const Real                          & bb,
             const Teuchos::RCP<Epetra_Operator> & BB)
{
  //Assert
  assert(AA->OperatorDomainMap().SameAs(BB->OperatorDomainMap()));
  assert(AA->OperatorRangeMap().SameAs(BB->OperatorRangeMap()));
  
  //Control logic
  matrixOk = true;
  
  //Copy data
  a = aa;
  b = bb;
  
  A = AA;
  B = BB;
}

Teuchos::RCP<epetraOperatorLcomb>
epetraOperatorLcomb::
op(const Real                          & aa,
   const Teuchos::RCP<Epetra_Operator> & AA,
   const Real                          & bb,
   const Teuchos::RCP<Epetra_Operator> & BB)
{
  return(Teuchos::rcp(new epetraOperatorLcomb(aa,AA,bb,BB)));
}


//_________________________________________________________________________________________________
// INHERITED FUNCTIONS
//-------------------------------------------------------------------------------------------------
int
epetraOperatorLcomb::
SetUseTranspose(bool UseTranspose)
{
  assert(matrixOk);
  
  return(false);
}

int
epetraOperatorLcomb::
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
  
  Y.Update(b,Yb,a);
  
  return(0);
}

int
epetraOperatorLcomb::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  assert(matrixOk);
  assert(false);
  return(1);
}

double
epetraOperatorLcomb::
NormInf() const
{
  assert(matrixOk);
  
  return(a * A->NormInf() + b * B->NormInf());
}

const char *
epetraOperatorLcomb::
Label() const
{
  assert(matrixOk);
  
  return("operatorLcomb");
}

bool
epetraOperatorLcomb::
UseTranspose() const
{
  assert(matrixOk);
  
  return(false);
}

bool
epetraOperatorLcomb::
HasNormInf() const
{
  assert(matrixOk);
  
  return(A->HasNormInf() && B->HasNormInf());
}

const Epetra_Comm &
epetraOperatorLcomb::
Comm() const
{
  assert(matrixOk);
  
  return(A->Comm());
}

const Epetra_Map  &
epetraOperatorLcomb::
OperatorDomainMap() const
{
  assert(matrixOk);
  
  return(A->OperatorDomainMap());
}

const Epetra_Map  &
epetraOperatorLcomb::
OperatorRangeMap() const
{
  assert(matrixOk);
  
  return(A->OperatorRangeMap());
}
