/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "epetraOperatorMult.h"

//_________________________________________________________________________________________________
// CONSTRUCTORS AND FUNCTIONS
//-------------------------------------------------------------------------------------------------
epetraOperatorMult::
epetraOperatorMult(const Real & C)
{
  matrixOk = false;
  c = C;
}

epetraOperatorMult::
epetraOperatorMult(const Teuchos::RCP<Epetra_Operator> & AA,
                   const Teuchos::RCP<Epetra_Operator> & BB,
                   const Real & C)
{
  //Assert
  assert(AA->OperatorDomainMap().SameAs(BB->OperatorRangeMap()));
  
  //Control logic
  matrixOk = true;
  
  //Copy data
  A = AA;
  B = BB;
  c = C;
}
    
epetraOperatorMult::
epetraOperatorMult(const epetraOperatorMult & op)
{
  //Assert
  assert(op.A->OperatorDomainMap().SameAs(op.B->OperatorRangeMap()));
  
  //Control logic
  matrixOk = op.matrixOk;
  
  //Copy data
  A = op.A;
  B = op.B;
  c = op.c;
}

epetraOperatorMult
epetraOperatorMult::
operator=(const epetraOperatorMult & op)
{
  //Assert
  assert(op.A->OperatorDomainMap().SameAs(op.B->OperatorRangeMap()));
  
  //Control logic
  matrixOk = op.matrixOk;
  
  //Copy data
  A = op.A;
  B = op.B;
  c = op.c;
  
  //Return
  return(*this);
}    

void
epetraOperatorMult::
setOperators(const Teuchos::RCP<Epetra_Operator> & AA,
             const Teuchos::RCP<Epetra_Operator> & BB,
             const Real & C)
{
  //Assert
  assert(AA->OperatorDomainMap().SameAs(BB->OperatorRangeMap()));
  
  //Control logic
  matrixOk = true;
  
  //Copy data
  A = AA;
  B = BB;
  c = C;
}

Teuchos::RCP<epetraOperatorMult>
epetraOperatorMult::
op(const Teuchos::RCP<Epetra_Operator> & AA,
   const Teuchos::RCP<Epetra_Operator> & BB,
   const Real & C)
{
  return(Teuchos::rcp(new epetraOperatorMult(AA,BB,C)));
}

//_________________________________________________________________________________________________
// INHERITED FUNCTIONS
//-------------------------------------------------------------------------------------------------  
int
epetraOperatorMult::
SetUseTranspose(bool UseTranspose)
{
  assert(matrixOk);
  
  return(false);
}

int
epetraOperatorMult::
Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  //Control logic
  assert(matrixOk);
  assert(A->OperatorDomainMap().SameAs(B->OperatorRangeMap()));
 
  //Sum
  Epetra_MultiVector Z(A->OperatorDomainMap(),1);
  
  B->Apply(X,Z);
  A->Apply(Z,Y);
  
  Y.Scale(c);
  
  //Return
  return(0);
}

int
epetraOperatorMult::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  assert(matrixOk);
  assert(false);
  return(1);
}

double
epetraOperatorMult::
NormInf() const
{
  assert(matrixOk);
  
  return(A->NormInf() * B->NormInf());
}

const char *
epetraOperatorMult::
Label() const
{
  assert(matrixOk);
  
  return("operatorMult");
}

bool
epetraOperatorMult::
UseTranspose() const
{
  assert(matrixOk);
  
  return(false);
}

bool
epetraOperatorMult::
HasNormInf() const
{
  assert(matrixOk);
  
  return(A->HasNormInf() && B->HasNormInf());
}

const Epetra_Comm &
epetraOperatorMult::
Comm() const
{
  assert(matrixOk);
  
  return(A->Comm());
}

const Epetra_Map &
epetraOperatorMult::
OperatorDomainMap() const
{
  assert(matrixOk);
  
  return(B->OperatorDomainMap());
}

const Epetra_Map &
epetraOperatorMult::
OperatorRangeMap() const
{
  assert(matrixOk);
  
  return(A->OperatorRangeMap());
}
