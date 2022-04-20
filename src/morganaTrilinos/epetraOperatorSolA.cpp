/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "epetraOperatorSolA.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS AND FUNCTIONS
//-------------------------------------------------------------------------------------------------
epetraOperatorSolA::
epetraOperatorSolA()
{
  matrixOk = false;
}

epetraOperatorSolA::
epetraOperatorSolA(const Epetra_Operator    & AA,
                   const Epetra_MultiVector & VV)
{
  //Control logic
  matrixOk = true;
  
  //Copy data
  A = &AA;
  V = &VV;
}

epetraOperatorSolA::
epetraOperatorSolA(const epetraOperatorSolA & op)
{
  //Control logic
  matrixOk = op.matrixOk;
  
  //Copy data
  A = op.A;
  V = op.V;
}    

epetraOperatorSolA
epetraOperatorSolA::
operator=(const epetraOperatorSolA & op)
{
  //Control logic
  matrixOk = op.matrixOk;
  
  //Copy data
  A = op.A;
  V = op.V;
  
  //Return
  return(*this);
}

void
epetraOperatorSolA::
setOperators(const Epetra_Operator    & AA,
             const Epetra_MultiVector & VV)
{
  //Control logic
  matrixOk = true;
  
  //Copy data
  A = &AA;
  V = &VV;
}


//_________________________________________________________________________________________________
// INHERITED FUNCTIONS
//-------------------------------------------------------------------------------------------------
int
epetraOperatorSolA::
SetUseTranspose(bool UseTranspose)
{
  assert(matrixOk);
  
  return(false);
}

int
epetraOperatorSolA::
Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  //Control logic
  assert(matrixOk);
 
  //Sum 
  Real v;
  V->Dot(X,&v);
  A->Apply(X,Y);
  Y.Update(v,*V,1.0);
  
  return(0);
}

int
epetraOperatorSolA::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  assert(matrixOk);
  assert(false);
  return(1);
}

double
epetraOperatorSolA::
NormInf() const
{
  assert(matrixOk);
  
  Real normV;
  V->NormInf(&normV);
  
  return(A->NormInf() + normV * normV);
}

const char *
epetraOperatorSolA::
Label() const
{
  assert(matrixOk);
  
  return("operatorSolA");
}

bool
epetraOperatorSolA::
UseTranspose() const
{
  assert(matrixOk);
  
  return(false);
}

bool
epetraOperatorSolA::
HasNormInf() const
{
  assert(matrixOk);
  
  return(A->HasNormInf());
}

const Epetra_Comm &
epetraOperatorSolA::
Comm() const
{
  assert(matrixOk);
  
  return(A->Comm());
}

const Epetra_Map  &
epetraOperatorSolA::
OperatorDomainMap() const
{
  assert(matrixOk);
  
  return(A->OperatorDomainMap());
}

const Epetra_Map  &
epetraOperatorSolA::
OperatorRangeMap() const
{
  assert(matrixOk);
  
  return(A->OperatorRangeMap());
}
