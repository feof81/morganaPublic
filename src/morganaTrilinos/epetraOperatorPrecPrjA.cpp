/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "epetraOperatorPrecPrjA.h"
#include <boost/concept_check.hpp>


//_________________________________________________________________________________________________
// CONSTRUCTOR AND FUNCTIONS
//-------------------------------------------------------------------------------------------------
epetraOperatorPrecPrjA::
epetraOperatorPrecPrjA()
{
  vectOk = false;
  precOk = false;
  b      = 1.0;
}

void
epetraOperatorPrecPrjA::
setVector(const Real & B, const Teuchos::RCP<VECTOR> & VectorV)
{
  vectOk  = true;
  vectorV = VectorV;
  b       = B;
}

void
epetraOperatorPrecPrjA::
setPreconditioner(const Teuchos::RCP<PREC> & BelosPrec)
{
  precOk    = true;
  prec = BelosPrec;
}


//_________________________________________________________________________________________________
// INHERITED FUNCTIONS
//-------------------------------------------------------------------------------------------------
int
epetraOperatorPrecPrjA::
SetUseTranspose(bool UseTranspose)
{
  return(0);
}

int
epetraOperatorPrecPrjA::
Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  //Asserts
  assert(vectOk);
  assert(precOk);
  
  //Alloc
  Real vx, vv, vy;
  Epetra_MultiVector Z(X);
  
  vectorV->Norm2(&vv);
  vv = vv * vv;

  //Preconditioner
  prec->Apply(Z,Y);
  
  //Post-projection
  vectorV->Dot(Y,&vy);
  
  Y.Update(-vy / vv, *vectorV, 1.0);
  Y.Update(  b / vv, *vectorV, 1.0);
  
  return(0);
}

int
epetraOperatorPrecPrjA::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  return(1);
}

double
epetraOperatorPrecPrjA::
NormInf() const
{
  return(prec->NormInf());
}

const char *
epetraOperatorPrecPrjA::
Label() const
{
  return("operatorPrecPrjA");
}

bool
epetraOperatorPrecPrjA::
UseTranspose() const
{
  return(false);
}

bool
epetraOperatorPrecPrjA::
HasNormInf() const
{
  return(true);
}

const Epetra_Comm &
epetraOperatorPrecPrjA::
Comm() const
{
  assert(precOk);
  return(prec->Comm());
}

const Epetra_Map &
epetraOperatorPrecPrjA::
OperatorDomainMap() const
{
  return(prec->OperatorDomainMap());
}

const Epetra_Map &
epetraOperatorPrecPrjA::
OperatorRangeMap() const
{
  return(prec->OperatorRangeMap());
}
