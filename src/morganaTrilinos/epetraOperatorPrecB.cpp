/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "epetraOperatorPrecB.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS AND FUNCTIONS
//-------------------------------------------------------------------------------------------------
epetraOperatorPrecB::
epetraOperatorPrecB()
{
  matrixOk = false;
}

epetraOperatorPrecB::
epetraOperatorPrecB(const Teuchos::RCP<Belos::EpetraPrecOp> & AA,
                    const Teuchos::RCP<Epetra_MultiVector>  & VV)
{
  //Control logic
  matrixOk = true;
  
  //Copy data
  belosPrec = AA;
  V         = VV;
  
  //Build vector invA_V
  invA_V = Teuchos::rcp(new Epetra_MultiVector(VV->Map(),1));
  belosPrec->Apply(*V,*invA_V);
}

epetraOperatorPrecB::
epetraOperatorPrecB(const epetraOperatorPrecB & op)
{
  //Control logic
  matrixOk = op.matrixOk;
  
  //Copy data
  belosPrec = op.belosPrec;
  V         = op.V;
  
  //Build vector invA_V
  invA_V = Teuchos::rcp(new Epetra_MultiVector(V->Map(),1));
  belosPrec->Apply(*V,*invA_V);
}    

epetraOperatorPrecB
epetraOperatorPrecB::
operator=(const epetraOperatorPrecB & op)
{
  //Control logic
  matrixOk = op.matrixOk;
  
  //Copy data
  belosPrec = op.belosPrec;
  V         = op.V;
  
  //Build vector invA_V
  invA_V = Teuchos::rcp(new Epetra_MultiVector(V->Map(),1));
  belosPrec->Apply(*V,*invA_V);
  
  //Return
  return(*this);
}

void
epetraOperatorPrecB::
setOperators(const Teuchos::RCP<Belos::EpetraPrecOp> & A,
             const Teuchos::RCP<Epetra_MultiVector>  & V)
{
  //Control logic
  matrixOk = true;
  
  //Copy data
  belosPrec = A;
  V         = V;
  
  //Build vector invA_V
  invA_V = Teuchos::rcp(new Epetra_MultiVector(V->Map(),1));
  belosPrec->Apply(*V,*invA_V);
}


//_________________________________________________________________________________________________
// INHERITED FUNCTIONS
//-------------------------------------------------------------------------------------------------
int
epetraOperatorPrecB::
SetUseTranspose(bool UseTranspose)
{
  assert(matrixOk);
  
  return(false);
}

int
epetraOperatorPrecB::
Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  //Control logic
  assert(matrixOk);
 
  //Compute
  Teuchos::RCP<Epetra_MultiVector> invA_X(new Epetra_MultiVector(X));
  belosPrec->Apply(X,*invA_X);
  
  Real vt_invA_x, vt_invA_v;
  invA_X->Dot(*V,&vt_invA_x);
  invA_V->Dot(*V,&vt_invA_v);
  
  Y = *invA_X;
  Y.Update( vt_invA_x / (1.0 - vt_invA_v), *invA_V, 1.0);
  
  return(0);
}

int
epetraOperatorPrecB::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  assert(matrixOk);
  assert(false);
  return(1);
}

double
epetraOperatorPrecB::
NormInf() const
{
  assert(matrixOk);
  
  return(belosPrec->NormInf());
}

const char *
epetraOperatorPrecB::
Label() const
{
  assert(matrixOk);
  
  return("operatorPrecB");
}

bool
epetraOperatorPrecB::
UseTranspose() const
{
  assert(matrixOk);
  
  return(false);
}

bool
epetraOperatorPrecB::
HasNormInf() const
{
  assert(matrixOk);
  
  return(belosPrec->HasNormInf());
}

const Epetra_Comm &
epetraOperatorPrecB::
Comm() const
{
  assert(matrixOk);
  
  return(belosPrec->Comm());
}

const Epetra_Map  &
epetraOperatorPrecB::
OperatorDomainMap() const
{
  assert(matrixOk);
  
  return(belosPrec->OperatorDomainMap());
}

const Epetra_Map  &
epetraOperatorPrecB::
OperatorRangeMap() const
{
  assert(matrixOk);
  
  return(belosPrec->OperatorRangeMap());
}
