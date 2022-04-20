/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "epetraOperatorPrecDiad.h"


//_________________________________________________________________________________________________
// CONSTRUCTOR AND FUNCTIONS
//-------------------------------------------------------------------------------------------------
epetraOperatorPrecDiad::
epetraOperatorPrecDiad()
{
  vectOk = false;
  b      = 1.0;
}

epetraOperatorPrecDiad::
epetraOperatorPrecDiad(const Real & B,
                       const Teuchos::RCP<VECTOR> & VectorV,
                       const Teuchos::RCP<VECTOR> & VectorWt,
                       const Teuchos::RCP<MAP>    & RangeMap,
                       const Teuchos::RCP<MAP>    & DomainMap)
{
  vectOk = true;
  
  b = B;
  
  vectorV  = VectorV;
  vectorWt = VectorWt;
  
  rangeMap  = RangeMap;
  domainMap = DomainMap;
}

void
epetraOperatorPrecDiad::
setVector(const Real & B,
          const Teuchos::RCP<VECTOR> & VectorV,
          const Teuchos::RCP<VECTOR> & VectorWt,
          const Teuchos::RCP<MAP>    & RangeMap,
          const Teuchos::RCP<MAP>    & DomainMap)
{
  vectOk  = true;
  
  b = B;
  
  vectorV  = VectorV;
  vectorWt = VectorWt;
  
  rangeMap  = RangeMap;
  domainMap = DomainMap;
}


//_________________________________________________________________________________________________
// INHERITED FUNCTIONS
//-------------------------------------------------------------------------------------------------
int
epetraOperatorPrecDiad::
SetUseTranspose(bool UseTranspose)
{
  return(0);
}

int
epetraOperatorPrecDiad::
Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  //Asserts
  assert(vectOk);
  
  //Alloc
  Real vx;

  //Compute
  vectorWt->Dot(X,&vx);
  Y = *vectorV;
  Y.Scale(vx * b);
  
  return(0);
}

int
epetraOperatorPrecDiad::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  return(1);
}

double
epetraOperatorPrecDiad::
NormInf() const
{
  UInt numVect = vectorV->NumVectors();
  double norms[numVect];
  double maxVal = 0.0;
  
  vectorV->NormInf(&norms[0]);
  
  for(UInt i=0; i < numVect; ++i)
  { maxVal = std::max(maxVal,norms[i]); }

  return(maxVal);
}

const char *
epetraOperatorPrecDiad::
Label() const
{
  return("operatorDiad");
}

bool
epetraOperatorPrecDiad::
UseTranspose() const
{
  return(false);
}

bool
epetraOperatorPrecDiad::
HasNormInf() const
{
  return(true);
}

const Epetra_Comm &
epetraOperatorPrecDiad::
Comm() const
{
  return(vectorV->Comm());
}

const Epetra_Map &
epetraOperatorPrecDiad::
OperatorDomainMap() const
{
  return(*domainMap);
}

const Epetra_Map &
epetraOperatorPrecDiad::
OperatorRangeMap() const
{
  return(*rangeMap);
}
