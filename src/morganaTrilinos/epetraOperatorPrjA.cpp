/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "epetraOperatorPrjA.h"
#include "epetraMatrixManip.h"

//_________________________________________________________________________________________________
// CONSTRUCTOR AND FUNCTIONS
//-------------------------------------------------------------------------------------------------
epetraOperatorPrjA::
epetraOperatorPrjA()
{
  vectOk = false;
}

epetraOperatorPrjA::
epetraOperatorPrjA(const Epetra_MultiVector & VV,
                   const Epetra_Map & DomainMap,
                   const Epetra_Map & RangeMap)
{
  vectOk = true;
  
  V         = Teuchos::rcp(new Epetra_MultiVector(VV));
  domainMap = Teuchos::rcp(new Epetra_Map(DomainMap));
  rangeMap  = Teuchos::rcp(new Epetra_Map(RangeMap));
}

epetraOperatorPrjA::
epetraOperatorPrjA(const Teuchos::RCP<Epetra_MultiVector> & VV,
                   const Teuchos::RCP<Epetra_Map> & DomainMap,
                   const Teuchos::RCP<Epetra_Map> & RangeMap)
{
  vectOk = true;

  V         = VV;
  domainMap = DomainMap;
  rangeMap  = RangeMap;
}

void
epetraOperatorPrjA::
setVector(const Epetra_MultiVector & VV,
          const Epetra_Map & DomainMap,
          const Epetra_Map & RangeMap)
{
  vectOk = true;
  
  V         = Teuchos::rcp(new Epetra_MultiVector(VV));
  domainMap = Teuchos::rcp(new Epetra_Map(DomainMap));
  rangeMap  = Teuchos::rcp(new Epetra_Map(RangeMap));
}

void
epetraOperatorPrjA::
setVector(const Teuchos::RCP<Epetra_MultiVector> & VV,
          const Teuchos::RCP<Epetra_Map> & DomainMap,
          const Teuchos::RCP<Epetra_Map> & RangeMap)
{
  vectOk = true;
  
  V         = VV;
  domainMap = DomainMap;
  rangeMap  = RangeMap;
}


//_________________________________________________________________________________________________
// INHERITED FUNCTIONS
//-------------------------------------------------------------------------------------------------
int
epetraOperatorPrjA::
SetUseTranspose(bool UseTranspose)
{
  return(0);
}

int
epetraOperatorPrjA::
Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  assert(vectOk);

  Y = X;
  
  Real Vt_V, Vt_X;
  V->Dot(*V, &Vt_V);
  V->Dot( X, &Vt_X);

  Y.Update(- Vt_X / Vt_V, *V, 1.0);
  
  return(0);
}

int
epetraOperatorPrjA::
ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  return(1);
}

double
epetraOperatorPrjA::
NormInf() const
{
  Real val;
  V->NormInf(&val); 
  
  return(val);
}

const char *
epetraOperatorPrjA::
Label() const
{
  return("operatorPrjA");
}

bool
epetraOperatorPrjA::
UseTranspose() const
{
  return(false);
}

bool
epetraOperatorPrjA::
HasNormInf() const
{
  return(true);
}

const Epetra_Comm &
epetraOperatorPrjA::
Comm() const
{
  assert(vectOk);
  return(V->Comm());
}

const Epetra_Map &
epetraOperatorPrjA::
OperatorDomainMap() const
{
  assert(vectOk);
  return(*domainMap);
}

const Epetra_Map &
epetraOperatorPrjA::
OperatorRangeMap() const
{
  assert(vectOk);
  return(*rangeMap);
}
