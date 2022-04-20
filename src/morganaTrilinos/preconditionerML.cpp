/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "preconditionerML.h"


preconditionerML::
preconditionerML() : virtualPreconditioner<double, Epetra_MultiVector, Epetra_RowMatrix, Belos::EpetraPrecOp>()
{
}

bool
preconditionerML::
initialize()
{
  typedef virtualPreconditioner<double, Epetra_MultiVector, Epetra_RowMatrix, Belos::EpetraPrecOp> VIRTPRECON;
  assert(list.strong_count() != 0);
  
  Epetra_RowMatrix * matrixPtr = VIRTPRECON::A.get();
  prec = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*matrixPtr, *list, false));
  return(true);
}

bool
preconditionerML::
compute()
{
  prec->ComputePreconditioner();
  prec->AnalyzeCoarse();
  return(true);
}

Teuchos::RCP<Belos::EpetraPrecOp>
preconditionerML::
getPreconditioner()
{
  return(Teuchos::rcp(new Belos::EpetraPrecOp(prec)));
}
