/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include <string>
#include "preconditionerIfpack.h"

using namespace std;


preconditionerIfpack::
preconditionerIfpack() : virtualPreconditioner<double, Epetra_MultiVector, Epetra_RowMatrix, Belos::EpetraPrecOp>()
{
}

bool
preconditionerIfpack::
initialize()
{
  typedef virtualPreconditioner<double, Epetra_MultiVector, Epetra_RowMatrix, Belos::EpetraPrecOp> VIRTPRECON;
  
  assert(list.strong_count() != 0);
  assert(list->isParameter(string("precType")));
  assert(list->isParameter(string("overlap")));
  
  string precType = VIRTPRECON::list->get<string>( string("precType") );
  int    overlap  = VIRTPRECON::list->get<int>( string("overlap") );  
  
  Ifpack ifpackFactory;
  Epetra_RowMatrix * matrixPtr = VIRTPRECON::A.get();
  prec = Teuchos::rcp(ifpackFactory.Create(precType, matrixPtr, overlap));
  assert(prec != Teuchos::null);
  
  bool flagParams = prec->SetParameters(*VIRTPRECON::list);
  bool flagInit   = prec->Initialize();
  
  return(flagInit && flagParams);
}

bool
preconditionerIfpack::
compute()
{
  bool flagCompute = prec->Compute();  
  return(flagCompute);
}

Teuchos::RCP<Belos::EpetraPrecOp>
preconditionerIfpack::
getPreconditioner()
{
  return(Teuchos::rcp(new Belos::EpetraPrecOp(prec)));
}
