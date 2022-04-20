/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef PRECONDITIONERML_H
#define PRECONDITIONERML_H

#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"

#include "virtualPreconditioner.hpp"


/*! Ifpack preconditioner */
class preconditionerML : public virtualPreconditioner<double, Epetra_MultiVector, Epetra_RowMatrix, Belos::EpetraPrecOp>
{
    /*! @name Internal data */ //@{
  public:
    Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> prec;
    //@}
  
    /*! @name Constructor and functions */ //@{
  public:
    preconditionerML();
    bool initialize();
    bool compute();
    Teuchos::RCP<Belos::EpetraPrecOp> getPreconditioner();
    //@}
};

#endif
