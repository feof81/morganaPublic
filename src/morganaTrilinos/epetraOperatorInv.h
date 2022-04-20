/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef EPETRAOPERATORINV_H
#define EPETRAOPERATORINV_H

#include "typesInterface.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_RCP.hpp"

#include "Epetra_Map.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_RowMatrixTransposer.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"

#include "factoryPreconditioner.hpp"
#include "factoryLinearSolver.hpp"


/*! Projection operator */
class epetraOperatorInv : public Epetra_Operator
{
  /*! @name Typedefs */ //@{
  public:
    typedef factoryPreconditioner<double, Epetra_MultiVector, Epetra_RowMatrix, Belos::EpetraPrecOp> FACTORYPREC;
    typedef virtualPreconditioner<double, Epetra_MultiVector, Epetra_RowMatrix, Belos::EpetraPrecOp> PREC;
    
    typedef factoryLinearSolver<double, Epetra_MultiVector, Epetra_Operator> FACTORYSOLVER;
    typedef virtualLinearSolver<double, Epetra_MultiVector, Epetra_Operator> SOLVER;
    //@}
  
    /*! @name Data */ //@{
  public:
    bool  matrixOk;
    bool  solverOk;
    
    Teuchos::RCP<Epetra_CrsMatrix> A;
    
    Teuchos::ParameterList solverList;
    Teuchos::ParameterList preconList;
    
    Teuchos::RCP<PREC>                preconditioner;
    Teuchos::RCP<Belos::EpetraPrecOp> belosPrec;
    Teuchos::RCP<SOLVER>              linearSolver;
    //@}
    
    /*! @name Constructors and functions */ //@{
  public:
    epetraOperatorInv();
    epetraOperatorInv(const Teuchos::RCP<Epetra_CrsMatrix> & A);
    
    epetraOperatorInv(const Teuchos::RCP<Epetra_CrsMatrix> & A,
                      const Teuchos::ParameterList & solverList,
                      const Teuchos::ParameterList & preconditionerList);
    
    epetraOperatorInv(const epetraOperatorInv & op);
    
    epetraOperatorInv operator=(const epetraOperatorInv & op);
    
    void setOperator(const Teuchos::RCP<Epetra_CrsMatrix> & A);
    void setSolver(const Teuchos::ParameterList & solverList,
                   const Teuchos::ParameterList & preconditionerList);
    
    static Teuchos::RCP<epetraOperatorInv> op(const Teuchos::RCP<Epetra_CrsMatrix> & AA);
    //@}
    
    /*! @name Inhierited functions */ //@{
  public:
    int                 SetUseTranspose(bool UseTranspose);
    int                 Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const;
    int                 ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const;
    double              NormInf() const;
    const char        * Label() const;
    bool                UseTranspose() const;
    bool                HasNormInf() const;
    const Epetra_Comm & Comm() const;
    const Epetra_Map  & OperatorDomainMap() const;
    const Epetra_Map  & OperatorRangeMap() const;
    //@}
};

#endif
