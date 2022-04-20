/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef EPETRAOPERATORPRECDIAD_H
#define EPETRAOPERATORPRECDIAD_H

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


/*! Projected preconditioner operator */
class epetraOperatorPrecDiad : public Epetra_Operator
{
    /*! @name Typedefs */ //@{
  public:
    typedef Epetra_MultiVector VECTOR;
    typedef Epetra_Map         MAP;
    //@}
    
    /*! @name Data */ //@{
  public:
    bool vectOk;
    Real b;
    Teuchos::RCP<VECTOR> vectorV;
    Teuchos::RCP<VECTOR> vectorWt;
    Teuchos::RCP<MAP>    rangeMap, domainMap;
    //@}
    
    /*! @name Constructors and functions */ //@{
  public:
    epetraOperatorPrecDiad();
    epetraOperatorPrecDiad(const Real & b,
                           const Teuchos::RCP<VECTOR> & VectorV,
                           const Teuchos::RCP<VECTOR> & VectorWt,
                           const Teuchos::RCP<MAP>    & RangeMap,
                           const Teuchos::RCP<MAP>    & DomainMap);
    
    void setVector(const Real & b,
                   const Teuchos::RCP<VECTOR> & VectorV,
                   const Teuchos::RCP<VECTOR> & VectorWt,
                   const Teuchos::RCP<MAP>    & RangeMap,
                   const Teuchos::RCP<MAP>    & DomainMap);
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
