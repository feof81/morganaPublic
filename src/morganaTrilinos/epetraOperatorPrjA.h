/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef EPETRAOPERATORPRJA_H
#define EPETRAOPERATORPRJA_H

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
class epetraOperatorPrjA : public Epetra_Operator
{  
    /*! @name Data */ //@{
  public:
    bool vectOk;
    Teuchos::RCP<Epetra_MultiVector> V;
    Teuchos::RCP<Epetra_Map> domainMap, rangeMap;
    //@}
    
    /*! @name Constructors and functions */ //@{
  public:
    epetraOperatorPrjA();
    epetraOperatorPrjA(const Epetra_MultiVector & VV,
                       const Epetra_Map & DomainMap,
                       const Epetra_Map & RangeMap);
    
    epetraOperatorPrjA(const Teuchos::RCP<Epetra_MultiVector> & VV,
                       const Teuchos::RCP<Epetra_Map> & DomainMap,
                       const Teuchos::RCP<Epetra_Map> & RangeMap);
    
    void setVector(const Epetra_MultiVector & VV,
                   const Epetra_Map & DomainMap,
                   const Epetra_Map & RangeMap);
    
    void setVector(const Teuchos::RCP<Epetra_MultiVector> & VV,
                   const Teuchos::RCP<Epetra_Map> & DomainMap,
                   const Teuchos::RCP<Epetra_Map> & RangeMap);
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
