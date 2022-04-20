/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef EPETRAOPERATORBLOCK_H
#define EPETRAOPERATORBLOCK_H

#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_RCP.hpp"

#include "Epetra_ConfigDefs.h"
#include "Epetra_MpiComm.h"

#include "Epetra_DataAccess.h"
#include "Epetra_Map.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_RowMatrixTransposer.h"
#include "Epetra_Operator.h"

#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

#include "simpleFormats.hpp"
#include "typesInterface.hpp"

#include "epetraBlockVectorAssemble.h"
#include "epetraVectorStrip.h"


/*! Block matrix assembler */
class epetraOperatorBlock : public Epetra_Operator
{
    /*! @name Typedefs */ //@{
  public:
    typedef Teuchos::RCP<Epetra_Map>         RCP_MAP;
    typedef Teuchos::RCP<Epetra_BlockMap>    RCP_BLOCKMAP;
    typedef Teuchos::RCP<Epetra_MpiComm>     RCP_MPICOMM;
    typedef Teuchos::RCP<Epetra_MultiVector> RCP_MULTIVECT;
    typedef Teuchos::RCP<Epetra_Operator>    RCP_OPERATOR;
    //@}
    
    /*! @name Logics */ //@{
  public:
    bool dimLoaded;
    bool commDevLoaded;
    sVect<bool> rowLoaded;
    sVect<bool> colLoaded;
    bool startupOk;
    //@}
    
    /*! @name Links */ //@{
  public:
    RCP_MPICOMM commDev;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    UInt numRows, numCols;
    sVect<RCP_BLOCKMAP> rowMaps, colMaps;
    sVect<UInt> offsetRow, offsetCol;
    sVect<sVect<RCP_OPERATOR> > operatorMatrix;
    sVect<sVect<bool> >         boolMatrix;
    sVect<RCP_MULTIVECT> subRowVect, subColVect;
    RCP_MAP outRowMap, outColMap;
    //@}
    
    /*! @name Internal structures */ //@{
  public:
    epetraBlockVectorAssemble blockVect;
    epetraVectorStrip         vectStrip;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    epetraOperatorBlock();
    epetraOperatorBlock(const UInt & NumRows, const UInt & NumCols);
    void setDim(const UInt & NumRows, const UInt & NumCols);
    void setEpetraComm(const RCP_MPICOMM & CommDev);
    void setEpetraComm(const Epetra_MpiComm & CommDev);
    //@}
    
    /*! @name Topology */ //@{
  public:
    void setRowMap(const UInt & i, const Epetra_BlockMap & map);
    void setRowMap(const UInt & i, const RCP_BLOCKMAP & map);
    void setColMap(const UInt & i, const Epetra_BlockMap & map);
    void setColMap(const UInt & i, const RCP_BLOCKMAP & map);
    void startup();
    void setOperator(const RCP_OPERATOR & sourceOperator, const UInt & i, const UInt & j);
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
