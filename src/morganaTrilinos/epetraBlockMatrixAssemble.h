/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef EPETRABLOCKMATRIXASSEMBLE_H
#define EPETRABLOCKMATRIXASSEMBLE_H

#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_RCP.hpp"

#include "Epetra_ConfigDefs.h"
#include "Epetra_MpiComm.h"

#include "Epetra_DataAccess.h"
#include "Epetra_Map.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_RowMatrixTransposer.h"

#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

#include "simpleFormats.hpp"
#include "typesInterface.hpp"


/*! Block matrix assembler */
class epetraBlockMatrixAssemble
{
    /*! @name Typedefs */ //@{
  public:
    typedef Teuchos::RCP<Epetra_Map>         RCP_MAP;
    typedef Teuchos::RCP<Epetra_BlockMap>    RCP_BLOCKMAP;
    typedef Teuchos::RCP<Epetra_CrsMatrix>   RCP_MATRIX;
    typedef Teuchos::RCP<Epetra_FECrsMatrix> RCP_FEMATRIX;
    typedef Teuchos::RCP<Epetra_MpiComm>     RCP_MPICOMM;
    //@}
  
    /*! @name Logics */ //@{
  public:
    bool dimLoaded;
    bool commDevLoaded;
    sVect<bool> rowLoaded;
    sVect<bool> colLoaded;
    bool startupOk;
    bool assembleOK;
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
    RCP_MAP outRowMap, outColMap;
    RCP_FEMATRIX outMatrix;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    epetraBlockMatrixAssemble();
    epetraBlockMatrixAssemble(const UInt & NumRows, const UInt & NumCols);
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
    //@}
    
    /*! @name Assemble */ //@{
  public:
    void setMatrix(const Epetra_CrsMatrix & sourceMatrix, const UInt & i, const UInt & j);
    void setMatrix(const RCP_MATRIX & sourceMatrix, const UInt & i, const UInt & j);
    void setMatrixTransp(const Epetra_CrsMatrix & sourceMatrix, const UInt & i, const UInt & j);
    void setMatrixTransp(const RCP_MATRIX & sourceMatrix, const UInt & i, const UInt & j);
    void assemble();
    RCP_MAP getRowMap() const;
    RCP_MAP getColMap() const;
    RCP_FEMATRIX getMatrix() const;
    void printToFile(const string & name) const;
    //@}
};


#endif
