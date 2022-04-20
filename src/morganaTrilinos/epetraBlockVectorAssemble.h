/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef EPETRABLOCKVECTORASSEMBLE_H
#define EPETRABLOCKVECTORASSEMBLE_H

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


class epetraBlockVectorAssemble
{
    /*! @name Typedefs */ //@{
  public:
    typedef Teuchos::RCP<Epetra_Map>         RCP_MAP;
    typedef Teuchos::RCP<Epetra_BlockMap>    RCP_BLOCKMAP;
    typedef Teuchos::RCP<Epetra_MultiVector> RCP_VECTOR;
    typedef Teuchos::RCP<Epetra_FEVector>    RCP_FEVECTOR;
    typedef Teuchos::RCP<Epetra_MpiComm>     RCP_MPICOMM;
    //@}
  
    /*! @name Logics */ //@{
  public:
    bool dimLoaded;
    bool commDevLoaded;
    sVect<bool> rowLoaded;
    bool startupOk;
    //@}
    
    /*! @name Links */ //@{
  public:
    RCP_MPICOMM commDev;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    UInt numRows;
    sVect<RCP_BLOCKMAP> rowMaps;
    sVect<UInt>  offsetRow;
    RCP_MAP      outRowMap;
    RCP_FEVECTOR outVector;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    epetraBlockVectorAssemble();
    epetraBlockVectorAssemble(const UInt & NumRows);
    void setDim(const UInt & NumRows);
    void setEpetraComm(const RCP_MPICOMM & CommDev);
    void setEpetraComm(const Epetra_MpiComm & CommDev);
    //@}
    
    /*! @name Topology */ //@{
  public:
    void setRowMap(const UInt & i, const Epetra_BlockMap & map);
    void setRowMap(const UInt & i, const RCP_BLOCKMAP & map);
    void startup();
    //@}
    
    /*! @name Assemble */ //@{
  public:
    void setVector(const Epetra_MultiVector & sourceVector, const UInt & i);
    void setVector(const RCP_VECTOR & sourceVector, const UInt & i);
    RCP_MAP getRowMap() const;
    RCP_FEVECTOR getVector() const;
    //@}
};


#endif
