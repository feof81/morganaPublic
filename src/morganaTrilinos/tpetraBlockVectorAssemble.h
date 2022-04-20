/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef TPETRABLOCKVECTORASSEMBLE_H
#define TPETRABLOCKVECTORASSEMBLE_H

#include "typesInterface.hpp"
#include "simpleFormats.hpp"

#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map_decl.hpp"
#include "Tpetra_MultiVector_decl.hpp"


/*! Block vector assembler */
class tpetraBlockVectorAssemble
{
    /*! @name Typedefs */ //@{
  public:
    typedef Teuchos::MpiComm<int> COMM;
    typedef Tpetra::Map<>         MAP;
    typedef Tpetra::MultiVector<> VECTOR;
    
    typedef Tpetra::global_size_t             TPETRA_GLOBAL_TYPE;
    typedef typename MAP::global_ordinal_type TPETRA_GLOBAL_ORDINAL;
    typedef typename VECTOR::scalar_type      TPETRA_SCALAR;
    //@}
    
    /*! @name Typedefs - RCP */ //@{
  public:
    typedef Teuchos::RCP<MAP>    RCP_MAP;
    typedef Teuchos::RCP<VECTOR> RCP_VECTOR;
    typedef Teuchos::RCP<COMM>   RCP_MPICOMM;
    //@}
    
    /*! @name Logics */ //@{
  public:
    bool dimLoaded;
    bool commDevLoaded;
    sVect<bool> rowLoaded;
    bool startupOk;
    bool assembleOK;
    //@}
    
    /*! @name Links */ //@{
  public:
    RCP_MPICOMM commDev;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    UInt numRows;
    sVect<RCP_MAP> rowMaps;
    sVect<UInt>  offsetRow;
    RCP_MAP      outRowMap;
    RCP_VECTOR   outVector;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    tpetraBlockVectorAssemble();
    tpetraBlockVectorAssemble(const UInt & NumRows);
    void setDim(const UInt & NumRows);
    void setTpetraComm(const RCP_MPICOMM & CommDev);
    void setTpetraComm(const COMM        & CommDev);
    //@}
    
    /*! @name Topology */ //@{
  public:
    void setRowMap(const UInt & i, const MAP     & map);
    void setRowMap(const UInt & i, const RCP_MAP & map);
    void startup();
    //@}
    
    /*! @name Assemble */ //@{
  public:
    void setVector(const VECTOR     & sourceVector, const UInt & i);
    void setVector(const RCP_VECTOR & sourceVector, const UInt & i);
    RCP_MAP    getRowMap() const;
    RCP_VECTOR getVector() const;
    //@}
};


#endif
