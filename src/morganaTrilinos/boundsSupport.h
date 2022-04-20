/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef BOUNDSSUPPORT_H
#define BOUNDSSUPPORT_H

#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_RCP.hpp"

#include "Epetra_ConfigDefs.h"
#include "Epetra_MpiComm.h"

#include "Epetra_DataAccess.h"
#include "Epetra_Map.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_RowMatrixTransposer.h"

#include "typesInterface.hpp"

 

/*! Support for the implementation of bounds */
class boundsSupport
{
    /*! @name Typedefs */ //@{
  public:
    typedef Teuchos::RCP<Epetra_MpiComm>     RCP_MPICOMM;
    typedef Teuchos::RCP<Epetra_Map>         RCP_BLOCKMAP;
    typedef Teuchos::RCP<Epetra_MultiVector> RCP_VECTOR;
    typedef Teuchos::RCP<Epetra_FECrsMatrix>   RCP_MATRIX;
    //@}
    
    /*! @name Links and internal data */ //@{
  public:
    int pid;
    RCP_MPICOMM commDev;
    bool commDevLoaded;
    //@}
    
    /*! @name Constructors and set */ //@{
  public:
    boundsSupport();
    boundsSupport(const Epetra_MpiComm & CommDev);
    boundsSupport(const RCP_MPICOMM & CommDev);
    void setCommDev(const Epetra_MpiComm & CommDev);
    void setCommDev(const RCP_MPICOMM & CommDev);
    //@}
    
    /*! @name Functions */ //@{
  public:
    void increment();
    int  getPid();
    Epetra_Map          getMapRef();
    RCP_BLOCKMAP        getMapRcp();
    Epetra_MultiVector  getSourceRef(const Real & value);
    RCP_VECTOR          getSourceRcp(const Real & value);
    Epetra_FECrsMatrix  vector_to_matrix(const Epetra_MultiVector & vector);
    RCP_MATRIX          vector_to_matrix(const RCP_VECTOR & vector);
    //@}
};


#endif
