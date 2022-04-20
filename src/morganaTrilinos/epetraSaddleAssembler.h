/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef EPETRASADDLEASSEMBLER_H
#define EPETRASADDLEASSEMBLER_H

#include "epetraBlockMatrixAssemble.h"

/*! Saddle matrix assembler with diagonal contributions and extra-diagonal contributions */
class epetraSaddleAssembler : public epetraBlockMatrixAssemble
{
    /*! @name Typedefs */ //@{
  public:
    typedef Teuchos::RCP<Epetra_Map>         RCP_MAP;
    typedef Teuchos::RCP<Epetra_BlockMap>    RCP_BLOCKMAP;
    typedef Teuchos::RCP<Epetra_CrsMatrix>   RCP_MATRIX;
    typedef Teuchos::RCP<Epetra_FECrsMatrix> RCP_FEMATRIX;
    typedef Teuchos::RCP<Epetra_MpiComm>     RCP_MPICOMM;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    epetraSaddleAssembler();
    epetraSaddleAssembler(const UInt & Dim);
    void setDim(const UInt & Dim);
    void setEpetraComm(const RCP_MPICOMM & CommDev);
    void setEpetraComm(const Epetra_MpiComm & CommDev);
    //@}
    
    /*! @name Topology */ //@{
  public:
    void setMap(const UInt & i, const Epetra_BlockMap & map);
    void setMap(const UInt & i, const RCP_BLOCKMAP & map);
    void startup();
    //@}
    
    /*! @name Assemble */ //@{
  public:
    /*! Put \c sourceMatrix on the diagonal, position (i,i) */
    void setDiagMatrix(const Epetra_CrsMatrix & sourceMatrix, const UInt & i);
    
    /*! Put \c sourceMatrix on the diagonal, position (i,i) */
    void setDiagMatrix(const RCP_MATRIX & sourceMatrix, const UInt & i);
    
    /*! Given a matrix M, put the matrix (M^T * M * scale) on the diagonal, position (i,i) */
    void setMTM(const Epetra_CrsMatrix & sourceMatrix, const UInt & i, const Real & scale = 1.0);
    
    /*! Given a matrix M, put the matrix (M^T * M * scale) on the diagonal, position (i,i) */
    void setMTM(const RCP_MATRIX & sourceMatrix, const UInt & i, const Real & scale = 1.0);
    
    /*! Given a matrix M, put the matrix (M * M^T * scale) on the diagonal, position (i,i) */
    void setMMT(const Epetra_CrsMatrix & sourceMatrix, const UInt & i, const Real & scale = 1.0);
    
    /*! Given a matrix M, put the matrix (M * M^T * scale) on the diagonal, position (i,i) */
    void setMMT(const RCP_MATRIX & sourceMatrix, const UInt & i, const Real & scale = 1.0);
    
    /*! Put the matrix (scale * I) on the diagonal, position (i,i) */
    void setEye(const UInt & i, const Real & scale = 1.0);
    
    /*! Put the \c sourceMatrix in position (i,j) and its transpose in position (j,i) */
    void setSymmMatrix(const Epetra_CrsMatrix & sourceMatrix, const UInt & i, const UInt & j);
    
    /*! Put the \c sourceMatrix in position (i,j) and its transpose in position (j,i) */
    void setSymmMatrix(const RCP_MATRIX & sourceMatrix, const UInt & i, const UInt & j);
    
    /*! Global assemble */
    void assemble();
    //@}
    
    /*! @name Return */ //@{
  public:
    RCP_MAP      getMap() const;
    RCP_FEMATRIX getMatrix() const;
    //@}
};

#endif
