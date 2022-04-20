/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef TPETRASADDLEASSEMBLER_H
#define TPETRASADDLEASSEMBLER_H

#include "tpetraBlockMatrixAssemble.h"

/*! Saddle matrix assembler with diagonal contributions and extra-diagonal contributions */
class tpetraSaddleAssembler : public tpetraBlockMatrixAssemble
{
    /*! @name Typedefs */ //@{
  public:
    typedef Teuchos::MpiComm<int> COMM;
    typedef Tpetra::Map<>         MAP;
    typedef Tpetra::CrsMatrix<>   MATRIX;
    
    typedef Tpetra::global_size_t             TPETRA_GLOBAL_TYPE;
    typedef typename MAP::global_ordinal_type TPETRA_GLOBAL_ORDINAL;
    typedef typename MATRIX::scalar_type      TPETRA_SCALAR;
    //@}
    
    /*! @name Typedefs - RCP */ //@{
  public:
    typedef Teuchos::RCP<MAP>    RCP_MAP;
    typedef Teuchos::RCP<MATRIX> RCP_MATRIX;
    typedef Teuchos::RCP<COMM>   RCP_MPICOMM;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    tpetraSaddleAssembler();
    tpetraSaddleAssembler(const UInt & Dim);
    void setDim(const UInt & Dim);
    void setEpetraComm(const RCP_MPICOMM & CommDev);
    void setEpetraComm(const COMM        & CommDev);
    //@}
    
    /*! @name Topology */ //@{
  public:
    void setMap(const UInt & i, const MAP     & map);
    void setMap(const UInt & i, const RCP_MAP & map);
    void startup();
    //@}
    
    /*! @name Assemble */ //@{
  public:
    /*! Put \c sourceMatrix on the diagonal, position (i,i) */
    void setDiagMatrix(const MATRIX & sourceMatrix, const UInt & i);
    
    /*! Put \c sourceMatrix on the diagonal, position (i,i) */
    void setDiagMatrix(const RCP_MATRIX & sourceMatrix, const UInt & i);
    
    /*! Given a matrix M, put the matrix (M^T * M * scale) on the diagonal, position (i,i) */
    void setMTM(const MATRIX & sourceMatrix, const UInt & i, const Real & scale = 1.0);
    
    /*! Given a matrix M, put the matrix (M^T * M * scale) on the diagonal, position (i,i) */
    void setMTM(const RCP_MATRIX & sourceMatrix, const UInt & i, const Real & scale = 1.0);
    
    /*! Put the matrix (scale * I) on the diagonal, position (i,i) */
    void setEye(const UInt & i, const Real & scale = 1.0);
    
    /*! Put the \c sourceMatrix in position (i,j) and its transpose in position (j,i) */
    void setSymmMatrix(const MATRIX & sourceMatrix, const UInt & i, const UInt & j);
    
    /*! Put the \c sourceMatrix in position (i,j) and its transpose in position (j,i) */
    void setSymmMatrix(const RCP_MATRIX & sourceMatrix, const UInt & i, const UInt & j);
    
    /*! Global assemble */
    void assemble();
    //@}
    
    /*! @name Return */ //@{
  public:
    RCP_MAP    getMap() const;
    RCP_MATRIX getMatrix() const;
    //@}
};

#endif
