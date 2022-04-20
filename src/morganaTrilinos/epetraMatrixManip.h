/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef EPETRAMATRIXMANIP_H
#define EPETRAMATRIXMANIP_H

#include "typesInterface.hpp"

#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_RCP.hpp"

#include "Epetra_Map.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_SerialComm.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_RowMatrixTransposer.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"

#include "EpetraExt_MatrixMatrix.h"


/*! Simle Matrix manipulation class */
class epetraMatrixManip
{
    /*! @name Typedefs */ //@{
  public:
    typedef Teuchos::RCP<Epetra_CrsMatrix>       RCP;
    typedef Teuchos::RCP<Epetra_Vector>          RCP_VECTOR;
    typedef Teuchos::RCP<const Epetra_CrsMatrix> RCP_CONST;
    //@}  
  
    /*! @name Constructors */ //@{
  public:
    epetraMatrixManip();
    //@}
    
    /*! @name Output functions */ //@{
  public:
    /*! Transpose */
    Epetra_CrsMatrix transpose(const Epetra_CrsMatrix & A) const;
    
    /*! Transpose */
    RCP transpose(const RCP_CONST & A) const;
    
    /*! diag(A) */
    Epetra_CrsMatrix diag(const Epetra_CrsMatrix & A) const;
    
    /*! diag(A) */
    RCP diag(const RCP_CONST & A) const;
    
    /*! diag(A)^{-1} */
    Epetra_CrsMatrix diagInv(const Epetra_CrsMatrix & A) const;
    
    /*! diag(A)^{-1} */
    RCP diagInv(const RCP_CONST & A) const;
    
    /*! Scaled eye matrix */
    Epetra_CrsMatrix getEyeRef(const Real       & scale,
                               const Epetra_Map & Map) const;
    
    /*! Scaled eye matrix */
    RCP getEyeRcp(const Real       & scale,
                  const Epetra_Map & Map) const;
    
    /*! Returns A * B */
    Epetra_CrsMatrix multiply(const Epetra_CrsMatrix & A,
                              const Epetra_CrsMatrix & B) const;
    
    /*! Returns A * B  */
    RCP multiply(const RCP_CONST & A,
                 const RCP_CONST & B) const;
    
    /*! Returns A * B * C  */
    Epetra_CrsMatrix multiply(const Epetra_CrsMatrix & A,
                              const Epetra_CrsMatrix & B,
                              const Epetra_CrsMatrix & C) const;
    
    /*! Returns A * B * C */
    RCP multiply(const RCP_CONST & A,
                 const RCP_CONST & B,
                 const RCP_CONST & C) const;
    
    
    /*! Returns B * diag(A)^{-1} * B^T */
    Epetra_CrsMatrix shurDiag(const Epetra_CrsMatrix & A,
                              const Epetra_CrsMatrix & B) const;
    
    /*! Returns B * diag(A)^{-1} * B^T */
    RCP shurDiag(const RCP_CONST & A,
                 const RCP_CONST & B) const;
    
    /*! Returns B * diag(A)^{-1} * C^T */
    Epetra_CrsMatrix shurDiag(const Epetra_CrsMatrix & A,
                              const Epetra_CrsMatrix & B,
                              const Epetra_CrsMatrix & C) const;
    
    /*! Returns A * diag(B)^{-1} * C^T */
    RCP shurDiag(const RCP_CONST & A,
                 const RCP_CONST & B,
                 const RCP_CONST & C) const;
    
    /*! Returns A * inv(B) * C^T */
    Epetra_CrsMatrix shur(const Epetra_CrsMatrix & A,
                          const Epetra_CrsMatrix & invB,
                          const Epetra_CrsMatrix & C) const;
    
    /*! Returns B * inv(A) * C^T */
    RCP shur(const RCP_CONST & A,
             const RCP_CONST & invB,
             const RCP_CONST & C) const;

    /*! Returns B * C^T */
    Epetra_CrsMatrix BCT(const Epetra_CrsMatrix & B,
                         const Epetra_CrsMatrix & C) const;
    
    /*! Returns B * C^T */
    RCP BCT(const RCP_CONST & B,
            const RCP_CONST & C) const;
    
    /*! The inverse of a 4x4 block matrix */
    Epetra_CrsMatrix inv4x4(const Epetra_CrsMatrix & A) const;
    
    /*! The inverse of a 4x4 block matrix */
    RCP inv4x4(const RCP_CONST & A) const;
    
    /*! MatrixSum */
    Epetra_CrsMatrix sum(const Epetra_CrsMatrix & A,
                         const Epetra_CrsMatrix & B);
    
    /*! MatrixSum */
    RCP sum(const RCP_CONST & A,
            const RCP_CONST & B) const;
    
    /*! MatrixSum */
    Epetra_CrsMatrix sum(const Epetra_CrsMatrix & A,
                         const Epetra_CrsMatrix & B,
                         const Epetra_CrsMatrix & C);
    
    /*! MatrixSum */
    RCP sum(const RCP_CONST & A,
            const RCP_CONST & B,
            const RCP_CONST & C) const;
            
    /*! Matrix LV */
    Epetra_CrsMatrix matrixLV(const Epetra_CrsMatrix & LV,
                              const Real & toll = 1.0e-10);
    
    /*! Matrix LV */
    RCP matrixLV(const RCP_CONST & LV,
                 const Real & toll = 1.0e-10) const;
                 
    Epetra_CrsMatrix taylorInv(const Epetra_CrsMatrix & A,
                               const UInt & it);
    
    RCP taylorInv(const RCP_CONST & A,
                  const UInt & it);
    //@}
};

#endif
