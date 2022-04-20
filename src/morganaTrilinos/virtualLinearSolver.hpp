/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef VIRTUALLINEARSOLVER_HPP
#define VIRTUALLINEARSOLVER_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_RCPDecl.hpp"


/*! Virtual interface for linear solvers */
template<typename SCALARTYPE, typename MV, typename OP>
class virtualLinearSolver
{
    /*! @name Internal data */ //@{
  public:
    Teuchos::RCP<const OP> A;
    Teuchos::RCP<MV>       X;
    Teuchos::RCP<const MV> T;
    Teuchos::RCP<const OP> LP;
    Teuchos::RCP<const OP> RP;
    Teuchos::RCP<Teuchos::ParameterList> list;
    //@}
  
    /*! @name Constructor and functions */ //@{
  public:
    virtualLinearSolver();
    virtual ~virtualLinearSolver() {};
    void setOperator(const Teuchos::RCP<const OP> & AA);
    void setLHS(const Teuchos::RCP<MV> & XX) ;
    void setRHS(const Teuchos::RCP<const MV> & TT);
    void setProblem(const Teuchos::RCP<const OP> & AA, const Teuchos::RCP<MV> & XX, const Teuchos::RCP<const MV> & TT);
    void setLeftPrec(const Teuchos::RCP<const OP> & LPP);
    void setRightPrec(const Teuchos::RCP<const OP> & RPP);
    void setParameters(const Teuchos::RCP<Teuchos::ParameterList> & List);
    virtual bool solve() = 0;
    //@}
};


template<typename SCALARTYPE, typename MV, typename OP>
virtualLinearSolver<SCALARTYPE,MV,OP>::
virtualLinearSolver()
{
}

template<typename SCALARTYPE, typename MV, typename OP>
void
virtualLinearSolver<SCALARTYPE,MV,OP>::
setOperator(const Teuchos::RCP<const OP> & AA)
{
  A = AA;
}

template<typename SCALARTYPE, typename MV, typename OP>
void
virtualLinearSolver<SCALARTYPE,MV,OP>::
setLHS(const Teuchos::RCP<MV> & XX)
{
  X = XX;
}

template<typename SCALARTYPE, typename MV, typename OP>
void
virtualLinearSolver<SCALARTYPE,MV,OP>::
setRHS(const Teuchos::RCP<const MV> & TT)
{
  T = TT;
}

template<typename SCALARTYPE, typename MV, typename OP>
void
virtualLinearSolver<SCALARTYPE,MV,OP>::
setProblem(const Teuchos::RCP<const OP> & AA, const Teuchos::RCP<MV> & XX, const Teuchos::RCP<const MV> & TT)
{
  A = AA;
  X = XX;
  T = TT;
}

template<typename SCALARTYPE, typename MV, typename OP>
void
virtualLinearSolver<SCALARTYPE,MV,OP>::
setLeftPrec(const Teuchos::RCP<const OP> & LPP)
{
  LP = LPP;
}

template<typename SCALARTYPE, typename MV, typename OP>
void
virtualLinearSolver<SCALARTYPE,MV,OP>::
setRightPrec(const Teuchos::RCP<const OP> & RPP)
{
  RP = RPP;
}

template<typename SCALARTYPE, typename MV, typename OP>
void
virtualLinearSolver<SCALARTYPE,MV,OP>::
setParameters(const Teuchos::RCP<Teuchos::ParameterList> & List)
{
  list = List;
}

#endif
