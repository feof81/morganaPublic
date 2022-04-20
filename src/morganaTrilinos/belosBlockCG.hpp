/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef BELOSBLOCKCG_HPP
#define BELOSBLOCKCG_HPP

#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosEpetraAdapter.hpp"

#include "virtualLinearSolver.hpp"


/*! Belos gmRes solver */
template<typename SCALARTYPE, typename MV, typename OP>
class belosBlockCG : public virtualLinearSolver<SCALARTYPE,MV,OP>
{
    /*! @name Constructor and functions */ //@{
  public:
    belosBlockCG();
    bool solve();
    //@}
};


template<typename SCALARTYPE, typename MV, typename OP>
belosBlockCG<SCALARTYPE,MV,OP>::
belosBlockCG() : virtualLinearSolver<SCALARTYPE,MV,OP>()
{
}

template<typename SCALARTYPE, typename MV, typename OP>
bool
belosBlockCG<SCALARTYPE,MV,OP>::
solve()
{
  typedef virtualLinearSolver<SCALARTYPE,MV,OP>     VIRTUALSOLVER;
  typedef Belos::LinearProblem<SCALARTYPE,MV,OP>    BELOS_LINEARPROBLEM;
  typedef Belos::BlockCGSolMgr<SCALARTYPE,MV,OP> BELOS_SOLVER;
  
  //build the linear problem
  Teuchos::RCP<BELOS_LINEARPROBLEM> problem = Teuchos::rcp(new BELOS_LINEARPROBLEM(VIRTUALSOLVER::A, VIRTUALSOLVER::X, VIRTUALSOLVER::T));
  
  //link the preconditioner
  cout << VIRTUALSOLVER::LP.strong_count() << endl;
  
  if(VIRTUALSOLVER::LP.strong_count() != 0)
  { problem->setLeftPrec(VIRTUALSOLVER::LP); }
  
  if(VIRTUALSOLVER::RP.strong_count() != 0)
  { problem->setRightPrec(VIRTUALSOLVER::RP); }
  
  //Consistency check
  bool flagSetProblem = problem->setProblem();
  assert(flagSetProblem);
  
  //Solve
  BELOS_SOLVER belosSolver(problem, VIRTUALSOLVER::list);
  
  try
  {
    Belos::ReturnType ret = belosSolver.solve();
    //return(ret == Belos::Converged);
    return(true);
  }
  catch(std::exception err)
  {
    cout << "Exception Belos" << endl;
    return(false);
  }
}

#endif
