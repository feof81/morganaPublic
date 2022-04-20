/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FACTORYLINEARSOLVER_HPP
#define FACTORYLINEARSOLVER_HPP

#include "virtualLinearSolver.hpp"
#include "belosBlockGmres.hpp"
#include "belosBlockCG.hpp"
#include "belosBlockBiCGStab.hpp"
#include "belosLSQR.hpp"


enum morganaLinSolver{ belosGmres=1, belosCG=2, belosBiCGStab=3, belosLsqr=4 };


/*! Factory class for the linear solvers */
template<typename SCALARTYPE, typename MV, typename OP>
class factoryLinearSolver
{
    /*! @name Typedefs */ //@{
  public:
    typedef virtualLinearSolver<SCALARTYPE,MV,OP> SOLVER;
    //@}
  
    /*! @name Constructor and functions */ //@{
  public:
    factoryLinearSolver();
    Teuchos::RCP<SOLVER> create(const Teuchos::RCP<Teuchos::ParameterList> & List);
    //@}
};


template<typename SCALARTYPE, typename MV, typename OP>
factoryLinearSolver<SCALARTYPE,MV,OP>::
factoryLinearSolver()
{
}

template<typename SCALARTYPE, typename MV, typename OP>
Teuchos::RCP<typename factoryLinearSolver<SCALARTYPE,MV,OP>::SOLVER>
factoryLinearSolver<SCALARTYPE,MV,OP>::
create(const Teuchos::RCP<Teuchos::ParameterList> & List)
{
  assert(List.strong_count() != 0);
  assert(List->isParameter(string("Solver")));
  
  morganaLinSolver solverTag = List->template get<morganaLinSolver>( string("Solver") );
  Teuchos::RCP<SOLVER> pointer;
  
  
  switch(solverTag)
  {
    case(belosGmres) :
      pointer = Teuchos::rcp(new belosBlockGmres<SCALARTYPE,MV,OP>());
      pointer->setParameters(List);
      return(pointer);
      
   case(belosCG) :
      pointer = Teuchos::rcp(new belosBlockCG<SCALARTYPE,MV,OP>());
      pointer->setParameters(List);
      return(pointer);
      
   case(belosBiCGStab) :
     pointer = Teuchos::rcp(new belosBlockBiCGStab<SCALARTYPE,MV,OP>());
      pointer->setParameters(List);
      return(pointer);
      
   case(belosLsqr) :
     pointer = Teuchos::rcp(new belosLSQR<SCALARTYPE,MV,OP>());
      pointer->setParameters(List);
      return(pointer);
      
    default:
      cout << "LINEAR SOLVER NOT SPECIFIED!" << endl;
      assert(0);
  }
  
  return(pointer);
}

#endif