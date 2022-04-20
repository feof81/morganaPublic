/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FACTORYPRECONDITIONER_HPP
#define FACTORYPRECONDITIONER_HPP

#include "virtualPreconditioner.hpp"
#include "preconditionerIfpack.h"
#include "preconditionerML.h"


enum morganaPreconditioner{ precIfpack=1, precMl=2 };


/*! Factory class for the preconditioners */
template<typename SCALARTYPE, typename MV, typename INOP, typename OUTOP>
class factoryPreconditioner
{
  /*! @name Typedefs */ //@{
  public:
    typedef virtualPreconditioner<SCALARTYPE,MV,INOP,OUTOP> PREC;
    //@}
  
    /*! @name Constructor and functions */ //@{
  public:
    factoryPreconditioner();
    Teuchos::RCP<PREC> create(const Teuchos::RCP<Teuchos::ParameterList> & List);
    Teuchos::RCP<PREC> create(const Teuchos::RCP<Teuchos::ParameterList> & List, const morganaPreconditioner & precTag);
    //@}
};


template<typename SCALARTYPE, typename MV, typename INOP, typename OUTOP>
factoryPreconditioner<SCALARTYPE,MV,INOP,OUTOP>::
factoryPreconditioner()
{
}


template<typename SCALARTYPE, typename MV, typename INOP, typename OUTOP>
Teuchos::RCP<typename factoryPreconditioner<SCALARTYPE,MV,INOP,OUTOP>::PREC>
factoryPreconditioner<SCALARTYPE,MV,INOP,OUTOP>::
create(const Teuchos::RCP<Teuchos::ParameterList> & List)
{
  assert(List.total_count() != 0);
  assert(List->isParameter(std::string("precClass")));
  
  morganaPreconditioner precTag = List->template get<morganaPreconditioner>( std::string("precClass") );
  Teuchos::RCP<PREC> pointer;
  
  switch(precTag)
  {
    case(precIfpack) :
    {
      pointer = Teuchos::rcp(new preconditionerIfpack());
      pointer->setParameters(List);
      return(pointer);
    }
      
    case(precMl) :
    {
      Teuchos::RCP<Teuchos::ParameterList> MLList(new Teuchos::ParameterList(List->sublist("mlList")));
      pointer = Teuchos::rcp(new preconditionerML());
      pointer->setParameters(MLList);
      return(pointer);
    }
      
    default:
      std::cout << "PRECONDITIONER NOT SPECIFIED!" << std::endl;
      assert(0);
  }
  
  return(pointer);
};

template<typename SCALARTYPE, typename MV, typename INOP, typename OUTOP>
Teuchos::RCP<typename factoryPreconditioner<SCALARTYPE,MV,INOP,OUTOP>::PREC>
factoryPreconditioner<SCALARTYPE,MV,INOP,OUTOP>::
create(const Teuchos::RCP<Teuchos::ParameterList> & List, const morganaPreconditioner & precTag)
{
  assert(List.strong_count() != 0);
  Teuchos::RCP<PREC> pointer;
  
  switch(precTag)
  {
    case(precIfpack) :
      pointer = Teuchos::rcp(new preconditionerIfpack());
      pointer->setParameters(List);
      return(pointer);
      
    case(precMl) :
      pointer = Teuchos::rcp(new preconditionerML());
      pointer->setParameters(List);
      return(pointer);
      
    default:
      std::cout << "PRECONDITIONER NOT SPECIFIED!" << std::endl;
      assert(0);
  }
  
  return(pointer);
}

#endif
