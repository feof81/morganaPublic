/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef VIRTUALPRECONDITIONER_HPP
#define VIRTUALPRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"


/*! Virtual interface for preconditioners */
template<typename SCALARTYPE, typename MV, typename INOP, typename OUTOP>
class virtualPreconditioner
{
    /*! @name Internal data */ //@{
  public:
    Teuchos::RCP<const INOP> A;
    Teuchos::RCP<Teuchos::ParameterList> list;
    //@}
    
    /*! @name Constructor and functions */ //@{
  public:
    virtualPreconditioner();
    virtual ~virtualPreconditioner() {};
    void setOperator(const Teuchos::RCP<const INOP> & AA);
    void setParameters(const Teuchos::RCP<Teuchos::ParameterList> & List);
    virtual bool initialize() = 0;
    virtual bool compute() = 0;
    virtual Teuchos::RCP<OUTOP> getPreconditioner() = 0;
    //@}
};


template<typename SCALARTYPE, typename MV, typename INOP, typename OUTOP>
virtualPreconditioner<SCALARTYPE,MV,INOP,OUTOP>::
virtualPreconditioner()
{
}

template<typename SCALARTYPE, typename MV, typename INOP, typename OUTOP>
void
virtualPreconditioner<SCALARTYPE,MV,INOP,OUTOP>::
setOperator(const Teuchos::RCP<const INOP> & AA)
{
  A = AA;
}

template<typename SCALARTYPE, typename MV, typename INOP, typename OUTOP>
void
virtualPreconditioner<SCALARTYPE,MV,INOP,OUTOP>::
setParameters(const Teuchos::RCP<Teuchos::ParameterList> & List)
{
  list = List;
}


#endif
