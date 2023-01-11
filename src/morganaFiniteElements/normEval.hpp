/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef NORMEVAL_HPP
#define NORMEVAL_HPP

#include "vectorBuilder.hpp"


template<typename LOCVECTOR>
class normEval : public vectorBuilder<LOCVECTOR>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename LOCVECTOR::FUNCTIONAL  FUNCTIONAL;
    typedef typename LOCVECTOR::PMAPTYPE    PMAPTYPE;
    //@}
    
    /*! @name Constructors and functions */ //@{
  public:
    normEval();
    normEval(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op);
    normEval(communicator & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op);
    Real eval();
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename LOCVECTOR>
normEval<LOCVECTOR>::
normEval() : vectorBuilder<LOCVECTOR>()
{
}

template<typename LOCVECTOR>
normEval<LOCVECTOR>::
normEval(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op) : vectorBuilder<LOCVECTOR>(CommDev,Op)
{
}

template<typename LOCVECTOR>
normEval<LOCVECTOR>::
normEval(communicator & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op) : vectorBuilder<LOCVECTOR>(CommDev,Op)
{
}

template<typename LOCVECTOR>
Real
normEval<LOCVECTOR>::
eval()
{
  typedef vectorBuilder<LOCVECTOR> VECTORBUILDER;
  
  //Evaluate vector
  Teuchos::RCP<Epetra_FEVector> vector;
  VECTORBUILDER::buildEpetraVector(vector);
  
  //Sum vector components
  double val;
  
  Epetra_MultiVector ones(vector->Map(),1);
  ones.PutScalar(1.0);
  ones.Dot(*vector,&val);
  
  return(val);
}

#endif
