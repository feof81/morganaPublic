/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef MATRIXBUILDER_HPP
#define MATRIXBUILDER_HPP

#include "traitsMatrixBuilder.hpp"


/*! Matrix builder */
template<typename LOCMATRIX>
class matrixBuilder : public traitsMatrixBuilder<LOCMATRIX, typename LOCMATRIX::PMAPTYPE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename LOCMATRIX::OPERATOR  OPERATOR;
    typedef typename LOCMATRIX::PMAPTYPE  PMAPTYPE;
    //@}
    
    /*! @name Constructors and set */ //@{
  public:
    matrixBuilder();
    matrixBuilder(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<OPERATOR> & Op);
    matrixBuilder(communicator & CommDev, const Teuchos::RCP<OPERATOR> & Op);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename LOCMATRIX>
matrixBuilder<LOCMATRIX>::
matrixBuilder() : traitsMatrixBuilder<LOCMATRIX, typename LOCMATRIX::PMAPTYPE>()
{
}

template<typename LOCMATRIX>
matrixBuilder<LOCMATRIX>::
matrixBuilder(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<OPERATOR> & Op) : traitsMatrixBuilder<LOCMATRIX, typename LOCMATRIX::PMAPTYPE>(CommDev,Op)
{
}


template<typename LOCMATRIX>
matrixBuilder<LOCMATRIX>::
matrixBuilder(communicator & CommDev, const Teuchos::RCP<OPERATOR> & Op) : traitsMatrixBuilder<LOCMATRIX, typename LOCMATRIX::PMAPTYPE>(CommDev,Op)
{
}

#endif
