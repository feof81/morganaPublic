/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef VECTORBUILDER_HPP
#define VECTORBUILDER_HPP

#include "traitsVectorBuilder.hpp"


/*! Vector builder */
template<typename LOCVECTOR>
class vectorBuilder : public traitsVectorBuilder<LOCVECTOR, typename LOCVECTOR::PMAPTYPE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename LOCVECTOR::FUNCTIONAL  FUNCTIONAL;
    typedef typename LOCVECTOR::PMAPTYPE    PMAPTYPE;
    //@}
    
    /*! @name Constructors and set */ //@{
  public:
    vectorBuilder();
    vectorBuilder(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op);
    vectorBuilder(communicator & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename LOCVECTOR>
vectorBuilder<LOCVECTOR>::
vectorBuilder() : traitsVectorBuilder<LOCVECTOR, typename LOCVECTOR::PMAPTYPE>()
{
}

template<typename LOCVECTOR>
vectorBuilder<LOCVECTOR>::
vectorBuilder(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op) : traitsVectorBuilder<LOCVECTOR, typename LOCVECTOR::PMAPTYPE>(CommDev,Op)
{
}

template<typename LOCVECTOR>
vectorBuilder<LOCVECTOR>::
vectorBuilder(communicator & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op) : traitsVectorBuilder<LOCVECTOR, typename LOCVECTOR::PMAPTYPE>(CommDev,Op)
{
}


#endif
