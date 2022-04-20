/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef TPETRAVECTOR_TO_PVECT_H
#define TPETRAVECTOR_TO_PVECT_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map_decl.hpp"
#include "Tpetra_MultiVector_decl.hpp"

#include "pMap.hpp"
#include "pVect.hpp"
#include "pMapGlobalManip.h"
#include "pVectGlobalManip.hpp"


//_________________________________________________________________________________________________
// DUMMY EMPTY
//-------------------------------------------------------------------------------------------------
/*! Conversion from a \c Tpetra_MultiVector to \c pVect (DATA = Real).
 Empty specialization*/
template<typename MAP>
class tpetraVector_to_pVect
{
};


//_________________________________________________________________________________________________
// PMAPITEM SPECIALIZATION
//-------------------------------------------------------------------------------------------------
/*! Conversion from a \c Tpetra_MultiVector to \c pVect (DATA = Real).
 pMapItem specialization */
template<>
class tpetraVector_to_pVect<pMapItem>
{
    /*! @name Typedefs - Morgana */ //@{
  public:
    typedef pMapItem        MAP;
    typedef pMap<MAP>       PMAP;
    typedef pVect<Real,MAP> PVECT;
    //@}
    
    /*! @name Typedefs - Tpetra */ //@{
  public:
    typedef Teuchos::MpiComm<int> TPETRA_COMM;
    typedef Tpetra::Map<>         TPETRA_MAP;
    typedef Tpetra::MultiVector<> TPETRA_VECTOR;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructors and set */ //@{
  public:
    tpetraVector_to_pVect();
    tpetraVector_to_pVect(const Teuchos::RCP<communicator> & CommDev);
    tpetraVector_to_pVect(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Functions */ //@{
  public:
    PVECT convert(const Teuchos::RCP<TPETRA_VECTOR> & epVect) const;
    PVECT convert(const TPETRA_VECTOR & epVect) const;
    PVECT convert(const Teuchos::RCP<TPETRA_VECTOR> & epVect, const Teuchos::RCP<PMAP> & refMap) const;
    PVECT convert(const TPETRA_VECTOR & epVect, const PMAP & refMap) const;
    PVECT convertAndChange(const Teuchos::RCP<TPETRA_VECTOR> & epVect, const Teuchos::RCP<PMAP> & newMap) const;
    PVECT convertAndChange(const TPETRA_VECTOR & epVect, PMAP & newMap) const;
    //@}
};



//_________________________________________________________________________________________________
// PMAPITEMSHARE pMapItemShare
//-------------------------------------------------------------------------------------------------
/*! Conversion from a \c Tpetra_MultiVector to \c pVect (DATA = Real) */
template<>
class tpetraVector_to_pVect<pMapItemShare>
{
    /*! @name Typedefs - Morgana */ //@{
  public:
    typedef pMapItemShare   MAP;
    typedef pMap<MAP>       PMAP;
    typedef pVect<Real,MAP> PVECT;
    //@}
    
    /*! @name Typedefs - Tpetra */ //@{
  public:
    typedef Teuchos::MpiComm<int> TPETRA_COMM;
    typedef Tpetra::Map<>         TPETRA_MAP;
    typedef Tpetra::MultiVector<> TPETRA_VECTOR;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructors and set */ //@{
  public:
    tpetraVector_to_pVect();
    tpetraVector_to_pVect(const Teuchos::RCP<communicator> & CommDev);
    tpetraVector_to_pVect(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Functions */ //@{
  public:
    PVECT convert(const Teuchos::RCP<TPETRA_VECTOR> & epVect) const;
    PVECT convert(const TPETRA_VECTOR & epVect) const;
    PVECT convert(const Teuchos::RCP<TPETRA_VECTOR> & epVect, const Teuchos::RCP<PMAP> & refMap) const;
    PVECT convert(const TPETRA_VECTOR & epVect, const PMAP & refMap) const;
    PVECT convertAndChange(const Teuchos::RCP<TPETRA_VECTOR> & epVect, const Teuchos::RCP<PMAP> & newMap) const;
    PVECT convertAndChange(const TPETRA_VECTOR & epVect, PMAP & newMap) const;
    //@}
};


#endif
