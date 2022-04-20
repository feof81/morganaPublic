/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef EPETRAVECTOR_TO_PVECT_H
#define EPETRAVECTOR_TO_PVECT_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_FEVector.h"

#include "pMap.hpp"
#include "pVect.hpp"
#include "pMapGlobalManip.h"
#include "pVectGlobalManip.hpp"


//_________________________________________________________________________________________________
// DUMMY EMPTY
//-------------------------------------------------------------------------------------------------
/*! Conversion from a \c Epetra_MultiVector to \c pVect (DATA = Real).
 Empty specialization*/
template<typename MAP>
class epetraVector_to_pVect
{
};


//_________________________________________________________________________________________________
// PMAPITEM SPECIALIZATION
//-------------------------------------------------------------------------------------------------
/*! Conversion from a \c Epetra_MultiVector to \c pVect (DATA = Real).
 pMapItem specialization */
template<>
class epetraVector_to_pVect<pMapItem>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItem        MAP;
    typedef pMap<MAP>       PMAP;
    typedef pVect<Real,MAP> PVECT;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructors and set */ //@{
  public:
    epetraVector_to_pVect();
    epetraVector_to_pVect(const Teuchos::RCP<communicator> & CommDev);
    epetraVector_to_pVect(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Functions */ //@{
  public:
    PVECT convert(const Teuchos::RCP<Epetra_MultiVector> & epVect) const;
    PVECT convert(const Epetra_MultiVector & epVect) const;
    PVECT convert(const Teuchos::RCP<Epetra_MultiVector> & epVect, const Teuchos::RCP<PMAP> & refMap) const;
    PVECT convert(const Epetra_MultiVector & epVect, const PMAP & refMap) const;
    PVECT convertAndChange(const Teuchos::RCP<Epetra_MultiVector> & epVect, const Teuchos::RCP<PMAP> & newMap) const;
    PVECT convertAndChange(const Epetra_MultiVector & epVect, PMAP & newMap) const;
    //@}
};



//_________________________________________________________________________________________________
// PMAPITEMSHARE pMapItemShare
//-------------------------------------------------------------------------------------------------
/*! Conversion from a \c Epetra_MultiVector to \c pVect (DATA = Real) */
template<>
class epetraVector_to_pVect<pMapItemShare>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItemShare   MAP;
    typedef pMap<MAP>       PMAP;
    typedef pVect<Real,MAP> PVECT;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructors and set */ //@{
  public:
    epetraVector_to_pVect();
    epetraVector_to_pVect(const Teuchos::RCP<communicator> & CommDev);
    epetraVector_to_pVect(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Functions */ //@{
  public:
    PVECT convert(const Teuchos::RCP<Epetra_MultiVector> & epVect) const;
    PVECT convert(const Epetra_MultiVector & epVect) const;
    PVECT convert(const Teuchos::RCP<Epetra_MultiVector> & epVect, const Teuchos::RCP<PMAP> & refMap) const;
    PVECT convert(const Epetra_MultiVector & epVect, const PMAP & refMap) const;
    PVECT convertAndChange(const Teuchos::RCP<Epetra_MultiVector> & epVect, const Teuchos::RCP<PMAP> & newMap) const;
    PVECT convertAndChange(const Epetra_MultiVector & epVect, PMAP & newMap) const;
    //@}
};


#endif
