/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "epetraVector_to_pVect.h"


//_________________________________________________________________________________________________
// PMAPITEM SPECIALIZATION
//-------------------------------------------------------------------------------------------------
epetraVector_to_pVect<pMapItem>::
epetraVector_to_pVect()
{
  commDevLoaded = false;
}

epetraVector_to_pVect<pMapItem>::
epetraVector_to_pVect(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

epetraVector_to_pVect<pMapItem>::
epetraVector_to_pVect(communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

void
epetraVector_to_pVect<pMapItem>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

void
epetraVector_to_pVect<pMapItem>::
setCommunicator(communicator & CommDev)
{
   commDevLoaded = true;
   commDev       = Teuchos::rcpFromRef(CommDev);
}

epetraVector_to_pVect<pMapItem>::PVECT
epetraVector_to_pVect<pMapItem>::
convert(const Teuchos::RCP<Epetra_MultiVector> & epVect) const
{
  assert(commDevLoaded);
  
  //Map extraction
  Epetra_BlockMap   epMap = epVect->Map();
  int         numElements = epMap.NumMyElements();
  int                base = epMap.IndexBase();
  int * newGlobalElements = epMap.MyGlobalElements();
  
  //Extract data
  const Epetra_MultiVector & vectRef = *epVect;
  double * vectPtr = vectRef[0];
  
  //pVect
  MAP item;
  item.setPid(commDev->rank());
  
  pVect<Real,MAP> mvect(numElements);
  
  for(UInt i=1; i <= UInt(numElements); ++i)
  {
    item.setLid(i);
    item.setGid( newGlobalElements[i-1] - base + 1 );
    
    mvect.getMapL(i)  = item;
    mvect.getDataL(i) = vectPtr[i-1];
  }
  
  mvect.updateFinder();
  
  //Check
  pMapGlobalManip<MAP> checker(commDev);
  assert(checker.check(mvect.getMapRef()));
  
  return(mvect);
}

epetraVector_to_pVect<pMapItem>::PVECT
epetraVector_to_pVect<pMapItem>::
convert(const Epetra_MultiVector & epVect) const
{
  assert(commDevLoaded);
  
  //Map extraction
  Epetra_BlockMap   epMap = epVect.Map();
  int         numElements = epMap.NumMyElements();
  int                base = epMap.IndexBase();
  int * newGlobalElements = epMap.MyGlobalElements();
  
  //Extract data
  const Epetra_MultiVector & vectRef = epVect;
  double * vectPtr = vectRef[0];
  
  //pVect
  MAP item;
  item.setPid(commDev->rank());
  
  pVect<Real,MAP> mvect(numElements);
  
  for(UInt i=1; i <= UInt(numElements); ++i)
  {
    item.setLid(i);
    item.setGid( newGlobalElements[i-1] - base + 1 );
    
    mvect.getMapL(i)  = item;
    mvect.getDataL(i) = vectPtr[i-1];
  }
  
  mvect.updateFinder();
  
  //Check
  pMapGlobalManip<MAP> checker(commDev);
  assert(checker.check(mvect.getMapRef()));
  
  return(mvect);
}

epetraVector_to_pVect<pMapItem>::PVECT
epetraVector_to_pVect<pMapItem>::
convert(const Teuchos::RCP<Epetra_MultiVector> & epVect, const Teuchos::RCP<PMAP> & refMap) const
{
  assert(commDevLoaded);
  
  //Map extraction
  Epetra_BlockMap   epMap = epVect->Map();
  int         numElements = epMap.NumMyElements();
  int                base = epMap.IndexBase();
  int * newGlobalElements = epMap.MyGlobalElements();
  
  assert(int(refMap->size()) == numElements);
  
  //Extract data
  const Epetra_MultiVector & vectRef = *epVect;
  double * vectPtr = vectRef[0];
  
  //pVect
  MAP item;
  
  pVect<Real,MAP> mvect(numElements);
  
  for(UInt i=1; i <= UInt(numElements); ++i)
  {
    item = refMap->get(i);
    
    assert(int(item.getPid()) == commDev->rank());
    assert(item.getLid() == i);
    assert(item.getGid() == (newGlobalElements[i-1] - base + 1) );
    
    mvect.getMapL(i)  = item;
    mvect.getDataL(i) = vectPtr[i-1];
  }
  
  mvect.updateFinder();
  
  //Check
  pMapGlobalManip<MAP> checker(commDev);
  assert(checker.check(mvect.getMapRef()));
  
  return(mvect);
}

epetraVector_to_pVect<pMapItem>::PVECT
epetraVector_to_pVect<pMapItem>::
convert(const Epetra_MultiVector & epVect, const PMAP & refMap) const
{
  assert(commDevLoaded);
  
  //Map extraction
  Epetra_BlockMap   epMap = epVect.Map();
  int         numElements = epMap.NumMyElements();
  int                base = epMap.IndexBase();
  int * newGlobalElements = epMap.MyGlobalElements();
  
  assert(refMap.size() == UInt(numElements));
  
  //Extract data
  const Epetra_MultiVector & vectRef = epVect;
  double * vectPtr = vectRef[0];
  
  //pVect
  MAP item;
  
  pVect<Real,MAP> mvect(numElements);
  
  for(UInt i=1; i <= UInt(numElements); ++i)
  {
    item = refMap(i);
    
    assert(item.getPid() == UInt(commDev->rank()));
    assert(item.getLid() == i);
    assert(item.getGid() == UInt(newGlobalElements[i-1] - base + 1) );
    
    mvect.getMapL(i)  = item;
    mvect.getDataL(i) = vectPtr[i-1];
  }
  
  mvect.updateFinder();
  
  //Check
  pMapGlobalManip<MAP> checker(commDev);
  assert(checker.check(mvect.getMapRef()));
  
  return(mvect);
}

epetraVector_to_pVect<pMapItem>::PVECT
epetraVector_to_pVect<pMapItem>::
convertAndChange(const Teuchos::RCP<Epetra_MultiVector> & epVect, const Teuchos::RCP<PMAP> & newMap) const
{
  PVECT tempVect = convert(epVect);
  
  pVectGlobalManip<Real,MAP> changer(commDev);
  changer.changeMap(tempVect, *newMap);
  
  assert(changer.check(tempVect));
  
  return(tempVect);
}

typename epetraVector_to_pVect<pMapItem>::PVECT
epetraVector_to_pVect<pMapItem>::
convertAndChange(const Epetra_MultiVector & epVect, PMAP & newMap) const
{
  PVECT tempVect = convert(epVect);
  
  pVectGlobalManip<Real,MAP> changer(commDev);
  changer.changeMap(tempVect, newMap);
  
  assert(changer.check(tempVect));
  
  return(tempVect);
}



//_________________________________________________________________________________________________
// PMAPITEMSHARE pMapItemShare
//-------------------------------------------------------------------------------------------------
epetraVector_to_pVect<pMapItemShare>::
epetraVector_to_pVect()
{
  commDevLoaded = false;
}

epetraVector_to_pVect<pMapItemShare>::
epetraVector_to_pVect(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

epetraVector_to_pVect<pMapItemShare>::
epetraVector_to_pVect(communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

void
epetraVector_to_pVect<pMapItemShare>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

void
epetraVector_to_pVect<pMapItemShare>::
setCommunicator(communicator & CommDev)
{
   commDevLoaded = true;
   commDev       = Teuchos::rcpFromRef(CommDev);
}

epetraVector_to_pVect<pMapItemShare>::PVECT
epetraVector_to_pVect<pMapItemShare>::
convert(const Teuchos::RCP<Epetra_MultiVector> & epVect) const
{
  assert(commDevLoaded);
  
  //Map extraction
  Epetra_BlockMap   epMap = epVect->Map();
  int         numElements = epMap.NumMyElements();
  int                base = epMap.IndexBase();
  int * newGlobalElements = epMap.MyGlobalElements();
  
  //Extract data
  const Epetra_MultiVector & vectRef = *epVect;
  double * vectPtr = vectRef[0];
  
  //pVect
  MAP item;
  item.setPid(commDev->rank());
  item.setShared(false);
  item.setOwned(true);
  
  pVect<Real,MAP> mvect(numElements);
  
  for(UInt i=1; i <= UInt(numElements); ++i)
  {
    item.setLid(i);
    item.setGid( newGlobalElements[i-1] - base + 1 );
    
    mvect.getMapL(i)  = item;
    mvect.getDataL(i) = vectPtr[i-1];
  }
  
  mvect.updateFinder();
  
  //Check
  pMapGlobalManip<MAP> checker(commDev);
  assert(checker.check(mvect.getMapRef()));
  
  return(mvect);
}

epetraVector_to_pVect<pMapItemShare>::PVECT
epetraVector_to_pVect<pMapItemShare>::
convert(const Epetra_MultiVector & epVect) const
{
  assert(commDevLoaded);
  
  //Map extraction
  Epetra_BlockMap   epMap = epVect.Map();
  int         numElements = epMap.NumMyElements();
  int                base = epMap.IndexBase();
  int * newGlobalElements = epMap.MyGlobalElements();
  
  //Extract data
  const Epetra_MultiVector & vectRef = epVect;
  double * vectPtr = vectRef[0];
  
  //pVect
  MAP item;
  item.setPid(commDev->rank());
  item.setShared(false);
  item.setOwned(true);
  
  pVect<Real,MAP> mvect(numElements);
  
  for(UInt i=1; i <= UInt(numElements); ++i)
  {
    item.setLid(i);
    item.setGid( newGlobalElements[i-1] - base + 1 );
    
    mvect.getMapL(i)  = item;
    mvect.getDataL(i) = vectPtr[i-1];
  }
  
  mvect.updateFinder();
  
  //Check
  pMapGlobalManip<MAP> checker(commDev);
  assert(checker.check(mvect.getMapRef()));
  
  return(mvect);
}

epetraVector_to_pVect<pMapItemShare>::PVECT
epetraVector_to_pVect<pMapItemShare>::
convert(const Teuchos::RCP<Epetra_MultiVector> & epVect, const Teuchos::RCP<PMAP> & refMap) const
{
  assert(commDevLoaded);
  
  //Map extraction
  Epetra_BlockMap   epMap = epVect->Map();
  int         numElements = epMap.NumMyElements();
  int                base = epMap.IndexBase();
  int * newGlobalElements = epMap.MyGlobalElements();
  
  assert(int(refMap->size()) == numElements);
  
  //Extract data
  const Epetra_MultiVector & vectRef = *epVect;
  double * vectPtr = vectRef[0];
  
  //pVect
  MAP item;
  
  pVect<Real,MAP> mvect(numElements);
  
  for(UInt i=1; i <= UInt(numElements); ++i)
  {
    item = refMap->get(i);
    
    assert(int(item.getPid()) == commDev->rank());
    assert(item.getLid() == i);
    assert(item.getGid() == (newGlobalElements[i-1] - base + 1) );
    
    mvect.getMapL(i)  = item;
    mvect.getDataL(i) = vectPtr[i-1];
  }
  
  mvect.updateFinder();
  
  //Check
  pMapGlobalManip<MAP> checker(commDev);
  assert(checker.check(mvect.getMapRef()));
  
  return(mvect);
}

epetraVector_to_pVect<pMapItemShare>::PVECT
epetraVector_to_pVect<pMapItemShare>::
convert(const Epetra_MultiVector & epVect, const PMAP & refMap) const
{
  assert(commDevLoaded);
  
  //Map extraction
  Epetra_BlockMap   epMap = epVect.Map();
  int         numElements = epMap.NumMyElements();
  int                base = epMap.IndexBase();
  int * newGlobalElements = epMap.MyGlobalElements();
  
  assert(refMap.size() == UInt(numElements));
  
  //Extract data
  const Epetra_MultiVector & vectRef = epVect;
  double * vectPtr = vectRef[0];
  
  //pVect
  MAP item;
  
  pVect<Real,MAP> mvect(numElements);
  
  for(UInt i=1; i <= UInt(numElements); ++i)
  {
    item = refMap(i);
    
    assert(item.getPid() == UInt(commDev->rank()));
    assert(item.getLid() == i);
    assert(item.getGid() == UInt(newGlobalElements[i-1] - base + 1) );
    
    mvect.getMapL(i)  = item;
    mvect.getDataL(i) = vectPtr[i-1];
  }
  
  mvect.updateFinder();
  
  //Check
  pMapGlobalManip<MAP> checker(commDev);
  assert(checker.check(mvect.getMapRef()));
  
  return(mvect);
}

epetraVector_to_pVect<pMapItemShare>::PVECT
epetraVector_to_pVect<pMapItemShare>::
convertAndChange(const Teuchos::RCP<Epetra_MultiVector> & epVect, const Teuchos::RCP<PMAP> & newMap) const
{
  PVECT tempVect = convert(epVect);
  
  pVectGlobalManip<Real,MAP> changer(commDev);
  changer.changeMap(tempVect, *newMap);
  
  assert(changer.check(tempVect));
  
  return(tempVect);
}

epetraVector_to_pVect<pMapItemShare>::PVECT
epetraVector_to_pVect<pMapItemShare>::
convertAndChange(const Epetra_MultiVector & epVect, PMAP & newMap) const
{
  PVECT tempVect = convert(epVect);
  
  pVectGlobalManip<Real,MAP> changer(commDev);
  changer.changeMap(tempVect, newMap);
  
  assert(changer.check(tempVect));
  
  return(tempVect);
}

