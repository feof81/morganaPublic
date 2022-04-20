/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "tpetraVector_to_pVect.h"


//_________________________________________________________________________________________________
// PMAPITEM SPECIALIZATION
//-------------------------------------------------------------------------------------------------
tpetraVector_to_pVect<pMapItem>::
tpetraVector_to_pVect()
{
  commDevLoaded = false;
}

tpetraVector_to_pVect<pMapItem>::
tpetraVector_to_pVect(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

tpetraVector_to_pVect<pMapItem>::
tpetraVector_to_pVect(communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

void
tpetraVector_to_pVect<pMapItem>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

void
tpetraVector_to_pVect<pMapItem>::
setCommunicator(communicator & CommDev)
{
   commDevLoaded = true;
   commDev       = Teuchos::rcpFromRef(CommDev);
}

tpetraVector_to_pVect<pMapItem>::PVECT
tpetraVector_to_pVect<pMapItem>::
convert(const Teuchos::RCP<TPETRA_VECTOR> & epVect) const
{
  return(convert(*epVect));
}

tpetraVector_to_pVect<pMapItem>::PVECT
tpetraVector_to_pVect<pMapItem>::
convert(const TPETRA_VECTOR & epVect) const
{
  assert(commDevLoaded);
  
  //Map extraction
  Teuchos::RCP<TPETRA_MAP> tpMap = epVect.getMap();
  int numElements = tpMap->getNodeNumElements();
  int        base = tpMap->getIndexBase();
  
  //Extract data
  Teuchos::ArrayRCP<const double> vectPtr = epVect.getData(0);
  
  //pVect
  MAP item;
  item.setPid(commDev->rank());
  
  pVect<Real,MAP> mvect(numElements);
  
  for(UInt i=1; i <= UInt(numElements); ++i)
  {
    item.setLid(i);
    item.setGid( tpMap->getGlobalElement(i-1) - base + 1 );
    
    mvect.getMapL(i)  = item;
    mvect.getDataL(i) = vectPtr[i-1];
  }
  
  mvect.updateFinder();
  
  //Check
  pMapGlobalManip<MAP> checker(commDev);
  assert(checker.check(mvect.getMapRef()));
  
  return(mvect);
}

tpetraVector_to_pVect<pMapItem>::PVECT
tpetraVector_to_pVect<pMapItem>::
convert(const Teuchos::RCP<TPETRA_VECTOR> & epVect, const Teuchos::RCP<PMAP> & refMap) const
{
  return(convert(epVect,refMap));
}

tpetraVector_to_pVect<pMapItem>::PVECT
tpetraVector_to_pVect<pMapItem>::
convert(const TPETRA_VECTOR & epVect, const PMAP & refMap) const
{
  assert(commDevLoaded);
  
  //Map extraction
  Teuchos::RCP<TPETRA_MAP> tpMap = epVect.getMap();
  int numElements = tpMap->getNodeNumElements();
  int        base = tpMap->getIndexBase();
  
  assert(int(refMap.size()) == numElements);
  
  //Extract data
  Teuchos::ArrayRCP<const double> vectPtr = epVect.getData(0);
  
  //pVect
  MAP item;
  
  pVect<Real,MAP> mvect(numElements);
  
  for(UInt i=1; i <= UInt(numElements); ++i)
  {
    item = refMap.get(i);
    
    assert(int(item.getPid()) == commDev->rank());
    assert(item.getLid() == i);
    assert(item.getGid() == (tpMap->getGlobalElement(i-1) - base + 1) );
    
    mvect.getMapL(i)  = item;
    mvect.getDataL(i) = vectPtr[i-1];
  }
  
  mvect.updateFinder();
}

tpetraVector_to_pVect<pMapItem>::PVECT
tpetraVector_to_pVect<pMapItem>::
convertAndChange(const Teuchos::RCP<TPETRA_VECTOR> & epVect, const Teuchos::RCP<PMAP> & newMap) const
{
  PVECT tempVect = convert(epVect);
  
  pVectGlobalManip<Real,MAP> changer(commDev);
  changer.changeMap(tempVect, *newMap);
  
  assert(changer.check(tempVect));
  
  return(tempVect);
}

typename tpetraVector_to_pVect<pMapItem>::PVECT
tpetraVector_to_pVect<pMapItem>::
convertAndChange(const TPETRA_VECTOR & epVect, PMAP & newMap) const
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
tpetraVector_to_pVect<pMapItemShare>::
tpetraVector_to_pVect()
{
  commDevLoaded = false;
}

tpetraVector_to_pVect<pMapItemShare>::
tpetraVector_to_pVect(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

tpetraVector_to_pVect<pMapItemShare>::
tpetraVector_to_pVect(communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

void
tpetraVector_to_pVect<pMapItemShare>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

void
tpetraVector_to_pVect<pMapItemShare>::
setCommunicator(communicator & CommDev)
{
   commDevLoaded = true;
   commDev       = Teuchos::rcpFromRef(CommDev);
}

tpetraVector_to_pVect<pMapItemShare>::PVECT
tpetraVector_to_pVect<pMapItemShare>::
convert(const Teuchos::RCP<TPETRA_VECTOR> & epVect) const
{
  return(convert(*epVect));
}

tpetraVector_to_pVect<pMapItemShare>::PVECT
tpetraVector_to_pVect<pMapItemShare>::
convert(const TPETRA_VECTOR & epVect) const
{
  assert(commDevLoaded);
  
  //Map extraction
  Teuchos::RCP<TPETRA_MAP> tpMap = epVect.getMap();
  int numElements = tpMap->getNodeNumElements();
  int        base = tpMap->getIndexBase();
  
  //Extract data
  Teuchos::ArrayRCP<const double> vectPtr = epVect.getData(0);
  
  //pVect
  MAP item;
  item.setPid(commDev->rank());
  item.setShared(false);
  item.setOwned(true);
  
  pVect<Real,MAP> mvect(numElements);
  
  for(UInt i=1; i <= UInt(numElements); ++i)
  {
    item.setLid(i);
    item.setGid(tpMap->getGlobalElement(i-1) - base + 1 );
    
    mvect.getMapL(i)  = item;
    mvect.getDataL(i) = vectPtr[i-1];
  }
  
  mvect.updateFinder();
  
  //Check
  pMapGlobalManip<MAP> checker(commDev);
  assert(checker.check(mvect.getMapRef()));
  
  return(mvect);
}

tpetraVector_to_pVect<pMapItemShare>::PVECT
tpetraVector_to_pVect<pMapItemShare>::
convert(const Teuchos::RCP<TPETRA_VECTOR> & epVect, const Teuchos::RCP<PMAP> & refMap) const
{
  return(convert(*epVect,*refMap));
}

tpetraVector_to_pVect<pMapItemShare>::PVECT
tpetraVector_to_pVect<pMapItemShare>::
convert(const TPETRA_VECTOR & epVect, const PMAP & refMap) const
{
  assert(commDevLoaded);
  
  //Map extraction
  Teuchos::RCP<TPETRA_MAP> tpMap = epVect.getMap();
  int numElements = tpMap->getNodeNumElements();
  int        base = tpMap->getIndexBase();
  
  assert(refMap.size() == UInt(numElements));
  
  //Extract data
  Teuchos::ArrayRCP<const double> vectPtr = epVect.getData(0);
  
  //pVect
  MAP item;
  
  pVect<Real,MAP> mvect(numElements);
  
  for(UInt i=1; i <= UInt(numElements); ++i)
  {
    item = refMap(i);
    
    assert(item.getPid() == UInt(commDev->rank()));
    assert(item.getLid() == i);
    assert(item.getGid() == UInt(tpMap->getGlobalElement(i-1) - base + 1) );
    
    mvect.getMapL(i)  = item;
    mvect.getDataL(i) = vectPtr[i-1];
  }
  
  mvect.updateFinder();
  
  //Check
  pMapGlobalManip<MAP> checker(commDev);
  assert(checker.check(mvect.getMapRef()));
  
  return(mvect);
}

tpetraVector_to_pVect<pMapItemShare>::PVECT
tpetraVector_to_pVect<pMapItemShare>::
convertAndChange(const Teuchos::RCP<TPETRA_VECTOR> & epVect, const Teuchos::RCP<PMAP> & newMap) const
{
  PVECT tempVect = convert(epVect);
  
  pVectGlobalManip<Real,MAP> changer(commDev);
  changer.changeMap(tempVect, *newMap);
  
  assert(changer.check(tempVect));
  
  return(tempVect);
}

tpetraVector_to_pVect<pMapItemShare>::PVECT
tpetraVector_to_pVect<pMapItemShare>::
convertAndChange(const TPETRA_VECTOR & epVect, PMAP & newMap) const
{
  PVECT tempVect = convert(epVect);
  
  pVectGlobalManip<Real,MAP> changer(commDev);
  changer.changeMap(tempVect, newMap);
  
  assert(changer.check(tempVect));
  
  return(tempVect);
}
