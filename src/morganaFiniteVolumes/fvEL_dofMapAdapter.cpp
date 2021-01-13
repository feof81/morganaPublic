/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/



#include "fvEL_dofMapAdapter.h"

fvEL_dofMapAdapter::
fvEL_dofMapAdapter()
{
  dataLoaded = false;
}

fvEL_dofMapAdapter::
fvEL_dofMapAdapter(const sVect<bool> & ItemIsActive, const sVect<UInt> & NewItemLid)
{
  assert(ItemIsActive.size() == NewItemLid.size());
  
  dataLoaded = true;
  itemIsActive = ItemIsActive;
  newItemLid   = NewItemLid;
  
  startup();
}

void
fvEL_dofMapAdapter::
setData(const sVect<bool> & ItemIsActive, const sVect<UInt> & NewItemLid)
{
  assert(ItemIsActive.size() == NewItemLid.size());
  
  dataLoaded = true;
  itemIsActive = ItemIsActive;
  newItemLid   = NewItemLid;
  
  startup();
}

void
fvEL_dofMapAdapter::
startup()
{
  //Count the new elements
  numNewItems = 0;
  
  for(UInt i=1; i <= itemIsActive.size(); ++i)
  { numNewItems += UInt(itemIsActive(i)); }
  
  
  //New to old map
  newToOld.resize(numNewItems);
  
  for(UInt i=1; i <= itemIsActive.size(); ++i)
  {
    if(itemIsActive(i))
    {
      newToOld( newItemLid(i) ) = i;
    }
  }
}

bool
fvEL_dofMapAdapter::
getOldIsActive(const UInt & i) const
{
  assert(dataLoaded == true);
  assert(i >= 1);
  assert(i <= itemIsActive.size());
  
  return(itemIsActive(i));
}

UInt
fvEL_dofMapAdapter::
getOldToNew(const UInt & i) const
{
  assert(dataLoaded == true);
  assert(i >= 1);
  assert(i <= itemIsActive.size());
  assert(itemIsActive(i));
  
  return(newItemLid(i));
}

UInt
fvEL_dofMapAdapter::
getNewToOld(const UInt & i) const
{
  assert(dataLoaded == true);
  assert(i >= 1);
  assert(i <= numNewItems);
  
  return(newToOld(i));
}

UInt
fvEL_dofMapAdapter::
getNumNewItems() const
{
  assert(dataLoaded == true);
  return(numNewItems);
}

