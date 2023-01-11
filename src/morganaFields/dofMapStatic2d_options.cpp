/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "dofMapStatic2d_options.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
dofMapStatic2d_options::
dofMapStatic2d_options()
{
}

dofMapStatic2d_options::
dofMapStatic2d_options(const dofMapStatic2d_options & O)
{
  typedef typename set<UInt>::iterator ITERATOR;
  ITERATOR iter;
  
  activeIds.clear();
  
  for(iter = O.activeIds.begin(); iter != O.activeIds.end(); ++iter)
  { activeIds.insert(*iter); }
  
  blockIds = O.blockIds;
}

dofMapStatic2d_options
dofMapStatic2d_options::
operator=(const dofMapStatic2d_options & O)
{
  typedef typename set<UInt>::iterator ITERATOR;
  ITERATOR iter;
  
  activeIds.clear();
  
  for(iter = O.activeIds.begin(); iter != O.activeIds.end(); ++iter)
  { activeIds.insert(*iter); }
  
  blockIds = O.blockIds;
  
  return(*this);
}

void
dofMapStatic2d_options::
clear()
{
  activeIds.clear();
}


//_________________________________________________________________________________________________
// GEOID FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
dofMapStatic2d_options::
addGeoId(const UInt & Id)
{
  activeIds.insert(Id);
}

UInt
dofMapStatic2d_options::
numGeoIds() const
{
  return(activeIds.size());
}

bool
dofMapStatic2d_options::
isGeoId(const UInt & Id) const
{
  return(activeIds.count(Id) > 0);
}

sVect<UInt>
dofMapStatic2d_options::
geoIdsList() const
{
  sVect<UInt> outList(activeIds.size());
  
  typedef typename set<UInt>::iterator ITERATOR;
  UInt i=1;
  
  for(ITERATOR iter = activeIds.begin(); iter != activeIds.end(); ++iter)
  {
    outList(i) = *iter;
    ++i;
  }
  
  return(outList);
}

set<UInt>
dofMapStatic2d_options::
getIdsSet() const
{
  return(activeIds);
}


//_________________________________________________________________________________________________
// BLOCKIDS FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
dofMapStatic2d_options::
setBlockNum(const UInt & numBlocks)
{
  blockIds.resize(numBlocks);
}

void
dofMapStatic2d_options::
addBlockGeoId(const SET & Ids)
{
  blockIds.push_back(Ids);
}

void
dofMapStatic2d_options::
addBlockGeoId(const UInt & k, const UInt & Id)
{
  assert(k >= 1);
  assert(k <= blockIds.size());
  
  blockIds(k).insert(Id);
}

void
dofMapStatic2d_options::
setBlockGeoId(const UInt & k, const SET & Ids)
{
  assert(k >= 1);
  assert(k <= blockIds.size());
  
  blockIds(k) = Ids;
}

void
dofMapStatic2d_options::
setBlockGeoId(const UInt & k, const sVect<UInt> & Ids)
{
  assert(k >= 1);
  assert(k <= blockIds.size());
  
  blockIds(k).clear();
  
  for(UInt i=1; i <= Ids.size(); ++i)
  { blockIds(k).insert(Ids(i)); }
}

bool
dofMapStatic2d_options::
isBlockGeoId(const UInt & k, const UInt & Id) const
{
  assert(k >= 1);
  assert(k <= blockIds.size());
  
  return(blockIds(k).count(Id) > 0);
}

sVect<UInt>
dofMapStatic2d_options::
getBlockIdsList(const UInt & k) const
{
  assert(k >= 1);
  assert(k <= blockIds.size());
  
  typedef typename set<UInt>::iterator ITERATOR;
  
  UInt i=1;
  sVect<UInt> outList(blockIds(k).size());
  
  for(ITERATOR iter = blockIds(k).begin(); iter != blockIds(k).end(); ++iter)
  {
    outList(i) = *iter;
    ++i;
  }
  
  return(outList);
}

set<UInt>
dofMapStatic2d_options::
getBlockIdsSet(const UInt & k) const
{
  assert(k >= 1);
  assert(k <= blockIds.size());
  
  return(blockIds(k));
}

UInt
dofMapStatic2d_options::
getNumBlocks() const
{
  return(blockIds.size());
}


//_________________________________________________________________________________________________
// OUTSTREAM FUNCTIONS
//-------------------------------------------------------------------------------------------------
ostream &
operator<<(ostream& f, const dofMapStatic2d_options & O)
{
  typedef typename set<UInt>::iterator ITERATOR;

  //GeoIds-----------------------------------------------------------
  f << "GeoIds" << endl;
  
  for(ITERATOR iter = O.activeIds.begin(); iter != O.activeIds.end(); ++iter)
  { f << *iter << endl; }
  
  //BlockIds---------------------------------------------------------
  f << "BlockIds" << endl;
  
  for(UInt k=1; k <= O.blockIds.size(); ++k)
  {
    for(ITERATOR iter = O.blockIds(k).begin(); iter != O.blockIds(k).end(); ++iter)
    { f << *iter << endl; }
  }
  
  return(f);
}

