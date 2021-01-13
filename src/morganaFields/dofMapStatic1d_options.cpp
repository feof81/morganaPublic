/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "dofMapStatic1d_options.h"


dofMapStatic1d_options::
dofMapStatic1d_options()
{
}

dofMapStatic1d_options::
dofMapStatic1d_options(const dofMapStatic1d_options & O)
{
  typedef typename set<UInt>::iterator ITERATOR;
  ITERATOR iter;
  
  activeIds.clear();
  
  for(iter = O.activeIds.begin(); iter != O.activeIds.end(); ++iter)
  { activeIds.insert(*iter); }
}

dofMapStatic1d_options
dofMapStatic1d_options::
operator=(const dofMapStatic1d_options & O)
{
  typedef typename set<UInt>::iterator ITERATOR;
  ITERATOR iter;
  
  activeIds.clear();
  
  for(iter = O.activeIds.begin(); iter != O.activeIds.end(); ++iter)
  { activeIds.insert(*iter); }
  
  return(*this);
}

void
dofMapStatic1d_options::
addGeoId(const UInt & Id)
{
  activeIds.insert(Id);
}

UInt
dofMapStatic1d_options::
numGeoIds() const
{
  return(activeIds.size());
}

bool
dofMapStatic1d_options::
isGeoId(const UInt & Id) const
{
  return(activeIds.count(Id) > 0);
}

void
dofMapStatic1d_options::
clear()
{
  activeIds.clear();
}

sVect<UInt>
dofMapStatic1d_options::
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
dofMapStatic1d_options::
getIdsSet() const
{
  return(activeIds);
}

ostream &
operator<<(ostream& f, const dofMapStatic1d_options & O)
{
  typedef typename set<UInt>::iterator ITERATOR;

  for(ITERATOR iter = O.activeIds.begin(); iter != O.activeIds.end(); ++iter)
  {
    f << *iter << endl;
  }
  
  return(f);
}
