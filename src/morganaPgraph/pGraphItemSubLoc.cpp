/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

morgana is free software: you can redistribute it and/or modify it under the terms of the GNU GeneralPublic License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "pGraphItemSubLoc.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
pGraphItemSubLoc::
pGraphItemSubLoc() : pGraphItem()
{ }

pGraphItemSubLoc::
pGraphItemSubLoc(const UInt & N) : pGraphItem(N)
{
  subLocId.resize(N);
}

pGraphItemSubLoc::
pGraphItemSubLoc(const sVect<UInt> & values) : pGraphItem(values)
{
  subLocId.resize(values.size());
}

pGraphItemSubLoc::
pGraphItemSubLoc(const pGraphItemSubLoc &G) : pGraphItem(G)
{
  subLocId = G.subLocId;
}



//_________________________________________________________________________________________________
// SIZING FUNCTIONS
//-------------------------------------------------------------------------------------------------
UInt 
pGraphItemSubLoc::
size() const
{
  assert(connected.size() == subLocId.size());
  return(connected.size());
}

void 
pGraphItemSubLoc::
resize(const UInt & dim)
{
  connected.resize(dim);
  subLocId.resize(dim);
}

void
pGraphItemSubLoc::
reserve(const UInt & n)
{
  assert(connected.size() == subLocId.size());
  connected.reserve(n);
  subLocId.reserve(n);
}

void
pGraphItemSubLoc::
push_back(const UInt & value, const UInt & subIndex)
{
  connected.push_back(value);
  subLocId.push_back(subIndex);
}



//_________________________________________________________________________________________________
// GET-SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
pGraphItemSubLoc::
setSubLocId(const UInt & i, const UInt & value)
{
  assert(i >= 1);
  assert(i <= subLocId.size());
  
  subLocId(i) = value;
}

void
pGraphItemSubLoc::
setCids(const sVect<UInt> & Connected)
{
  connected = Connected;
  subLocId.resize(connected.size());
}

void
pGraphItemSubLoc::
setData(const UInt & i, const UInt & Id, const UInt & subId)
{
  assert(connected.size() == subLocId.size());
  assert(i >= 1);
  assert(i<=subLocId.size());
  
  connected(i) = Id;
  subLocId(i)  = subId;
}

void
pGraphItemSubLoc::
setData(const sVect<UInt> & Connected, const sVect<UInt> & subIds)
{
  assert(Connected.size() == subIds.size());
  
  connected = Connected;
  subLocId  = subIds;
}

UInt &
pGraphItemSubLoc::
getSubLocId(const UInt & i)
{
  assert(i >= 1);
  assert(i<=subLocId.size());
  
  return(subLocId(i));
}

const UInt &
pGraphItemSubLoc::
getSubLocId(const UInt & i) const
{
  assert(i >= 1);
  assert(i<=subLocId.size());
  
  return(subLocId(i));
}



//_________________________________________________________________________________________________
// OPERATORS FUNCTIONS
//-------------------------------------------------------------------------------------------------
pGraphItemSubLoc & 
pGraphItemSubLoc::
operator=(const pGraphItemSubLoc & E)
{
  connected = E.connected;
  subLocId  = E.subLocId;
  
  return *this;
}



//_________________________________________________________________________________________________
// OUTSTREAM FUNCTIONS
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const pGraphItemSubLoc & G)
{
  f << "Num Points : " << G.size() << endl;
  f << "Nodes Id's : " << endl;
  
  for(UInt i=1; i<=G.size(); ++i)
  {
    f << "Id: " << i << " ConnectedId: " << G.getCid(i) << " SubLoc: " << G.getSubLocId(i) << endl;
  }
  f << endl;
  
  return(f);
}

size_t
pGraphItemSubLoc::
memSize() const
{
  return( sizeof(UInt) * (connected.size() + orderedNodes.size() + subLocId.size())
        + sizeof(nodesOrdered) );
}


//_________________________________________________________________________________________________
// ASSOCIATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
size_t dynamicSizeOf(const pGraphItemSubLoc & A)
{ return(A.memSize()); }
