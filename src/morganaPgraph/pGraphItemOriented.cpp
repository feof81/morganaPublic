/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "pGraphItemOriented.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
pGraphItemOriented::
pGraphItemOriented() : pGraphItem()
{ }

pGraphItemOriented::
pGraphItemOriented(const UInt & N) : pGraphItem(N)
{
  orientation.resize(N);
}

pGraphItemOriented::
pGraphItemOriented(const sVect<UInt> & values) : pGraphItem(values)
{
  orientation.resize(values.size());
}

pGraphItemOriented::
pGraphItemOriented(const pGraphItemOriented & G) : pGraphItem(G)
{
  orientation = G.orientation;
}



//_________________________________________________________________________________________________
// SIZING FUNCTIONS
//-------------------------------------------------------------------------------------------------
UInt 
pGraphItemOriented::
size() const
{
  assert(connected.size() == orientation.size());
  return(connected.size());
}

void
pGraphItemOriented::
resize(const UInt & dim)
{
  assert(connected.size() == orientation.size());
  connected.resize(dim);
  orientation.resize(dim);
}

void
pGraphItemOriented::
reserve(const UInt & n)
{
  assert(connected.size() == orientation.size());
  connected.reserve(n);
  orientation.reserve(n);
}

void
pGraphItemOriented::
push_back(const UInt & value, bool orient)
{
  assert(connected.size() == orientation.size());
  connected.push_back(value);
  orientation.push_back(orient);
}



//_________________________________________________________________________________________________
// GET-SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
pGraphItemOriented::
setOrientation(const UInt & i, const bool & value)
{
  assert(connected.size() == orientation.size());
  assert(i >= 1);
  assert(i<=orientation.size());
  
  orientation(i) = value;
}

void
pGraphItemOriented::
setCids(const sVect<UInt> & Connected)
{
  connected = Connected;
  orientation.resize(connected.size());
}
    
void
pGraphItemOriented::
setData(const UInt & i, const UInt & Id, const bool & logic)
{
  assert(connected.size() == orientation.size());
  assert(i >= 1);
  assert(i<=orientation.size());
  
  connected(i)   = Id;
  orientation(i) = logic;
}

void
pGraphItemOriented::
setData(const sVect<UInt> & Connected, const sVect<bool> & logics)
{
  assert(Connected.size() == logics.size());
  
  connected   = Connected;
  orientation = logics;
}

vector<bool>::const_reference
pGraphItemOriented::
getOrientation(const UInt & i) const 
{
  assert(connected.size() == orientation.size());
  assert(i >= 1);
  assert(i <= orientation.size());
  
  return(orientation(i));
}




//_________________________________________________________________________________________________
// OPERATORS FUNCTIONS
//-------------------------------------------------------------------------------------------------
pGraphItemOriented &
pGraphItemOriented::
operator=(const pGraphItemOriented & E)
{
  connected   = E.connected;
  orientation = E.orientation;
  
  return *this;
}



//_________________________________________________________________________________________________
// OUTSTREAM FUNCTIONS
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const pGraphItemOriented & G)
{
  f << "Num Points : " << G.size() << endl;
  f << "Nodes Id's : " << endl;
  
  for(UInt i=1; i <= G.size(); ++i)
  {
    f << "Id: " << i << " " << G.getCid(i) << "; " << "Orient: " << G.getOrientation(i) << endl;
  }
  f << endl;
  
  return(f);
}
