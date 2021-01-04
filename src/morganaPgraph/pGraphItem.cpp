/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "pGraphItem.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
pGraphItem::
pGraphItem()
{
  updateSorting();
}

pGraphItem::
pGraphItem(UInt n)
{
  connected.resize(n);
  updateSorting();
}

pGraphItem::
pGraphItem(const sVect<UInt> & values)
{
  connected = values;
  updateSorting();
}

pGraphItem::
pGraphItem(const pGraphItem & G)
{
  connected = G.connected;
  updateSorting();
}



//_________________________________________________________________________________________________
// SIZING FUNCTIONS
//-------------------------------------------------------------------------------------------------
UInt
pGraphItem::
size() const
{
  return(connected.size());
}

void
pGraphItem::
resize(const UInt & dim)
{
  connected.resize(dim);
}

void
pGraphItem::
reserve(const UInt & n)
{
  connected.reserve(n);
}

void
pGraphItem::
push_back(const UInt & value, const bool & update)
{
  connected.push_back(value);
  nodesOrdered = update;
  if(update) {updateSorting();}
}

void
pGraphItem::
merge(const pGraphItem & G)
{
  typedef set<UInt>::iterator ITERATOR;

  //Alloc
  set<UInt> tempSet;
  ITERATOR tempIter;
  pair<ITERATOR,bool> pair;

  //Local nodes
  for(UInt i=1; i <= connected.size(); ++i)
  { tempSet.insert(connected(i)); }
  
  //Other nodes
  for(UInt i=1; i <= G.connected.size(); ++i)
  { 
    pair = tempSet.insert(G.connected(i));
    
    if(pair.second) 
    { connected.push_back(G.connected(i)); }
  }
  
  //Updating
  updateSorting();
}

void
pGraphItem::
clear()
{
  connected.clear();
  orderedNodes.clear();
  nodesOrdered = false;
}



//_________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
UInt &
pGraphItem::
getCid(const UInt & i)
{
  assert(i >= 1);
  assert(i <= connected.size());
  
  nodesOrdered = false;
  
  return(connected(i));
}

const UInt &
pGraphItem::
getCid(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= connected.size());
  return(connected(i));
}

const sVect<UInt> &
pGraphItem::
getCids() const
{
  return(connected);
}



//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
pGraphItem::
setCid(const UInt & i, const UInt & value)
{
  assert(i >= 1);
  assert(i <= connected.size());
  
  connected(i) = value;
  updateSorting();
}

void
pGraphItem::
setCids(const sVect<UInt> & Connected)
{
  connected = Connected;
  updateSorting();
}



//_________________________________________________________________________________________________
// NODE SORTING FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
pGraphItem::
updateSorting()
{
  nodesOrdered = true;

  set<UInt> tempSet;
  set<UInt>::iterator tempIter;

  for(UInt i=1; i<=connected.size(); ++i)
  { tempSet.insert(connected(i)); }
  
  UInt k=1;
  orderedNodes.resize(connected.size());
   
  for(tempIter = tempSet.begin(); tempIter != tempSet.end(); ++tempIter)
  { 
    orderedNodes(k) = *tempIter;
    ++k;
  }
}
    
const UInt &
pGraphItem::
getSorted(const UInt & i) const
{
   assert(i >= 1);
   assert(i <= connected.size());
   assert(connected.size() == orderedNodes.size());

   assert(nodesOrdered);
   return(orderedNodes(i));
}

const bool &
pGraphItem::
getNodesOrdered() const
{
  return(nodesOrdered);
}


//_________________________________________________________________________________________________
// OPERATOR FUNCTIONS
//-------------------------------------------------------------------------------------------------
pGraphItem &
pGraphItem::
operator=(const pGraphItem & E)
{
  connected = E.connected;
  
  updateSorting();
  return *this;
}

bool
pGraphItem::
operator<(const pGraphItem & E) const
{
  assert(nodesOrdered);
  
  assert(E.orderedNodes.size() == orderedNodes.size());
  assert(orderedNodes.size()   == connected.size());
  
  if(E.connected.size() != connected.size())
  {
  }
  
  //Selection____________________________________________________________________
  for(UInt i=1; i <= orderedNodes.size(); ++i)
  {
    if(orderedNodes(i) < E.orderedNodes(i))
    { return(true); }
    
    if(orderedNodes(i) > E.orderedNodes(i))
    { return(false); }
  }
  
  return(false);
}

bool
pGraphItem::
operator!=(const pGraphItem & E) const
{
  bool flag = true;
  
  assert(nodesOrdered);
  
  assert(E.connected.size()    == connected.size());
  assert(E.orderedNodes.size() == orderedNodes.size());
  assert(orderedNodes.size()   == connected.size());
  
  
  //Selection___________________________________________________________________
  for(UInt i=1; i <= orderedNodes.size(); ++i)
  { flag = flag && (orderedNodes(i) == E.orderedNodes(i)); }
  
  return(!flag);
}


//_________________________________________________________________________________________________
// OUTSTREAM OPERATOR
//-------------------------------------------------------------------------------------------------
void
pGraphItem::
printSorted()
{
  cout << "ordered ids" << endl;
  for(UInt i=1; i <= orderedNodes.size(); ++i)
  {
    cout << orderedNodes(i) << " ";
  }
  cout << endl;
}

ostream & operator<<(ostream & f, const pGraphItem & G)
{
  f << "Num Connected  : " << G.size() << endl;
  f << "Connected Id's : ";
  
  for(UInt i=1; i <= G.size(); ++i)
  {
    f << G(i) << " ";
  }
  f << endl;
  
  return(f);
}


