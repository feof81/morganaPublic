/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "boundCard.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
boundCard::
boundCard()
{
}

boundCard::
boundCard(const set<UInt> & GeoIds, const UInt & Flag, const Real & Value)
{
  geoIds = GeoIds;
  flag   = Flag;
  value  = Value;
}

boundCard::
boundCard(const boundCard & B)
{
  geoIds = B.geoIds;
  flag   = B.flag;
  value  = B.value;
}

boundCard
boundCard::
operator=(const boundCard & B)
{
  geoIds = B.geoIds;
  flag   = B.flag;
  value  = B.value;
  
  return(*this);
}

bool
boundCard::
operator!=(const boundCard & B) const
{
  typedef set<UInt>::iterator ITER;
  
  if(geoIds.size() == B.geoIds.size())
  {
    bool flag = true;
    ITER iterExt = B.geoIds.begin();
    
    for(ITER iter = geoIds.begin(); iter != geoIds.end(); iter++)
    {
      flag = flag && (*iter == *iterExt);    
      iterExt++;
    }

    return(!flag);
  }
  else
  { return(true); }
}

bool
boundCard::
operator==(const boundCard & B) const
{
  return(!operator!=(B));
}


//_________________________________________________________________________________________________
// FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
boundCard::
setGeoIds(const set<UInt> & GeoIds)
{
  geoIds = GeoIds;
}

void
boundCard::
addGeoId(const UInt & GeoId)
{
  geoIds.insert(GeoId);
}

void
boundCard::
resetGeoIds()
{
  geoIds.clear(); 
}

void
boundCard::
setValue(const Real & Value)
{
  value = Value;
}

void
boundCard::
setFlag(const UInt & Flag)
{
  flag = Flag;
}

bool
boundCard::
isGeoId(const UInt & GeoId)
{
  return(geoIds.count(GeoId) == 1);
}

const set<UInt> &
boundCard::
getGeoIds() const
{
  return(geoIds);
}

sVect<UInt>
boundCard::
getVectGeoIds() const
{
  typedef set<UInt>::iterator ITER;
  ITER iter;
  
  sVect<UInt> list;
  
  for(iter = geoIds.begin(); iter != geoIds.end(); ++iter)
  { list.push_back(*iter); }
  
  return(list);
}

Real &
boundCard::
getValue()
{
  return(value);
}

const Real &
boundCard::
getValue() const
{
  return(value);
}

UInt &
boundCard::
getFlag()
{
  return(flag);
}

const UInt &
boundCard::
getFlag() const
{
  return(flag);
}


//_________________________________________________________________________________________________
// OUTSTREAM
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const boundCard & B)
{
  typedef set<UInt>::iterator ITERATOR;
  
  f << "geoIds: ";
  
  for(ITERATOR iter = B.geoIds.begin(); iter != B.geoIds.end(); iter++)
  {
    f << *iter << " ";
  }
  
  f << endl;
  f << "flag  : " << B.flag << endl;
  f << "Value : " << B.value << endl << endl;
  return(f);
}



//_________________________________________________________________________________________________
// FUNCTIONS
//-------------------------------------------------------------------------------------------------
sVect<UInt> boundCard_functions::boundList(const sVect<boundCard> & cards)
{
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> out, temp;
  sVect<UInt> final;
  
  for(UInt i=1; i <= cards.size(); ++i)
  {
    temp = cards(i).getGeoIds();
    
    for(ITERATOR iter = temp.begin(); iter != temp.end(); iter++)
    { out.insert(*iter); }
  }
  
  final.reserve(out.size());
  
  for(ITERATOR iter = out.begin(); iter != out.end(); iter++)
  { final.push_back(*iter); }
  
  return(final);
}

sVect<UInt> boundCard_functions::unBoundList(const sVect<UInt> mainList, const sVect<boundCard> & cards)
{
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> out, temp;
  sVect<UInt> final;
  
  for(UInt i=1; i <= cards.size(); ++i)
  {
    temp = cards(i).getGeoIds();
    
    for(ITERATOR iter = temp.begin(); iter != temp.end(); iter++)
    { out.insert(*iter); }
  }
  
  for(UInt i=1; i <= mainList.size(); ++i)
  {
    if(out.count(mainList(i)) != 1)
    { final.push_back(mainList(i)); }
  }
  
  return(final);
}

bool boundCard_functions::isCompatible(const sVect<UInt> mainList, const sVect<boundCard> & cards)
{
  //Set of the original list
  set<UInt> mainSet;
  
  for(UInt i=1; i <= mainList.size(); ++i)
  { mainSet.insert(mainList(i)); }
  
  //Sublist
  sVect<UInt> bList = boundCard_functions::boundList(cards);
  
  //Compatibility
  bool flag = true;
  
  for(UInt i=1; i <= bList.size(); ++i)
  { flag = flag & (mainSet.count(bList(i)) == 1); }
  
  return(flag);
}




