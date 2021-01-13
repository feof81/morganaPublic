/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "feRt0LT2d_card.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
feRt0LT2d_card::
feRt0LT2d_card()
{
  loaded = false;
  edgesOrientation.resize(3);
}

feRt0LT2d_card::
feRt0LT2d_card(const sVect<bool> & EdgesOrientation)
{
  loaded = true;
  edgesOrientation = EdgesOrientation;
}

feRt0LT2d_card::
feRt0LT2d_card(const feRt0LT2d_card & card)
{
  loaded           = card.loaded;
  edgesOrientation = card.edgesOrientation;
}

feRt0LT2d_card
feRt0LT2d_card::
operator=(const feRt0LT2d_card & card)
{
  loaded           = card.loaded;
  edgesOrientation = card.edgesOrientation;
  
  return(*this);
}



//_________________________________________________________________________________________________
// SET AND GET
//-------------------------------------------------------------------------------------------------
void
feRt0LT2d_card::
setEdgesOrientation(const sVect<bool> & EdgesOrientation)
{
  loaded = true;
  edgesOrientation = EdgesOrientation;
}

void
feRt0LT2d_card::
setEdgeOrientation(const UInt & i, const bool & EdgeOrientation)
{
  loaded = true;
  edgesOrientation(i) = EdgeOrientation;
}

const sVect<bool> &
feRt0LT2d_card::
getEdgesOrientation() const
{
  assert(loaded);
  return(edgesOrientation);
}

sVect<bool> &
feRt0LT2d_card::
getEdgesOrientation()
{
  assert(loaded);
  return(edgesOrientation);
}

bool
feRt0LT2d_card::
getEdgeOrientation(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= 3);
  assert(loaded);
  
  return(edgesOrientation(i));
}


//_________________________________________________________________________________________________
// PARALLEL SUPPORT
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
feRt0LT2d_card::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & loaded;
  ar & edgesOrientation;
}



//_________________________________________________________________________________________________
// PRINTOUT
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const feRt0LT2d_card & V)
{
  f << "loaded : " << V.loaded << endl;
  f << "Edges  : " << endl;
  f << V.edgesOrientation << endl;
  
  return(f);
}
