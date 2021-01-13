/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "feHY3d_card.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS AND FUNCTIONS
//-------------------------------------------------------------------------------------------------
feHY3d_card::
feHY3d_card()
{
  activeFaces.resize(4);
}

feHY3d_card::
feHY3d_card(const feHY3d_card & C)
{
  activeFaces = C.activeFaces;
}

feHY3d_card
feHY3d_card::
operator=(const feHY3d_card & C)
{
  activeFaces = C.activeFaces;
  return(*this);
}

bool
feHY3d_card::
operator!=(const feHY3d_card & C) const
{
  return(
  (activeFaces(1) == C.activeFaces(1)) &&
  (activeFaces(2) == C.activeFaces(2)) &&
  (activeFaces(3) == C.activeFaces(3)) &&
  (activeFaces(4) == C.activeFaces(4)) );  
}



//_________________________________________________________________________________________________
// SET AND GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
feHY3d_card::
setActiveFaces(const sVect<bool> & ActiveFaces)
{
  activeFaces = ActiveFaces;
}

void
feHY3d_card::
setActiveFace(const UInt & i, const bool & active)
{
  assert(i >= 1);
  assert(i <= 4);
  activeFaces(i) = active;
}

const sVect<bool> &
feHY3d_card::
getActiveFaces() const
{
  return(activeFaces);
}



//_________________________________________________________________________________________________
// PARALLEL SUPPORT
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const feHY3d_card & V)
{
  f << V.activeFaces(1) << " " << V.activeFaces(2) << " " << V.activeFaces(3) << " " << V.activeFaces(4) << endl;
  return(f);
}
