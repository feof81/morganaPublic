/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "feHY2d_card.h"
#include <iostream>


//_________________________________________________________________________________________________
// CONSTRUCTORS AND FUNCTIONS
//-------------------------------------------------------------------------------------------------
feHY2d_card::
feHY2d_card()
{
  activeEdges.resize(3);
}

feHY2d_card::
feHY2d_card(const feHY2d_card & C)
{
  activeEdges = C.activeEdges;
}

feHY2d_card
feHY2d_card::
operator=(const feHY2d_card & C)
{
  activeEdges = C.activeEdges;
  return(*this);
}

bool
feHY2d_card::
operator!=(const feHY2d_card & C) const
{
  return(
  (activeEdges(1) == C.activeEdges(1)) &&
  (activeEdges(2) == C.activeEdges(2)) &&
  (activeEdges(3) == C.activeEdges(3)) );  
}



//_________________________________________________________________________________________________
// SET AND GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
feHY2d_card::
setActiveEdges(const sVect<bool> & ActiveEdges)
{
  activeEdges = ActiveEdges;
}

void
feHY2d_card::
setActiveEdge(const UInt & i, const bool & active)
{
  assert(i >= 1);
  assert(i <= 4);
  activeEdges(i) = active;
}

const sVect<bool> &
feHY2d_card::
getActiveEdges() const
{
  return(activeEdges);
}


//_________________________________________________________________________________________________
// PARALLEL SUPPORT
//-------------------------------------------------------------------------------------------------
std::ostream & operator<<(ostream & f, const feHY2d_card & V)
{
  f << V.activeEdges(1) << " " << V.activeEdges(2) << " " << V.activeEdges(3) << endl;
  return(f);
}
