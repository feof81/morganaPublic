/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "feOscHYB3d_card.h"

feOscHYB3d_card::
feOscHYB3d_card()
{
}

feOscHYB3d_card::
feOscHYB3d_card(const UInt & N) : feOsc3d_card::feOsc3d_card(N), feHYB3d_card::feHYB3d_card()
{
}

feOscHYB3d_card::
feOscHYB3d_card(const feOscHYB3d_card & C) : feOsc3d_card::feOsc3d_card(C), feHYB3d_card::feHYB3d_card(C)
{
}

feOscHYB3d_card
feOscHYB3d_card::
operator=(const feOscHYB3d_card & C)
{
  activeFaces      = C.activeFaces;
  facesOrientation = C.facesOrientation;
  Yj = C.Yj;
  Hj = C.Hj;
  
  return(*this);
}

bool
feOscHYB3d_card::
operator!=(const feOscHYB3d_card & C) const
{
  return(feOsc3d_card::operator!=(C) || feHYB3d_card::operator!=(C));
}

//_________________________________________________________________________________________________
// PARALLEL SUPPORT
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const feOscHYB3d_card & V)
{
  assert(V.Yj.size() == V.Hj.size());
  
  f << "Osc card " << endl;
  
  f << "Active faces      : " << V.activeFaces(1)      << " " << V.activeFaces(2)      << " " << V.activeFaces(3)      << " " << V.activeFaces(4) << endl;
  f << "Faces orientation : " << V.facesOrientation(1) << " " << V.facesOrientation(2) << " " << V.facesOrientation(3) << " " << V.facesOrientation(4) << endl;
  
  for(UInt i=1; i <= V.Yj.size(); ++i)
  { 
    f << "Y" << i << ": " << V.Yj(i);
    f << "H" << i << ": " << V.Hj(i);
  }
  
  return(f);
}

