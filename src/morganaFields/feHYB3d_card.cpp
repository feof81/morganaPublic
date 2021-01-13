/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "feHYB3d_card.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS AND FUNCTIONS
//-------------------------------------------------------------------------------------------------
feHYB3d_card::
feHYB3d_card()
{
  activeFaces.resize(4);
  facesOrientation.resize(4);
}

feHYB3d_card::
feHYB3d_card(const feHYB3d_card & C)
{
  activeFaces      = C.activeFaces;
  facesOrientation = C.facesOrientation;
}

feHYB3d_card
feHYB3d_card::
operator=(const feHYB3d_card & C)
{
  activeFaces      = C.activeFaces;
  facesOrientation = C.facesOrientation;
  
  return(*this);
}

bool
feHYB3d_card::
operator!=(const feHYB3d_card & C) const
{
  return(
  (activeFaces(1) == C.activeFaces(1)) &&
  (activeFaces(2) == C.activeFaces(2)) &&
  (activeFaces(3) == C.activeFaces(3)) &&
  (activeFaces(4) == C.activeFaces(4)) &&
  
  (facesOrientation(1) == C.facesOrientation(1)) &&
  (facesOrientation(2) == C.facesOrientation(2)) &&
  (facesOrientation(3) == C.facesOrientation(3)) &&
  (facesOrientation(4) == C.facesOrientation(4)) );  
}



//_________________________________________________________________________________________________
// SET AND GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
feHYB3d_card::
setActiveFaces(const sVect<bool> & ActiveFaces)
{
  activeFaces = ActiveFaces;
}

void
feHYB3d_card::
setFacesOrientation(const sVect<bool> & FacesOrientation)
{
  facesOrientation = FacesOrientation;
}    

void
feHYB3d_card::
setActiveFace(const UInt & i, const bool & active)
{
  assert(i >= 1);
  assert(i <= 4);
  
  activeFaces(i) = active;
}

void
feHYB3d_card::
setFaceOrientation(const UInt & i, const bool & FaceOrientation)
{
  assert(i >= 1);
  assert(i <= 4);
  
  facesOrientation(i) = FaceOrientation;
}    

const sVect<bool> &
feHYB3d_card::
getActiveFaces() const
{
  return(activeFaces);
}

const sVect<bool> &
feHYB3d_card::
getFacesOrientation() const
{
  return(facesOrientation);
}    

sVect<bool> &
feHYB3d_card::
getActiveFaces()
{
  return(activeFaces);
}

sVect<bool> &
feHYB3d_card::
getFacesOrientation()
{
  return(facesOrientation);
}

bool
feHYB3d_card::
getActiveFaces(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= 4);
  
  return(activeFaces(i));
}

bool
feHYB3d_card::
getFaceOrientation(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= 4);
  
  return(facesOrientation(i));
}


//_________________________________________________________________________________________________
// PARALLEL SUPPORT
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const feHYB3d_card & V)
{
  f << "Active faces      : " << V.activeFaces(1)      << " " << V.activeFaces(2)      << " " << V.activeFaces(3)      << " " << V.activeFaces(4) << endl;
  f << "Faces orientation : " << V.facesOrientation(1) << " " << V.facesOrientation(2) << " " << V.facesOrientation(3) << " " << V.facesOrientation(4) << endl;
  
  return(f);
}
