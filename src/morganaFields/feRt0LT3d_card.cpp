/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "feRt0LT3d_card.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
feRt0LT3d_card::
feRt0LT3d_card()
{
  loaded = false;
  facesOrientation.resize(4);
}

feRt0LT3d_card::
feRt0LT3d_card(const sVect<bool> & FacesOrientation)
{
  loaded = true;
  facesOrientation = FacesOrientation;
}

feRt0LT3d_card::
feRt0LT3d_card(const feRt0LT3d_card & card)
{
  loaded           = card.loaded;
  facesOrientation = card.facesOrientation;
}

feRt0LT3d_card
feRt0LT3d_card::
operator=(const feRt0LT3d_card & card)
{
  loaded           = card.loaded;
  facesOrientation = card.facesOrientation;
  
  return(*this);
}



//_________________________________________________________________________________________________
// SET AND GET
//-------------------------------------------------------------------------------------------------
void
feRt0LT3d_card::
setFacesOrientation(const sVect<bool> & FacesOrientation)
{
  loaded = true;
  facesOrientation = FacesOrientation;
}

void
feRt0LT3d_card::
setFaceOrientation(const UInt & i, const bool & FaceOrientation)
{
  loaded = true;
  facesOrientation(i) = FaceOrientation;
}

const sVect<bool> &
feRt0LT3d_card::
getFacesOrientation() const
{
  assert(loaded);
  return(facesOrientation);
}

sVect<bool> &
feRt0LT3d_card::
getFacesOrientation()
{
  assert(loaded);
  return(facesOrientation);
}

bool
feRt0LT3d_card::
getFaceOrientation(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= 4);
  assert(loaded);
  
  return(facesOrientation(i));
}



//_________________________________________________________________________________________________
// PARALLEL SUPPORT
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
feRt0LT3d_card::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & loaded;
  ar & facesOrientation;
}



//_________________________________________________________________________________________________
// PRINTOUT
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const feRt0LT3d_card & V)
{
  f << "loaded : " << V.loaded << endl;
  f << "Faces  : " << endl;
  f << V.facesOrientation << endl;
  
  return(f);
}
