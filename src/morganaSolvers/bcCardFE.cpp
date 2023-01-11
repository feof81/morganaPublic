/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "bcCardFE.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
bcCardFE::
bcCardFE()
{
  geoId  = 0;
  bcFlag = 0;
}

bcCardFE::
bcCardFE(const UInt & GeoId, const UInt & BcFlag)
{
  geoId  = GeoId;
  bcFlag = BcFlag;
}

bcCardFE::
bcCardFE(const bcCardFE & B)
{
  geoId  = B.geoId;
  bcFlag = B.bcFlag;
}

bcCardFE
bcCardFE::
operator=(const bcCardFE & B)
{
  geoId  = B.geoId;
  bcFlag = B.bcFlag;
  
  return(*this);
}



//_________________________________________________________________________________________________
// FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
bcCardFE::
setGeoId(const UInt & GeoId)
{
  geoId = GeoId;
}

void
bcCardFE::
setBcFlag(const UInt & BcFlag)
{
  bcFlag = BcFlag;
}

const UInt &
bcCardFE::
getGeoId() const
{
  return(geoId);
}

const UInt &
bcCardFE::
getBcFlag() const
{
  return(bcFlag);
}

UInt &
bcCardFE::
getGeoId()
{
  return(geoId);
}

UInt &
bcCardFE::
getBcFlag()
{
  return(bcFlag);
}



//_________________________________________________________________________________________________
// ORDINAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
bool
bcCardFE::
operator<(const bcCardFE & B) const
{
  return(geoId < B.geoId);
}

bool
bcCardFE::
operator!=(const bcCardFE & B) const
{
  return(geoId != B.geoId);
}



//_________________________________________________________________________________________________
// OUTSTREAM OPERATOR
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const bcCardFE & B)
{
  f << "geoId: " << B.geoId << " flag: " << B.bcFlag << endl;
  return(f);
}

