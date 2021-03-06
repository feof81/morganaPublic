/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "feStaticDofCard1d.h"

feStaticDofCard1d::
feStaticDofCard1d()
{
}

feStaticDofCard1d::
feStaticDofCard1d(const ReferenceGeometry & GeoType, const UInt & Lev, const UInt & LocId)
{
  geoType = GeoType;
  lev     = Lev;
  locId   = LocId;
}

feStaticDofCard1d::
feStaticDofCard1d(const feStaticDofCard1d & F)
{
  geoType = F.geoType;
  lev     = F.lev;
  locId   = F.locId;
} 

feStaticDofCard1d
feStaticDofCard1d::
operator=(const feStaticDofCard1d & F)
{
  geoType = F.geoType;
  lev     = F.lev;
  locId   = F.locId;
  
  return(*this);
}

void
feStaticDofCard1d::
setGeoType(const ReferenceGeometry & GeoType)
{
  geoType = GeoType;
}

void
feStaticDofCard1d::
setLevel(const UInt & Lev)
{
  lev = Lev;
}

void
feStaticDofCard1d::
setLocalId(const UInt & LocId)
{
  locId = LocId;
}

void
feStaticDofCard1d::
setLocalElId(const UInt & LocElId)
{
  elId = LocElId;
}
    
ReferenceGeometry &
feStaticDofCard1d::
getGeoType()
{
  return(geoType);
}

UInt &
feStaticDofCard1d::
getLevel()
{
  return(lev);
}

UInt &
feStaticDofCard1d::
getLocalId()
{
  return(locId);
}

UInt &
feStaticDofCard1d::
getLocalElId()
{
  return(elId);
}

const ReferenceGeometry &
feStaticDofCard1d::
getGeoType() const
{
  return(geoType);
}

const UInt &
feStaticDofCard1d::
getLevel() const
{
  return(lev);
}

const UInt &
feStaticDofCard1d::
getLocalId() const
{
  return(locId);
}

const UInt &
feStaticDofCard1d::
getLocalElId() const
{
  return(elId);
}

ostream & operator<<(ostream & f, const feStaticDofCard1d & V)
{
  f << "GeoType : " << V.geoType << endl;
  f << "Level   : " << V.lev << endl;
  f << "LocId   : " << V.locId << endl;
  
  return(f);
}
