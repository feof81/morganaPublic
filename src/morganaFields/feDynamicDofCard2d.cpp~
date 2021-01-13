/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "feDynamicDofCard1d.h"

feDynamicDofCard1d::
feDynamicDofCard1d()
{
}

feDynamicDofCard1d::
feDynamicDofCard1d(const ReferenceGeometry & GeoType, const UInt & Lev, const UInt & LocId)
{
  geoType = GeoType;
  lev     = Lev;
  locId   = LocId;
}

feDynamicDofCard1d::
feDynamicDofCard1d(const feDynamicDofCard1d & F)
{
  geoType = F.geoType;
  lev     = F.lev;
  locId   = F.locId;
}

feDynamicDofCard1d
feDynamicDofCard1d::
operator=(const feDynamicDofCard1d & F)
{
  geoType = F.geoType;
  lev     = F.lev;
  locId   = F.locId;
  
  return(*this);
}

void
feDynamicDofCard1d::
setGeoType(const ReferenceGeometry & GeoType)
{
  geoType = GeoType;
}

void
feDynamicDofCard1d::
setLevel(const UInt & Lev)
{
  lev = Lev;
}

void
feDynamicDofCard1d::
setLocalId(const UInt & LocId)
{
  locId = LocId;
}

void
feDynamicDofCard1d::
setLocalElId(const UInt & LocElId)
{
  elId = LocElId;
}
    
ReferenceGeometry &
feDynamicDofCard1d::
getGeoType()
{
  return(geoType);
}

UInt &
feDynamicDofCard1d::
getLevel()
{
  return(lev);
}

UInt &
feDynamicDofCard1d::
getLocalId()
{
  return(locId);
}

UInt &
feDynamicDofCard1d::
getLocalElId()
{
  return(elId);
}

const ReferenceGeometry &
feDynamicDofCard1d::
getGeoType() const
{
  return(geoType);
}

const UInt &
feDynamicDofCard1d::
getLevel() const
{
  return(lev);
}

const UInt &
feDynamicDofCard1d::
getLocalId() const
{
  return(locId);
}

const UInt &
feDynamicDofCard1d::
getLocalElId() const
{
  return(elId);
}

ostream & operator<<(ostream & f, const feDynamicDofCard1d & V)
{
  f << "GeoType : " << V.geoType << endl;
  f << "Level   : " << V.lev << endl;
  f << "LocId   : " << V.locId << endl;
  
  return(f);
}
