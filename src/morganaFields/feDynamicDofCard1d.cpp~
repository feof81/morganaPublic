/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "feDynamicDofCard2d.h"

feDynamicDofCard2d::
feDynamicDofCard2d()
{
}

feDynamicDofCard2d::
feDynamicDofCard2d(const ReferenceGeometry & GeoType, const UInt & Lev, const UInt & LocId)
{
  geoType = GeoType;
  lev     = Lev;
  locId   = LocId;
}

feDynamicDofCard2d::
feDynamicDofCard2d(const feDynamicDofCard2d & F)
{
  geoType = F.geoType;
  lev     = F.lev;
  locId   = F.locId;
}

feDynamicDofCard2d
feDynamicDofCard2d::
operator=(const feDynamicDofCard2d & F)
{
  geoType = F.geoType;
  lev     = F.lev;
  locId   = F.locId;
  
  return(*this);
}

void
feDynamicDofCard2d::
setGeoType(const ReferenceGeometry & GeoType)
{
  geoType = GeoType;
}

void
feDynamicDofCard2d::
setLevel(const UInt & Lev)
{
  lev = Lev;
}

void
feDynamicDofCard2d::
setLocalId(const UInt & LocId)
{
  locId = LocId;
}

void
feDynamicDofCard2d::
setLocalElId(const UInt & LocElId)
{
  elId = LocElId;
}
    
ReferenceGeometry &
feDynamicDofCard2d::
getGeoType()
{
  return(geoType);
}

UInt &
feDynamicDofCard2d::
getLevel()
{
  return(lev);
}

UInt &
feDynamicDofCard2d::
getLocalId()
{
  return(locId);
}

UInt &
feDynamicDofCard2d::
getLocalElId()
{
  return(elId);
}

const ReferenceGeometry &
feDynamicDofCard2d::
getGeoType() const
{
  return(geoType);
}

const UInt &
feDynamicDofCard2d::
getLevel() const
{
  return(lev);
}

const UInt &
feDynamicDofCard2d::
getLocalId() const
{
  return(locId);
}

const UInt &
feDynamicDofCard2d::
getLocalElId() const
{
  return(elId);
}

ostream & operator<<(ostream & f, const feDynamicDofCard2d & V)
{
  f << "GeoType : " << V.geoType << endl;
  f << "Level   : " << V.lev << endl;
  f << "LocId   : " << V.locId << endl;
  
  return(f);
}
