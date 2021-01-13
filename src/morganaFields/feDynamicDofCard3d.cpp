/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "feDynamicDofCard3d.h"

feDynamicDofCard3d::
feDynamicDofCard3d()
{
}

feDynamicDofCard3d::
feDynamicDofCard3d(const ReferenceGeometry & GeoType, const UInt & Lev, const UInt & LocId)
{
  geoType = GeoType;
  lev     = Lev;
  locId   = LocId;
}

feDynamicDofCard3d::
feDynamicDofCard3d(const feDynamicDofCard3d & F)
{
  geoType = F.geoType;
  lev     = F.lev;
  locId   = F.locId;
}

feDynamicDofCard3d
feDynamicDofCard3d::
operator=(const feDynamicDofCard3d & F)
{
  geoType = F.geoType;
  lev     = F.lev;
  locId   = F.locId;
  
  return(*this);
}

void
feDynamicDofCard3d::
setGeoType(const ReferenceGeometry & GeoType)
{
  geoType = GeoType;
}

void
feDynamicDofCard3d::
setLevel(const UInt & Lev)
{
  lev = Lev;
}

void
feDynamicDofCard3d::
setLocalId(const UInt & LocId)
{
  locId = LocId;
}

void
feDynamicDofCard3d::
setLocalElId(const UInt & LocElId)
{
  elId = LocElId;
}
    
ReferenceGeometry &
feDynamicDofCard3d::
getGeoType()
{
  return(geoType);
}

UInt &
feDynamicDofCard3d::
getLevel()
{
  return(lev);
}

UInt &
feDynamicDofCard3d::
getLocalId()
{
  return(locId);
}

UInt &
feDynamicDofCard3d::
getLocalElId()
{
  return(elId);
}

const ReferenceGeometry &
feDynamicDofCard3d::
getGeoType() const
{
  return(geoType);
}

const UInt &
feDynamicDofCard3d::
getLevel() const
{
  return(lev);
}

const UInt &
feDynamicDofCard3d::
getLocalId() const
{
  return(locId);
}

const UInt &
feDynamicDofCard3d::
getLocalElId() const
{
  return(elId);
}

ostream & operator<<(ostream & f, const feDynamicDofCard3d & V)
{
  f << "GeoType : " << V.geoType << endl;
  f << "Level   : " << V.lev << endl;
  f << "LocId   : " << V.locId << endl;
  
  return(f);
}
