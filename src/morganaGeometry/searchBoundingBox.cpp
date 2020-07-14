/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "searchBoundingBox.h"
#include "morganaTypes.hpp"
#include "sOctTreeItem.h"


searchBoundingBox::
searchBoundingBox()
{
}

searchBoundingBox::
searchBoundingBox(const searchBoundingBox & B)
{
  elements = B.elements;
  Pmin     = B.Pmin;
  Pmax     = B.Pmax;
}

searchBoundingBox &
searchBoundingBox::
operator=(const searchBoundingBox & B)
{
  elements = B.elements;
  Pmin     = B.Pmin;
  Pmax     = B.Pmax;
  
  return(*this);
}
    
void
searchBoundingBox::
setElements(const sVect<UInt> & Elements)
{
  elements = Elements;
}

void
searchBoundingBox::
setBoundingBox(const point3d & PMIN, const point3d & PMAX)
{
  Pmin = PMIN;
  Pmax = PMAX;
}
    
const sVect<UInt> &
searchBoundingBox::
getElements() const
{
  return(elements);
}

const point3d &
searchBoundingBox::
getPmin() const
{
  return(Pmin);
}

const point3d &
searchBoundingBox::
getPmax() const
{
  return(Pmax);
}

bool
searchBoundingBox::
isInternal(const point3d & P)
{
  return(  
  (P.getX() >= (Pmin.getX() - geoToll)) &
  (P.getY() >= (Pmin.getY() - geoToll)) &
  (P.getZ() >= (Pmin.getZ() - geoToll)) &
  (P.getX() <= (Pmax.getX() + geoToll)) &
  (P.getY() <= (Pmax.getY() + geoToll)) &
  (P.getZ() <= (Pmax.getZ() + geoToll)) );
}
    
point3d
searchBoundingBox::
getPmin(const bool & ix, const bool & iy, const bool & iz) const
{
  point3d Pm = (Pmin + Pmax) / 2.0;
  Real X = Pmin.getX() * Real(!ix) + Pm.getX() * Real(ix);
  Real Y = Pmin.getY() * Real(!iy) + Pm.getY() * Real(iy);
  Real Z = Pmin.getZ() * Real(!iz) + Pm.getZ() * Real(iz);
  
  return(point3d(X,Y,Z));
}

point3d
searchBoundingBox::
getPmax(const bool & ix, const bool & iy, const bool & iz) const
{
  point3d Pm = (Pmin + Pmax) / 2.0;
  Real X = Pm.getX() * Real(!ix) + Pmax.getX() * Real(ix);
  Real Y = Pm.getY() * Real(!iy) + Pmax.getY() * Real(iy);
  Real Z = Pm.getZ() * Real(!iz) + Pmax.getZ() * Real(iz);
  
  return(point3d(X,Y,Z));
}

UInt
searchBoundingBox::
octPart(const point3d & P)
{
  point3d Pm = (Pmin + Pmax) / 2.0;
  
  bool ix = P.getX() >= Pm.getX();
  bool iy = P.getY() >= Pm.getY();
  bool iz = P.getZ() >= Pm.getZ();
  
  return(sOctTreeItem::octMap(ix,iy,iz) + 1);
}

void
searchBoundingBox::
octPart(bool & ix, bool & iy, bool & iz, const point3d & P)
{
  point3d Pm = (Pmin + Pmax) / 2.0;
  
  ix = P.getX() >= Pm.getX();
  iy = P.getY() >= Pm.getY();
  iz = P.getZ() >= Pm.getZ();
}

std::set<UInt>
searchBoundingBox::
octPart(const point3d & Bmin, const point3d & Bmax)
{
  UInt nx, ny, nz;
  bool ix, iy, iz;
  sVect<bool> indX(2), indY(2), indZ(2);
  std::set<UInt> out;
  point3d Pm = (Pmin + Pmax) / 2.0;
  
  ix = Bmin.getX() >= (Pm.getX() + geoToll);
  iy = Bmin.getY() >= (Pm.getY() + geoToll);
  iz = Bmin.getZ() >= (Pm.getZ() + geoToll);
  indX(1) = ix; indY(1) = iy; indZ(1) = iz;
  
  ix = Bmax.getX() >= (Pm.getX() - geoToll);
  iy = Bmax.getY() >= (Pm.getY() - geoToll);
  iz = Bmax.getZ() >= (Pm.getZ() - geoToll);
  indX(2) = ix; indY(2) = iy; indZ(2) = iz;
  
  nx = 1 + UInt(indX(1) != indX(2));
  ny = 1 + UInt(indY(1) != indY(2));
  nz = 1 + UInt(indZ(1) != indZ(2));
  
  for(UInt i=1; i <= nx; ++i)
  {
    for(UInt j=1; j <= ny; ++j)
    {
      for(UInt k=1; k <= nz; ++k)
      {
        out.insert( sOctTreeItem::octMap(indX(i), indY(j), indZ(k)) + 1 );
      }
    }
  }
  
  return(out);
}

Real
searchBoundingBox::
getMinDist(const point3d & P)
{
  point3d D;
  
  D.setX( (P.getX() > Pmax.getX()) * (P.getX() - Pmax.getX())  +  (P.getX() < Pmin.getX()) * (Pmin.getX() - P.getX()) );
  D.setY( (P.getY() > Pmax.getY()) * (P.getY() - Pmax.getY())  +  (P.getY() < Pmin.getY()) * (Pmin.getY() - P.getY()) );
  D.setZ( (P.getZ() > Pmax.getZ()) * (P.getZ() - Pmax.getZ())  +  (P.getZ() < Pmin.getZ()) * (Pmin.getZ() - P.getZ()) );
  
  return(D.norm2());
}

Real
searchBoundingBox::
getMaxDist(const point3d & P)
{
  point3d D;
  
  D.setX( std::max( abs(P.getX() - Pmax.getX()), abs(P.getX() - Pmin.getX()) ) );
  D.setY( std::max( abs(P.getY() - Pmax.getY()), abs(P.getY() - Pmin.getY()) ) );
  D.setZ( std::max( abs(P.getZ() - Pmax.getZ()), abs(P.getZ() - Pmin.getZ()) ) );
  
  return(D.norm2());
}

ostream & operator<<(ostream & f, const searchBoundingBox & B)
{
  f << "Elements: " << B.elements << endl;
  f << "Pmin    : " << B.Pmin << endl;
  f << "Pmax    : " << B.Pmax << endl;
  
  return(f);
}
