/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "polyDynamicSubCard.h"

polyDynamicSubCard::
polyDynamicSubCard()
{
}

polyDynamicSubCard::
polyDynamicSubCard(const UInt & Cx, const UInt & Cy, const UInt & Cz, const Real & C)
{
  pX = Cx;
  pY = Cy;
  pZ = Cz;
  c  = C;
}

polyDynamicSubCard::
polyDynamicSubCard(const polyDynamicSubCard & card)
{
  pX = card.pX;
  pY = card.pY;
  pZ = card.pZ;
  c  = card.c;
}
    
polyDynamicSubCard
polyDynamicSubCard::
operator=(const polyDynamicSubCard & card)
{
  pX = card.pX;
  pY = card.pY;
  pZ = card.pZ;
  c  = card.c;
  
  return(*this);
}

void
polyDynamicSubCard::
operator*=(const Real & val)
{
  c *= val;
}
    
const UInt &
polyDynamicSubCard::
getCx() const
{
  return(pX);
}

const UInt &
polyDynamicSubCard::
getCy() const
{
  return(pY);
}

const UInt &
polyDynamicSubCard::
getCz() const
{
  return(pZ);
}

const Real &
polyDynamicSubCard::
getCoeff() const
{
  return(c);
}

UInt &
polyDynamicSubCard::
getCx()
{
  return(pX);
}

UInt &
polyDynamicSubCard::
getCy()
{
  return(pY);
}

UInt &
polyDynamicSubCard::
getCz()
{
  return(pZ);
}

Real &
polyDynamicSubCard::
getCoeff()
{
  return(c);
}

void
polyDynamicSubCard::
setCx(const UInt & PX)
{
  pX = PX;
}

void
polyDynamicSubCard::
setCy(const UInt & PY)
{
  pY = PY;
}

void
polyDynamicSubCard::
setCz(const UInt & PZ)
{
  pZ = PZ;
}

void
polyDynamicSubCard::
setCoeff(const Real & Coeff)
{
  c = Coeff;
}

bool
polyDynamicSubCard::
operator<(const polyDynamicSubCard & E) const
{
  if(pX < E.pX) {return(true);}
  if(pX > E.pX) {return(false);}
  
  if(pY < E.pY) {return(true);}
  if(pY > E.pY) {return(false);}
  
  if(pZ <  E.pZ) {return(true);}
  else           {return(false);}

}
    
bool
polyDynamicSubCard::
operator!=(const polyDynamicSubCard & E) const
{
  return( !((pX == E.pX) && (pY == E.pY) && (pZ == E.pZ)) );
}

ostream & operator<<(ostream & f, const polyDynamicSubCard & G)
{
  f << "px: " << G.pX << " py: " << G.pY << " pz: " << G.pZ << " c: " << G.c << endl;
  
  return(f);
}
