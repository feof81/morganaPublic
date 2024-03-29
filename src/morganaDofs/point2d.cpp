/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "point2d.h"
#include <cmath>


//_________________________________________________________________________________________________
// CONSTRUCTORS AND DESTRUCTORS
//-------------------------------------------------------------------------------------------------
point2d::
point2d(Real xx, Real yy) : id(1)
{
  X[0] = xx;
  X[1] = yy;
}

point2d::
point2d(const point2d & V) : id(V.id)
{
  X[0] = V.X[0];
  X[1] = V.X[1];
}

point2d::
~point2d()
{
}


//_________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
UInt
point2d::
getId() const
{
  return(id);
}


//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
point2d::
setX(const Real  & xx)
{
  X[0] = xx;
}

void
point2d::
setY(const Real & yy)
{
  X[1] = yy;
}

void
point2d::
setI(const UInt & i, const Real & val)
{
  assert(i>=1);
  assert(i<=2);
  
  X[i-1] = val;
}

void
point2d::
setId(const UInt & Id)
{
  id = Id;
}



//_________________________________________________________________________________________________
// OPERATORS FUNCTIONS
//-------------------------------------------------------------------------------------------------
Real
point2d::
operator()(const UInt & i)
{
  assert(i>=1);
  assert(i<=2);
  
  return(X[i-1]);
}
    
const Real &
point2d::
operator()(const UInt & i) const
{
  assert(i>=1);
  assert(i<=2);
  
  return(X[i-1]);
}

point2d &
point2d::
operator=(const point2d &V)
{
  X[0] = V.X[0];
  X[1] = V.X[1];
  id   = V.id;
  
  return *this;
}

void
point2d::
operator+=(const point2d & V)
{
  X[0] += V.X[0];
  X[1] += V.X[1];
}

void
point2d::
operator-=(const point2d & V)
{
  X[0] -= V.X[0];
  X[1] -= V.X[1];
}

void 
point2d::
operator*=(const Real & a)
{
  X[0] *= a;
  X[1] *= a;
}

void
point2d::
operator/=(const Real & a)
{
  X[0] /= a;
  X[1] /= a;
}

Real 
point2d::
operator*(const point2d & V) const
{
  return(X[0]*V.X[0] + X[1]*V.X[1]);
}



//_________________________________________________________________________________________________
// ORDINAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
bool 
point2d::
operator<(const point2d & V) const
{
  /*for(UInt i=0; i<2; ++i)
  {
    if(X[i] < (V.X[i] - geoToll))
    { return(true); }
    
    if(X[i] > (V.X[i] + geoToll))
    { return(false); }
  }
  
  return(false);*/
  
  Real Vx = geoToll * round(V.X[0] / geoToll);
  Real Vy = geoToll * round(V.X[1] / geoToll);
  
  Real Xx = geoToll * round(X[0] / geoToll);
  Real Xy = geoToll * round(X[1] / geoToll);
  
  if(Xx < Vx) { return(true); }  if(Xx > Vx) { return(false); }
  if(Xy < Vy) { return(true); }  if(Xy > Vy) { return(false); }
  
  return(false);
}

bool
point2d::
operator!=(const point2d & V) const
{
  /*bool equal[2];

  equal[0] = (X[0] >= (V.X[0] - geoToll)) && (X[0] <= (V.X[0] + geoToll));
  equal[1] = (X[1] >= (V.X[1] - geoToll)) && (X[1] <= (V.X[1] + geoToll));

  return(!(equal[0] && equal[1] ) );*/
  
  Real Vx = geoToll * round(V.X[0] / geoToll);
  Real Vy = geoToll * round(V.X[1] / geoToll);
  
  Real Xx = geoToll * round(X[0] / geoToll);
  Real Xy = geoToll * round(X[1] / geoToll);
  
  return( !((Vx == Xx) && (Vy == Xy)) );
}


//_________________________________________________________________________________________________
// OTHER FUNCTIONS
//-------------------------------------------------------------------------------------------------
Real
point2d::
norm2() const
{
  return(sqrt(X[0]*X[0] + X[1]*X[1]) );
}

void
point2d::
clear()
{
  X[0] = 0.0;
  X[1] = 0.0;
}

ostream &
operator<<(ostream& f, const point2d & P)
{
  return f << P.X[0] << " " << P.X[1] << endl;
}

Real
point2d::
dot(const point2d & P1, const point2d & P2)
{
  return(
  P1.X[0] * P2.X[0] +
  P1.X[1] * P2.X[1]);
}

Real
point2d::
norm2(const point2d & P)
{
  return(sqrt(dot(P,P)));
}
    
size_t
point2d::
memSize() const
{
  return(sizeof(point2d));
}
