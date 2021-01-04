/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "point3d.h"
#include <cmath>


//_________________________________________________________________________________________________
// CONSTRUCTORS AND DESTRUCTORS
//-------------------------------------------------------------------------------------------------
point3d::
point3d(Real xx, Real yy, Real zz) : id(1)
{
  X[0] = xx;
  X[1] = yy;
  X[2] = zz;
}

point3d::
point3d(const point3d & V) : id(V.id)
{
  X[0] = V.X[0];
  X[1] = V.X[1];
  X[2] = V.X[2];
}

point3d::
~point3d()
{
}


//_________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
UInt
point3d::
getId() const
{
  return(id);
}


//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
point3d::
setX(const Real  & xx)
{
  X[0] = xx;
}

void
point3d::
setY(const Real & yy)
{
  X[1] = yy;
}

void
point3d::
setZ(const Real & zz)
{
  X[2] = zz;
}

void
point3d::
setI(const UInt & i, const Real & val)
{
  assert(i<=3);
  X[i-1] = val;
}

void
point3d::
set(const Real & xx, const Real & yy, const Real & zz)
{
  X[0] = xx;
  X[1] = yy;
  X[2] = zz;
}

void
point3d::
setId(const UInt & Id)
{
  id = Id;
}



//_________________________________________________________________________________________________
// OPERATORS FUNCTIONS
//-------------------------------------------------------------------------------------------------
Real
point3d::
operator()(const UInt & i)
{
  assert(i>=1);
  assert(i<=3);
  
  return(X[i-1]);
}
    
const Real &
point3d::
operator()(const UInt & i) const
{
  assert(i>=1);
  assert(i<=3);
  
  return(X[i-1]);
}

point3d &
point3d::
operator=(const point3d & V)
{
  X[0] = V.X[0];
  X[1] = V.X[1];
  X[2] = V.X[2];
  id   = V.id;
  
  return *this;
}

void
point3d::
operator+=(const point3d & V)
{
  X[0] += V.X[0];
  X[1] += V.X[1];
  X[2] += V.X[2];
}

void
point3d::
operator-=(const point3d & V)
{
  X[0] -= V.X[0];
  X[1] -= V.X[1];
  X[2] -= V.X[2];
}

void 
point3d::
operator*=(const Real & a)
{
  X[0] *= a;
  X[1] *= a;
  X[2] *= a;
}

void
point3d::
operator/=(const Real & a)
{
  X[0] /= a;
  X[1] /= a;
  X[2] /= a;
}

Real 
point3d::
operator*(const point3d & V) const
{
  return(X[0]*V.X[0] + X[1]*V.X[1] + X[2]*V.X[2]);
}

void
point3d::
rotor(const point3d & Px, const point3d & Py, const point3d & Pz)
{
  X[0] = Py.X[2]-Pz.X[1];
  X[1] = Pz.X[0]-Px.X[2];
  X[2] = Px.X[1]-Py.X[0];
}



//_________________________________________________________________________________________________
// ORDINAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
bool 
point3d::
operator<(const point3d & V) const
{
  /*if( sqrt( pow(X[0]-V.X[0],2) + pow(X[1]-V.X[1],2) + pow(X[2]-V.X[2],2) ) <= geoToll )
  { return(false); }
  
  for(UInt i=0; i<3; ++i)
  {
    if(X[i] < V.X[i])
    { return(true); }
    
    if(X[i] > V.X[i])
    { return(false); }
  }
  
  return(false);*/
  
  
  /*for(UInt i=0; i<3; ++i)
  {
    if(X[i] < (V.X[i] - geoToll))
    { return(true); }
    
    if(X[i] > (V.X[i] + geoToll))
    { return(false); }
  }
  
  return(false);*/
  
  
  
  Real Vx = geoToll * round(V.X[0] / geoToll);
  Real Vy = geoToll * round(V.X[1] / geoToll);
  Real Vz = geoToll * round(V.X[2] / geoToll);
  
  Real Xx = geoToll * round(X[0] / geoToll);
  Real Xy = geoToll * round(X[1] / geoToll);
  Real Xz = geoToll * round(X[2] / geoToll);
  
  if(Xx < Vx) { return(true); }  if(Xx > Vx) { return(false); }
  if(Xy < Vy) { return(true); }  if(Xy > Vy) { return(false); }
  if(Xz < Vz) { return(true); }  if(Xz > Vz) { return(false); }
  
  return(false);
}

bool
point3d::
operator!=(const point3d & V) const
{
  /*bool equal[3];

  equal[0] = (X[0] >= (V.X[0] - geoToll)) && (X[0] <= (V.X[0] + geoToll));
  equal[1] = (X[1] >= (V.X[1] - geoToll)) && (X[1] <= (V.X[1] + geoToll));
  equal[2] = (X[2] >= (V.X[2] - geoToll)) && (X[2] <= (V.X[2] + geoToll));

  return(!(equal[0] && equal[1] && equal[2] ) );*/
  
  
  //return( sqrt( pow(X[0]-V.X[0],2) + pow(X[1]-V.X[1],2) + pow(X[2]-V.X[2],2) ) > geoToll );
  
  
  Real Vx = geoToll * round(V.X[0] / geoToll);
  Real Vy = geoToll * round(V.X[1] / geoToll);
  Real Vz = geoToll * round(V.X[2] / geoToll);
  
  Real Xx = geoToll * round(X[0] / geoToll);
  Real Xy = geoToll * round(X[1] / geoToll);
  Real Xz = geoToll * round(X[2] / geoToll);
  
  return( !((Vx == Xx) && (Vy == Xy) && (Vz == Xz)) );
}

bool
point3d::
operator==(const point3d & V) const
{
  return( !this->operator!=(V) );
}


//_________________________________________________________________________________________________
// COMBINATORIAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
point3d
point3d::
combinationGS(const point3d & N) const
{
  Real a = X[0] * N.X[0] + X[1] * N.X[1] + X[2] * N.X[2];
  return(point3d(X[0] - N.X[0] * a,
		 X[1] - N.X[1] * a,
		 X[2] - N.X[2] * a));
}


//_________________________________________________________________________________________________
// OTHER FUNCTIONS
//-------------------------------------------------------------------------------------------------
Real
point3d::
norm2() const
{
  return(sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]) );
}

void
point3d::
clear()
{
  X[0] = 0.0;
  X[1] = 0.0;
  X[2] = 0.0;
}

std::pair<point3d,point3d>
point3d::
orthoBasis()
{
  assert(norm2(*this) > 0.0);
  
  point3d t1(  0,   X[2], -X[1]);
  point3d t2(-X[2],   0,   X[0]);
  
  Real g1 = t1.norm2();
  Real g2 = t2.norm2();
  
  t1.X[0] = t1.X[0] * Real(g1 > g2) + t2.X[0] * Real(g2 >= g1);
  t1.X[1] = t1.X[1] * Real(g1 > g2) + t2.X[1] * Real(g2 >= g1);
  t1.X[2] = t1.X[2] * Real(g1 > g2) + t2.X[2] * Real(g2 >= g1);
  
  t2.set(-X[1] * t1.X[2] + X[2] * t1.X[1],
         -X[2] * t1.X[0] + X[0] * t1.X[2],
         -X[0] * t1.X[1] + X[1] * t1.X[0] );
  
  std::pair<point3d,point3d> out;
  out.first  = t1;
  out.second = t2;
  
  return(out);
}

ostream &
operator<<(ostream& f, const point3d & P)
{
  return f << P.X[0] << " " << P.X[1] << " " << P.X[2] << endl;
}

Real
point3d::
dot(const point3d & P1, const point3d & P2)
{
  return(
  P1.X[0] * P2.X[0] +
  P1.X[1] * P2.X[1] +
  P1.X[2] * P2.X[2]);
}

Real
point3d::
norm2(const point3d & P)
{
  return(sqrt(dot(P,P)));
}


