/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "feOsc3d_card.h"

//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
feOsc3d_card::
feOsc3d_card()
{
}

feOsc3d_card::
feOsc3d_card(const UInt & N)
{
  Yj.resize(N);
  Hj.resize(N);
}

feOsc3d_card::
feOsc3d_card(const feOsc3d_card & C)
{
  Yj = C.Yj;
  Hj = C.Hj;
}

feOsc3d_card
feOsc3d_card::
operator=(const feOsc3d_card & C)
{
  Yj = C.Yj;
  Hj = C.Hj;
  
  return(*this);
}

bool
feOsc3d_card::
operator!=(const feOsc3d_card & C) const
{
  assert(Yj.size() == Hj.size());
  
  bool flag = true;
  
  for(UInt i=1; i <= Yj.size(); ++i)
  {
    flag = flag & (!(Yj(i) != C.Yj(i)));
    flag = flag & (!(Hj(i) != C.Hj(i)));
  }
  
  return(!flag);
}


//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
feOsc3d_card::
set(const sVect<point3d> & YJ, const sVect<point3d> & HJ)
{
  assert(YJ.size() == HJ.size());
  
  Yj = YJ;
  Hj = HJ;
}

void
feOsc3d_card::
set(const point3d & Y, const point3d & H, const UInt & i)
{
  assert(Yj.size() == Hj.size());
  assert(i <= Yj.size());
  assert(i >= 1);
  
  Yj(i) = Y;
  Hj(i) = H;
}

void
feOsc3d_card::
resize(const UInt & N)
{
  assert(Yj.size() == Hj.size());
  
  Yj.resize(N);
  Hj.resize(N);
}

void
feOsc3d_card::
conjugate()
{
  assert(Yj.size() == Hj.size());
  
  for(UInt i=1; i <= Hj.size(); ++i)
  { Hj(i) *= (-1.0); }
}

//_________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
const sVect<point3d> &
feOsc3d_card::
getY() const
{
  assert(Yj.size() == Hj.size());
  
  return(Yj);
}

const sVect<point3d> &
feOsc3d_card::
getH() const
{
  assert(Yj.size() == Hj.size());
  
  return(Hj);
}

const point3d &
feOsc3d_card::
getY(const UInt & i) const
{
  assert(Yj.size() == Hj.size());
  assert(i <= Yj.size());
  assert(i >= 1);
  
  return(Yj(i));
}

const point3d &
feOsc3d_card::
getH(const UInt & i) const
{
  assert(Yj.size() == Hj.size());
  assert(i <= Yj.size());
  assert(i >= 1);
  
  return(Hj(i));
}

UInt
feOsc3d_card::
size() const
{
  assert(Yj.size() == Hj.size());
  
  return(Yj.size());
}
    
sVect<point3d> &
feOsc3d_card::
getY()
{
  assert(Yj.size() == Hj.size());
  
  return(Yj);
}

sVect<point3d> &
feOsc3d_card::
getH()
{
  assert(Yj.size() == Hj.size());
  
  return(Hj);
}

point3d &
feOsc3d_card::
getY(const UInt & i)
{
  assert(Yj.size() == Hj.size());
  assert(i <= Yj.size());
  assert(i >= 1);
  
  return(Yj(i));
}

point3d &
feOsc3d_card::
getH(const UInt & i)
{
  assert(i <= Yj.size());
  assert(i >= 1);
  
  return(Hj(i));
}


//_________________________________________________________________________________________________
// PRINTOUT
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const feOsc3d_card & V)
{
  assert(V.Yj.size() == V.Hj.size());
  
  f << "Osc card " << endl;
  
  for(UInt i=1; i <= V.Yj.size(); ++i)
  { 
    f << "Y" << i << ": " << V.Yj(i);
    f << "H" << i << ": " << V.Hj(i);
  }
  
  return(f);
}
