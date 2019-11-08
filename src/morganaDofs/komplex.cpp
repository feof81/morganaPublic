/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "komplex.h"
#include <cmath>


//_________________________________________________________________________________________________
// CONSTRUCTORS AND DESTRUCTORS
//-------------------------------------------------------------------------------------------------
komplex::
komplex(Real Re, Real Im)
{
  re = Re;
  im = Im;
}

komplex::
komplex(const komplex & V)
{
  re = V.re;
  im = V.im;
}

komplex::
~komplex()
{
}

//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
komplex::
set(const Real & Re, const Real & Im)
{
  re = Re;
  im = Im;
}
    
void
komplex::
setReal(const Real & Re)
{
  re = Re;
}
    
void
komplex::
setImag(const Real & Im)
{
  im = Im;
}

void
komplex::
setI(const UInt & i, const Real & val)
{
  assert(i <= 2);
  assert(i >= 1);
  
  if(i == 1) { re = val; }
  else       { im = val; }
}


//_________________________________________________________________________________________________
// OPERATORS FUNCTIONS
//-------------------------------------------------------------------------------------------------
Real
komplex::
operator()(const UInt & i)
{
  assert(i <= 2);
  assert(i >= 1);
  
  if(i == 1) { return(re); }
  else       { return(im); }
}    

const Real &
komplex::
operator()(const UInt & i) const
{
  assert(i <= 2);
  assert(i >= 1);
  
  if(i == 1) { return(re); }
  else       { return(im); }
}   

komplex &
komplex::
operator=(const komplex & C)
{
  re = C.re;
  im = C.im;
  
  return *this;
}

void
komplex::
operator+=(const komplex & C)
{
  re += C.re;
  im += C.im;
}    

void
komplex::
operator-=(const komplex & C)
{
  re -= C.re;
  im -= C.im;
}  

void
komplex::
operator*=(const Real & a)
{
  re *= a;
  im *= a;
}    

void
komplex::
operator/=(const Real & a)
{
  re /= a;
  im /= a;
}    

void
komplex::
operator*=(const komplex & C)
{
  Real a = re;
  Real b = im;
  
  re = a * C.re - b * C.im;
  im = b * C.re + a * C.im;
}

void
komplex::
operator/=(const komplex & C)
{
  Real rhoA = sqrt(  re * re   +   im * im);
  Real rhoB = sqrt(C.re * C.re + C.im + C.im);
  
  Real thetaA = atan2(im,re);
  Real thetaB = atan2(C.im,C.re);
  
  re = (rhoA / rhoB) * cos(thetaA - thetaB);
  im = (rhoA / rhoB) * sin(thetaA - thetaB);
}


//_________________________________________________________________________________________________
// ORDINAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
bool
komplex::
operator<(const komplex & C) const
{
  if(re < (C.re - geoToll)) { return(true); }
  if(re > (C.re + geoToll)) { return(false); }
  
  if(im < (C.im - geoToll)) { return(true); }
  if(im > (C.im + geoToll)) { return(false); }
  
  return(false);
}

bool
komplex::
operator!=(const komplex & C) const
{
  bool equalRe, equalIm;

  equalRe = (re >= (C.re - geoToll)) && (re <= (C.re + geoToll));
  equalIm = (im >= (C.im - geoToll)) && (im <= (C.im + geoToll));

  return(!(equalRe && equalIm) );
}
    
bool
komplex::
operator==(const komplex & C) const
{
  return( !this->operator!=(C) );
}


//_________________________________________________________________________________________________
// OTHER FUNCTIONS
//-------------------------------------------------------------------------------------------------
Real
komplex::
norm() const
{
  return(sqrt(re * re + im * im));
}

Real
komplex::
phase() const
{
  return(atan2(im,re));
}
    
Real
komplex::
norm(const komplex & C)
{
  return(sqrt(C.re * C.re + C.im * C.im));
}

Real
komplex::
phase(const komplex & C)
{
  return(atan2(C.im,C.re));
}
    
Real
komplex::
cReal(const komplex & C)
{
  return(C.re);
}

Real
komplex::
cReal(const Real & C)
{
  return(C);
}
    
Real
komplex::
cImag(const komplex & C)
{
  return(C.im);
}

Real
komplex::
cImag(const Real & C)
{
  return(0.0);
}

Real
komplex::
cComp(const komplex & C, const UInt & i)
{
  assert(i <= 2);
  assert(i >= 1);
  
  if(i == 1) { return(C.re); }
  else       { return(C.im); }
}

Real
komplex::
cComp(const Real & C, const UInt & i)
{
  assert(i <= 2);
  assert(i >= 1);
  
  if(i == 1) { return(C); }
  else       { return(0.0); }
}

komplex
komplex::
iexp(const Real & theta)
{
  return(komplex(cos(theta),sin(theta)));
}

komplex
komplex::
conj(const komplex & C)
{
  return(komplex(C.re,-C.im));
}

ostream &
operator<<(ostream& f, const komplex & C)
{
  return f << "Real: " << C.re << " Imag: " << C.im << endl;
}

 