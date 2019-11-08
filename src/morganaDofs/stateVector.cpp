/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "stateVector.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
stateVector::
stateVector(const UInt & N) : Epetra_SerialDenseVector(N)
{ }

stateVector::
stateVector(Real *Comp, const UInt & N) : Epetra_SerialDenseVector(Copy,Comp,N)
{ }

stateVector::
stateVector(const stateVector &C) : Epetra_SerialDenseVector(C)
{ }



//_________________________________________________________________________________________________
// OPERATORS
//-------------------------------------------------------------------------------------------------
void 
stateVector::
operator+=(const stateVector &C)
{
  if(C.Length() != this->Length())
  { this->Resize(C.Length()); }
  
  for(int i=1; i <= this->Length(); ++i)
  { this->operator()(i) += C(i); }
}

void
stateVector::
operator-=(const stateVector & C)
{
  if(C.Length() != this->Length())
  { this->Resize(C.Length()); }
  
  for(int i=1; i <= this->Length(); ++i)
  { this->operator()(i) -= C(i); }
}

void
stateVector::
operator*=(const Real &a)
{
  for(int i=1; i <= this->Length(); ++i)
  { this->operator()(i) *= a; }
}

void 
stateVector::
operator/=(const Real &a)
{
  for(int i=1; i <= this->Length(); ++i)
  { this->operator()(i) /= a; }
}

Real
stateVector::
operator*(const stateVector &C) const
{
  assert(C.Length() == this->Length());
  Real tot=0.0;
  
  for(int i=1; i<=this->Length(); ++i)
  { tot += this->operator()(i) * C(i); }
  
  return(tot);
}



//_________________________________________________________________________________________________
// INTERNAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
Real 
stateVector::
sum() const
{
  Real tot = 0.0;
  
  for(int i=1; i<=this->Length(); ++i)
  { tot += this->operator()(i); }
  
  return(tot);
}

Real
stateVector::
norm1() const
{
   Real tot = 0.0;
  
  for(int i=1; i<=this->Length(); ++i)
  { tot += abs(this->operator()(i)); }
  
  return(tot);
}

void
stateVector::
cancel()
{
  for(int i=1; i <= this->Length(); ++i)
  { this->operator()(i) = 0.0; }
}

UInt
stateVector::
size() const
{
  return(UInt(this->Length()));
}


//_________________________________________________________________________________________________
// ORDINAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
bool
stateVector::
operator<(const stateVector & V) const
{
  assert(V.Length() == this->Length());
  
  for(int i=1; i <= V.Length(); ++i)
  {
    if( this->operator()(i) < (V.operator()(i) - geoToll) )
    {
      return(true);
    }
      
    if( this->operator()(i) > (V.operator()(i) - geoToll) )
    {
      return(false);
    }
  }
  
  return(false);
}

bool 
stateVector::
operator!=(const stateVector & V) const
{
  assert(V.Length() == this->Length());
  
  for(int i=1; i <= V.Length(); ++i)
  {
    if( ( this->operator()(i) > (V.operator()(i) + geoToll) ) || ( this->operator()(i) < (V.operator()(i) - geoToll) ) )
    {
      return(true);
    }
  }
  
  return(false);
}



//_________________________________________________________________________________________________
// OTHER FUNCTIONS
//-------------------------------------------------------------------------------------------------
ostream & 
operator<<( ostream & f, const stateVector & A)
{
  for(int i=1; i <= A.Length(); ++i)
  { f << A(i) << " "; }
  
  f << endl;
  
  return(f);
}
