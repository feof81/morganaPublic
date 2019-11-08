/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "stateMatrix.h"

//_________________________________________________________________________________________________
// CONSTRUCTORS FUNCTIONS
//-------------------------------------------------------------------------------------------------
stateMatrix::
stateMatrix() : Epetra_SerialDenseMatrix()
{ }

stateMatrix::
stateMatrix(const UInt & NRows, const UInt & NCols) : Epetra_SerialDenseMatrix(NRows,NCols)
{ }

stateMatrix::
stateMatrix(const stateMatrix &M) : Epetra_SerialDenseMatrix(M)
{ }



//_________________________________________________________________________________________________
// OPERATORS
//-------------------------------------------------------------------------------------------------
void
stateMatrix::
operator+=(const stateMatrix & M)
{
  if( (M.RowDim() != this->RowDim()) || (M.ColDim() != this->ColDim()) )
  { 
    this->Reshape(M.RowDim(),M.ColDim());
  }
  
  for(int i=1; i<=M.RowDim(); ++i)
  {
    for(int j=1; j<=M.ColDim(); ++j)
    {
      this->operator()(i,j) += M(i,j);
    }
  }
}

void 
stateMatrix::
operator-=(const stateMatrix & M)
{
  if( (M.RowDim() != this->RowDim()) || (M.ColDim() != this->ColDim()) )
  { 
    this->Reshape(M.RowDim(),M.ColDim());
  }
  
  for(int i=1; i<=M.RowDim(); ++i)
  {
    for(int j=1; j<=M.ColDim(); ++j)
    {
      this->operator()(i,j) -= M(i,j);
    }
  }
}

void
stateMatrix::
operator*=(const Real &a)
{
  for(int i=1; i<=this->RowDim(); ++i)
  {
    for(int j=1; j<=this->ColDim(); ++j)
    {
      this->operator()(i,j) = this->operator()(i,j) * a;
    }
  }
}

void
stateMatrix::
operator/=(const Real &a)
{
  for(int i=1; i<=this->RowDim(); ++i)
  {
    for(int j=1; j<=this->ColDim(); ++j)
    {
      this->operator()(i,j) = this->operator()(i,j) / a;
    }
  }
}


//_________________________________________________________________________________________________
// GET - SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
stateVector
stateMatrix::
getRow(const UInt & i)
{
  assert((int)i <= this->RowDim());
  stateVector out(this->ColDim());
  
  for(int j=1; j<=this->ColDim(); ++j)
  {
    out(j) = this->operator()(i,j);
  }
  
  return(out);
}

stateVector
stateMatrix::
getCol(const UInt & j)
{
  assert((int)j <= this->ColDim());
  stateVector out(this->RowDim());
  
  for(int i=1; i<=this->RowDim(); ++i)
  {
    out(i) = this->operator()(i,j);
  }
  
  return(out);
}

void
stateMatrix::
setRow(const UInt & i, const stateVector & V)
{
  assert((int)i <= this->RowDim());
  assert(V.Length() == this->ColDim());
  
  for(int j=1; j<=this->ColDim(); ++j)
  {
    this->operator()(i,j) = V(j);
  }
}

void
stateMatrix::
setCol(const UInt & j, const stateVector & V)
{
  assert((int)j <= this->ColDim());
  assert(V.Length() == this->RowDim());
  
  for(int i=1; i<=this->RowDim(); ++i)
  {
    this->operator()(i,j) = V(i);
  }
}



//_________________________________________________________________________________________________
// ORDINAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
bool
stateMatrix::
operator<(const stateMatrix & V) const
{
  assert(V.RowDim() != this->RowDim());
  assert(V.ColDim() != this->ColDim());
  
  
  for(int i=1; i <= this->RowDim(); ++i)
  {
    for(int j=1; j <= this->ColDim(); ++j)
    {
      if( this->operator()(i,j) < (V.operator()(i,j) - geoToll) )
      {
	return(true);
      }
      
      if( this->operator()(i,j) > (V.operator()(i,j) + geoToll) )
      {
	return(false);
      }
    }
  }
  
  return(false);
}

bool 
stateMatrix::
operator!=(const stateMatrix & V) const
{
  assert(V.RowDim() != this->RowDim());
  assert(V.ColDim() != this->ColDim());
  
  
  for(int i=1; i <= this->RowDim(); ++i)
  {
    for(int j=1; j <= this->ColDim(); ++j)
    {
      if( ( this->operator()(i,j) > (V.operator()(i,j) + geoToll) ) || ( this->operator()(i,j) < (V.operator()(i,j) - geoToll) ) )
      {
	return(true);
      }
    }
  }
  
  return(false);
}



//_________________________________________________________________________________________________
// OTHER FUNCTIONS
//-------------------------------------------------------------------------------------------------
ostream & 
operator<<( ostream & f, const stateMatrix & A)
{
  for(int i=1; i <= A.RowDim(); ++i)
  { 
    for(int j=1; j <= A.ColDim(); ++j)
    {
      f << A(i,j) << " ";
    }
    f << endl;
  }
  
  f << endl;
  
  return(f);
}
