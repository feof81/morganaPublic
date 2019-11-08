/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "traitsMultiply.h"


//_________________________________________________________________________________________________
// REAL-REAL TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<Real,Real>::DIRECTTYPE
traitsMultiply<Real,Real>::
multiply(const Real & A, const Real & B)
{
  return(A*B);
}

traitsMultiply<Real,Real>::SCALARTYPEA
traitsMultiply<Real,Real>::
scalarProductA(const Real & A, const Real & B)
{
  return(A*B);
}


//_________________________________________________________________________________________________
// REAL-POINT2D TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<Real,point2d>::DIRECTTYPE
traitsMultiply<Real,point2d>::
multiply(const Real & A, const point2d & B)
{
  return(B*A);
}


//_________________________________________________________________________________________________
// REAL-POINT3D TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<Real,point3d>::DIRECTTYPE
traitsMultiply<Real,point3d>::
multiply(const Real & A, const point3d & B)
{
  return(B*A);
}


//_________________________________________________________________________________________________
// REAL-TENSOR2D TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<Real,tensor2d>::DIRECTTYPE
traitsMultiply<Real,tensor2d>::
multiply(const Real & A, const tensor2d & B)
{
  return(B*A);
}


//_________________________________________________________________________________________________
// REAL-TENSOR3D TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<Real,tensor3d>::DIRECTTYPE
traitsMultiply<Real,tensor3d>::
multiply(const Real & A, const tensor3d & B)
{
  return(B*A);
}


//_________________________________________________________________________________________________
// REAL-STATEVECTOR TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<Real,stateVector>::DIRECTTYPE
traitsMultiply<Real,stateVector>::
multiply(const Real & A, const stateVector & B)
{
  return(B*A);
}


//_________________________________________________________________________________________________
// REAL-STATEMATRIX TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<Real,stateMatrix>::DIRECTTYPE
traitsMultiply<Real,stateMatrix>::
multiply(const Real & A, const stateMatrix & B)
{
  return(B*A);
}


//_________________________________________________________________________________________________
// REAL-POINT2D TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<Real,komplex>::DIRECTTYPE
traitsMultiply<Real,komplex>::
multiply(const Real & A, const komplex & B)
{
  return(B*A);
}


//_________________________________________________________________________________________________
// POINT2D-REAL TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<point2d,Real>::DIRECTTYPE 
traitsMultiply<point2d,Real>::
multiply(const point2d & A, const Real & B)
{
  return(A*B);
}


//_________________________________________________________________________________________________
// POINT2D-POINT2D TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<point2d,point2d>::DIRECTTYPE 
traitsMultiply<point2d,point2d>::
multiply(const point2d & A, const point2d & B)
{
  tensor2d T;
  T.setIJ(1,1,A.getX()*B.getX()); T.setIJ(1,2,A.getX()*B.getY());
  T.setIJ(2,1,A.getY()*B.getX()); T.setIJ(2,2,A.getY()*B.getY());
  
  return(T);
}

traitsMultiply<point2d,point2d>::CROSSTYPEA 
traitsMultiply<point2d,point2d>::
crossProductA(const point2d & A, const point2d & B)
{
  return(point3d(0.0,0.0, A.getX()*B.getY() - A.getY()*B.getX() ));
}

traitsMultiply<point2d,point2d>::SCALARTYPEA
traitsMultiply<point2d,point2d>::
scalarProductA(const point2d & A, const point2d & B)
{
  return(A*B);
}


//_________________________________________________________________________________________________
// POINT3D-REAL TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<point3d,Real>::DIRECTTYPE 
traitsMultiply<point3d,Real>::
multiply(const point3d & A, const Real & B)
{
  return(A*B);
}


//_________________________________________________________________________________________________
// POINT3D-POINT3D TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<point3d,point3d>::DIRECTTYPE 
traitsMultiply<point3d,point3d>::
multiply(const point3d & A, const point3d & B)
{
  tensor3d T;
  T.setIJ(1,1,A.getX() * B.getX()); T.setIJ(1,2,A.getX() * B.getY()); T.setIJ(1,3,A.getX() * B.getZ());
  T.setIJ(2,1,A.getY() * B.getX()); T.setIJ(2,2,A.getY() * B.getY()); T.setIJ(2,3,A.getY() * B.getZ());
  T.setIJ(3,1,A.getZ() * B.getX()); T.setIJ(3,2,A.getZ() * B.getY()); T.setIJ(3,3,A.getZ() * B.getZ());
  
  return(T);
}

traitsMultiply<point3d,point3d>::CROSSTYPEA
traitsMultiply<point3d,point3d>::
crossProductA(const point3d & A, const point3d & B)
{
  return(A^B);
}

traitsMultiply<point3d,point3d>::SCALARTYPEA
traitsMultiply<point3d,point3d>::
scalarProductA(const point3d & A, const point3d & B)
{
  return(A*B);
}


//_________________________________________________________________________________________________
// TENSOR2D-REAL TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<tensor2d,Real>::DIRECTTYPE
traitsMultiply<tensor2d,Real>::
multiply(const tensor2d & A, const Real & B)
{
  return(A*B);
}


//_________________________________________________________________________________________________
// TENSOR3D-REAL TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<tensor3d,Real>::DIRECTTYPE
traitsMultiply<tensor3d,Real>::
multiply(const tensor3d & A, const Real & B)
{
  return(A*B);
}


//_________________________________________________________________________________________________
// STATEVECTOR-REAL TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<stateVector,Real>::DIRECTTYPE
traitsMultiply<stateVector,Real>::
multiply(const stateVector & A, const Real & B)
{
  return(A*B);
}


//_________________________________________________________________________________________________
// STATEMATRIX-REAL TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<stateMatrix,Real>::DIRECTTYPE
traitsMultiply<stateMatrix,Real>::
multiply(const stateMatrix & A, const Real & B)
{
  return(A*B);
}


//_________________________________________________________________________________________________
// KOMPLEX-REAL TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<komplex,Real>::DIRECTTYPE 
traitsMultiply<komplex,Real>::
multiply(const komplex & A, const Real & B)
{
  return(A*B);
}


//_________________________________________________________________________________________________
// KOMPLEX-REAL TRAIT
//-------------------------------------------------------------------------------------------------
traitsMultiply<komplex,komplex>::DIRECTTYPE 
traitsMultiply<komplex,komplex>::
multiply(const komplex & A, const komplex & B)
{
  komplex C = A * B;
  
  return(C);
}
