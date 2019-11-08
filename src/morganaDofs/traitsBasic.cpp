/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "traitsBasic.h"


//_________________________________________________________________________________________________
// REAL TRAIT
//-------------------------------------------------------------------------------------------------
traitsBasic<Real>::
traitsBasic()
{ }

Real
traitsBasic<Real>::
getZero()
{
  return(0.0);
}

Real
traitsBasic<Real>::
getUnity()
{
  return(1.0);
}

Real
traitsBasic<Real>::
getUnityIJ(const UInt & i, const UInt & j)
{
  assert(i == 1);
  assert(j == 1);
  
  return(1.0);
}

Real
traitsBasic<Real>::
getIJ(const UInt & i, const UInt & j, const Real & V)
{
  assert(i == 1);
  assert(j == 1);
  return(V);
}

void
traitsBasic<Real>::
setIJ(const UInt & i, const UInt & j, const Real & value, Real & V)
{
  assert(i == 1);
  assert(j == 1);
  
  V = value;
}

void
traitsBasic<Real>::
setZero(Real & V)
{
  V = 0.0;
}

Real
traitsBasic<Real>::
norm(const Real & A)
{
  return(abs(A));
}



//_________________________________________________________________________________________________
// POINT2D TRAIT
//-------------------------------------------------------------------------------------------------
traitsBasic<point2d>::
traitsBasic()
{ }

point2d
traitsBasic<point2d>::
getZero()
{
  return(point2d(0.0,0.0));
}

point2d
traitsBasic<point2d>::
getUnity()
{
  return(point2d(1.0,1.0));
}

point2d
traitsBasic<point2d>::
getUnityIJ(const UInt & i, const UInt & j)
{
  assert(i >= 1);
  assert(i <= 2);
  assert(j == 1);
  
  point2d P;
  P.setI(i,1.0);
  return(P);
}

Real
traitsBasic<point2d>::
getIJ(const UInt & i, const UInt & j, const point2d & V)
{
  assert(i >= 1);
  assert(i <= 2);
  assert(j == 1);
  
  return(V.getI(i));
}

void
traitsBasic<point2d>::
setIJ(const UInt & i, const UInt & j, const Real & value, point2d & V)
{
  assert(i >= 1);
  assert(i <= 2);
  assert(j == 1);
  
  V.setI(i,value);
}

void 
traitsBasic<point2d>::
setZero(point2d & V)
{
  V.setX(0.0);
  V.setY(0.0);
}

Real 
traitsBasic<point2d>::
norm(const point2d & A)
{
  return(A.norm2());
}



//_________________________________________________________________________________________________
// POINT3D TRAIT
//-------------------------------------------------------------------------------------------------
traitsBasic<point3d>::
traitsBasic()
{ }

point3d
traitsBasic<point3d>::
getZero()
{
  return(point3d(0.0,0.0,0.0));
}

point3d
traitsBasic<point3d>::
getUnity()
{
  return(point3d(1.0,1.0,1.0));
}

point3d
traitsBasic<point3d>::
getUnityIJ(const UInt & i, const UInt & j)
{
  assert(i >= 1);
  assert(i <= 3);
  assert(j == 1);
  
  point3d P;
  P.setI(i,1.0);
  return(P);
}

Real
traitsBasic<point3d>::
getIJ(const UInt & i, const UInt & j, const point3d & V)
{
  assert(i >= 1);
  assert(i <= 3);
  assert(j == 1);
  
  return(V.getI(i));
}

void
traitsBasic<point3d>::
setIJ(const UInt & i, const UInt & j, const Real & value, point3d & V)
{
  assert(i >= 1);
  assert(i <= 3);
  assert(j == 1);
  
  V.setI(i,value);
}

void
traitsBasic<point3d>::
setZero(point3d & V)
{
  V.setX(0.0);
  V.setY(0.0);
  V.setZ(0.0);
}

Real
traitsBasic<point3d>::
norm(const point3d & A)
{
  return(A.norm2());
}



//_________________________________________________________________________________________________
// TENSOR2D TRAIT
//-------------------------------------------------------------------------------------------------
traitsBasic<tensor2d>::
traitsBasic()
{ }

tensor2d
traitsBasic<tensor2d>::
getZero()
{
  return(tensor2d());
}

tensor2d
traitsBasic<tensor2d>::
getUnity()
{
  tensor2d T;
  T.setIJ(1,1,1.0); T.setIJ(1,2,1.0);
  T.setIJ(2,1,1.0); T.setIJ(2,2,1.0);
  
  return(T);
}

tensor2d
traitsBasic<tensor2d>::
getUnityIJ(const UInt & i, const UInt & j)
{
  assert(i > 0);
  assert(i <= 2);
  assert(j > 0);
  assert(j <= 2);
  
  tensor2d T;
  T.setIJ(i,j,1.0);
  
  return(T);
}

Real
traitsBasic<tensor2d>::
getIJ(const UInt & i, const UInt & j, const tensor2d & V)
{
  assert(i >= 1); assert(i <= 2);
  assert(j >= 1); assert(j <= 2);
  
  return(V.getIJ(i,j));
}

void
traitsBasic<tensor2d>::
setIJ(const UInt & i, const UInt & j, const Real & value, tensor2d & V)
{
  assert(i >= 1); assert(i <= 2);
  assert(j >= 1); assert(j <= 2);
  
  V.setIJ(i,j,value);
}

void
traitsBasic<tensor2d>::
setZero(tensor2d & V)
{
  V.setIJ(1,1,0.0); V.setIJ(1,2,0.0);
  V.setIJ(2,1,0.0); V.setIJ(2,2,0.0);
}

Real
traitsBasic<tensor2d>::
norm(const tensor2d & A)
{
  return(A.getFrobenius());
}



//_________________________________________________________________________________________________
// TENSOR3D TRAIT
//-------------------------------------------------------------------------------------------------
traitsBasic<tensor3d>::
traitsBasic()
{ }

tensor3d
traitsBasic<tensor3d>::
getZero()
{
  return(tensor3d());
}

tensor3d
traitsBasic<tensor3d>::
getUnity()
{
  tensor3d T;
  T.setIJ(1,1,1.0); T.setIJ(1,2,1.0); T.setIJ(1,3,1.0);
  T.setIJ(2,1,1.0); T.setIJ(2,2,1.0); T.setIJ(2,3,1.0);
  T.setIJ(3,1,1.0); T.setIJ(3,2,1.0); T.setIJ(3,3,1.0);
  
  return(T);
}

tensor3d
traitsBasic<tensor3d>::
getUnityIJ(const UInt & i, const UInt & j)
{
  assert(i > 0);
  assert(i <= 3);
  assert(j > 0);
  assert(j <= 3);
  
  tensor3d T;
  T.setIJ(i,j,1.0);
  
  return(T);
}

Real
traitsBasic<tensor3d>::
getIJ(const UInt & i, const UInt & j, const tensor3d & V)
{
  assert(i >= 1); assert(i <= 3);
  assert(j >= 1); assert(j <= 3);
  
  return(V.getIJ(i,j));
}

void
traitsBasic<tensor3d>::
setIJ(const UInt & i, const UInt & j, const Real & value, tensor3d & V)
{
  assert(i >= 1); assert(i <= 3);
  assert(j >= 1); assert(j <= 3);
  
  V.setIJ(i,j,value);
}

void
traitsBasic<tensor3d>::
setZero(tensor3d & V)
{
  V.setIJ(1,1,0.0); V.setIJ(1,2,0.0); V.setIJ(1,3,0.0);
  V.setIJ(2,1,0.0); V.setIJ(2,2,0.0); V.setIJ(2,3,0.0);
  V.setIJ(3,1,0.0); V.setIJ(3,2,0.0); V.setIJ(3,3,0.0);
}

Real
traitsBasic<tensor3d>::
norm(const tensor3d & A)
{
  return(A.getFrobenius());
}


//_________________________________________________________________________________________________
// KOMPLEX TRAIT
//-------------------------------------------------------------------------------------------------
traitsBasic<komplex>::
traitsBasic()
{ }

komplex
traitsBasic<komplex>::
getZero()
{
  return(komplex(0.0,0.0));
}

komplex
traitsBasic<komplex>::
getUnity()
{
  return(komplex(1.0,1.0));
}

komplex
traitsBasic<komplex>::
getUnityIJ(const UInt & i, const UInt & j)
{
  assert(i >= 1);
  assert(i <= 2);
  assert(j == 1);
  
  komplex C;
  C.setI(i,1.0);
  return(C);
}

Real
traitsBasic<komplex>::
getIJ(const UInt & i, const UInt & j, const komplex & C)
{
  assert(i >= 1);
  assert(i <= 2);
  assert(j == 1);
  
  return(C.getI(i));
}

void
traitsBasic<komplex>::
setIJ(const UInt & i, const UInt & j, const Real & value, komplex & C)
{
  assert(i >= 1);
  assert(i <= 2);
  assert(j == 1);
  
  C.setI(i,value);
}

void 
traitsBasic<komplex>::
setZero(komplex & C)
{
  C.setReal(0.0);
  C.setImag(0.0);
}

Real 
traitsBasic<komplex>::
norm(const komplex & C)
{
  return(C.norm());
}
