/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef TRAITSBASIC_H
#define TRAITSBASIC_H

#include "point2d.h"
#include "point3d.h"

#include "tensor2d.h"
#include "tensor3d.h"

#include "staticVector.hpp"
#include "komplex.h"


/*! Generic traits - empty */
template<typename TYPE>
class traitsBasic;



//_________________________________________________________________________________________________
// REAL TRAIT
//-------------------------------------------------------------------------------------------------

/*! Real trait */
template<> class traitsBasic<Real>
{
  public:
    static const UInt numI = 1;
    static const UInt numJ = 1;
    static const morganaTypes myType = typeReal;
    
  public:
    traitsBasic();
    static Real getZero();
    static Real getUnity();
    static Real getUnityIJ(const UInt & i, const UInt & j);
    static Real getIJ(const UInt & i, const UInt & j, const Real & V);
    
  public:
    static void setIJ(const UInt & i, const UInt & j, const Real & value, Real & V);
    static void setZero(Real & V);
    
  public:
    static Real norm(const Real & A);
    static size_t memSize(const Real & A);
};



//_________________________________________________________________________________________________
// POINT2D TRAIT
//-------------------------------------------------------------------------------------------------
/*! Point2d trait */
template<> class traitsBasic<point2d>
{
  public:
    static const UInt numI = 2;
    static const UInt numJ = 1;
    static const morganaTypes myType = typePoint2d;
    
  public:
    traitsBasic();
    static point2d getZero();
    static point2d getUnity();
    static point2d getUnityIJ(const UInt & i, const UInt & j);
    static Real    getIJ(const UInt & i, const UInt & j, const point2d & V);
    
  public:
    static void setIJ(const UInt & i, const UInt & j, const Real & value, point2d & V);
    static void setZero(point2d & V);
    
  public:
    static Real norm(const point2d & A);
    static size_t memSize(const point2d & A);
};


//_________________________________________________________________________________________________
// POINT3D TRAIT
//-------------------------------------------------------------------------------------------------
/*! Point3d trait */
template<> class traitsBasic<point3d>
{
  public:
    static const UInt numI = 3;
    static const UInt numJ = 1;
    static const morganaTypes myType = typePoint3d;
    
  public:
    traitsBasic();
    static point3d getZero();
    static point3d getUnity();
    static point3d getUnityIJ(const UInt & i, const UInt & j);
    static Real    getIJ(const UInt & i, const UInt & j, const point3d & V);
    
  public:
    static void setIJ(const UInt & i, const UInt & j, const Real & value, point3d & V);
    static void setZero(point3d & V);
    
  public:
    static Real norm(const point3d & A);
    static size_t memSize(const point3d & A);
};



//_________________________________________________________________________________________________
// TENSOR2D TRAIT
//-------------------------------------------------------------------------------------------------
/*! Tensor2d trait */
template<> class traitsBasic<tensor2d>
{
  public:
    static const UInt numI = 2;
    static const UInt numJ = 2;
    static const morganaTypes myType = typeTensor2d;
    
  public:
    traitsBasic();
    static tensor2d getZero();
    static tensor2d getUnity();
    static tensor2d getUnityIJ(const UInt & i, const UInt & j);
    static Real     getIJ(const UInt & i, const UInt & j, const tensor2d & V);
    
  public:
    static void setIJ(const UInt & i, const UInt & j, const Real & value, tensor2d & V);
    static void setZero(tensor2d & V);
    
  public:
    static Real norm(const tensor2d & A);
    static size_t memSize(const tensor2d & A);
};



//_________________________________________________________________________________________________
// TENSOR3D TRAIT
//-------------------------------------------------------------------------------------------------
/*! Tensor3d trait */
template<> class traitsBasic<tensor3d>
{
  public:
    static const UInt numI = 3;
    static const UInt numJ = 3;
    static const morganaTypes myType = typeTensor3d;
    
  public:
    traitsBasic();
    static tensor3d getZero();
    static tensor3d getUnity();
    static tensor3d getUnityIJ(const UInt & i, const UInt & j);
    static Real     getIJ(const UInt & i, const UInt & j, const tensor3d & V);
    
  public:
    static void setIJ(const UInt & i, const UInt & j, const Real & value, tensor3d & V);
    static void setZero(tensor3d & V);
    
  public:
    static Real norm(const tensor3d & A);
    static size_t memSize(const tensor3d & A);
};



//_________________________________________________________________________________________________
// STATIC VECTOR TRAIT
//-------------------------------------------------------------------------------------------------
/*! StaticVector trait */
template<size_t N> class traitsBasic<staticVector<N> >
{
  public:
    static const UInt numI = N;
    static const UInt numJ = 1;
    static const morganaTypes myType = typeStaticVector;
    
  public:
    traitsBasic();
    static staticVector<N> getZero();
    static staticVector<N> getUnity();
    static staticVector<N> getUnityIJ(const UInt & i, const UInt & j);
    static Real            getIJ(const UInt & i, const UInt & j, const staticVector<N> & V);
    
  public:
    static void setIJ(const UInt & i, const UInt & j, const Real & value, staticVector<N> & V);
    static void setZero(staticVector<N> & V);
    
  public:
    static Real norm(const staticVector<N> & A);
    static size_t memSize(const staticVector<N> & A);
};

template<size_t N>
traitsBasic<staticVector<N> >::
traitsBasic()
{ }

template<size_t N>
staticVector<N>
traitsBasic<staticVector<N> >::
getZero()
{
  staticVector<N> temp;
  
  for(UInt k=1; k <= N; ++k)
  { temp(k) = 0.0; }
  
  return(temp);
}

template<size_t N>
staticVector<N>
traitsBasic<staticVector<N> >::
getUnity()
{
  staticVector<N> temp;
  
  for(UInt k=1; k <= N; ++k)
  { temp(k) = 1.0; }
  
  return(temp);
}

template<size_t N>
staticVector<N>
traitsBasic<staticVector<N> >::
getUnityIJ(const UInt & i, const UInt & j)
{
  staticVector<N> temp;
  
  for(UInt k=1; k <= N; ++k)
  { temp(k) = 0.0; }
  
  assert(i >= 0); assert(i <= N);
  assert(j==1);
  
  temp(i) = 1.0;
  
  return(temp);
}

template<size_t N>
Real
traitsBasic<staticVector<N> >::
getIJ(const UInt & i, const UInt & j, const staticVector<N> & V)
{
  assert(i >= 0); assert(i <= N);
  assert(j==1);
  
  return(V(i));
}

template<size_t N>
void 
traitsBasic<staticVector<N> >::
setIJ(const UInt & i, const UInt & j, const Real & value, staticVector<N> & V)
{
  assert(i <= N);
  assert(j==1);
  
  V(i) = value;
}

template<size_t N>
void 
traitsBasic<staticVector<N> >::
setZero(staticVector<N> & V)
{
  for(UInt k=1; k <= N; ++k)
  { V(k) = 0.0; }
}

template<size_t N>
Real 
traitsBasic<staticVector<N> >::
norm(const staticVector<N> & A)
{
  Real tot = 0.0;
  
  for(UInt k=1; k <= N; ++k)
  { tot += A(k); }
  
  return(tot);
}

template<size_t N>
size_t
traitsBasic<staticVector<N> >::
memSize(const staticVector<N> & A)
{
  return(A.memSize());
}



//_________________________________________________________________________________________________
// KOMPLEX TRAIT
//-------------------------------------------------------------------------------------------------
/*! Komplex trait */
template<> class traitsBasic<komplex>
{
  public:
    static const UInt numI = 2;
    static const UInt numJ = 1;
    static const morganaTypes myType = typeKomplex;
    
  public:
    traitsBasic();
    static komplex getZero();
    static komplex getUnity();
    static komplex getUnityIJ(const UInt & i, const UInt & j);
    static Real    getIJ(const UInt & i, const UInt & j, const komplex & C);
    
  public:
    static void setIJ(const UInt & i, const UInt & j, const Real & value, komplex & C);
    static void setZero(komplex & C);
    
  public:
    static Real norm(const komplex & C);
    static size_t memSize(const komplex & A);
};


#endif
