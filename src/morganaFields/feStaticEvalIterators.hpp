/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FESTATICEVALITERATORS_HPP
#define FESTATICEVALITERATORS_HPP

#include "morganaTypes.hpp"
#include "polyStatic.hpp"
#include "morganaFields.hpp"
#include "feStaticLists.hpp"


// REAL EVALUATION CLASS___________________________________________________________________________
template<FELabel T, UInt I>
struct feStaticRealEval
{
  typedef typename feStaticList<T,I>::BASE POLYCARD;
  
  static void eval(sVect<Real> & list, const point3d & Y)
  {
    list(I) = polyStatic<POLYCARD>::evaluateStatic(Y);
    feStaticRealEval<T,I-1>::eval(list,Y);
  }
};

template<FELabel T>
struct feStaticRealEval<T, 1>
{
  typedef typename feStaticList<T,1>::BASE POLYCARD;
  
  static void eval(sVect<Real> & list, const point3d & Y)
  {
    list(1) = polyStatic<POLYCARD>::evaluateStatic(Y);
  }
};



// REAL DX CLASS___________________________________________________________________________________
template<FELabel T, UInt I>
struct feStaticRealDx
{
  typedef typename feStaticList<T,I>::BASE POLYCARD;
  
  static void eval(sVect<Real> & list, const point3d & Y)
  {
    list(I) = polyStatic<POLYCARD>::evaluateGradientX(Y);
    feStaticRealDx<T,I-1>::eval(list,Y);
  }
};

template<FELabel T>
struct feStaticRealDx<T, 1>
{
  typedef typename feStaticList<T,1>::BASE POLYCARD;
  
  static void eval(sVect<Real> & list, const point3d & Y)
  {
    list(1) = polyStatic<POLYCARD>::evaluateGradientX(Y);
  }
};



// REAL DY CLASS___________________________________________________________________________________
template<FELabel T, UInt I>
struct feStaticRealDy
{
  typedef typename feStaticList<T,I>::BASE POLYCARD;
  
  static void eval(sVect<Real> & list, const point3d & Y)
  {
    list(I) = polyStatic<POLYCARD>::evaluateGradientY(Y);
    feStaticRealDy<T,I-1>::eval(list,Y);
  }
};

template<FELabel T>
struct feStaticRealDy<T, 1>
{
  typedef typename feStaticList<T,1>::BASE POLYCARD;
  
  static void eval(sVect<Real> & list, const point3d & Y)
  {
    list(1) = polyStatic<POLYCARD>::evaluateGradientY(Y);
  }
};



// REAL DZ CLASS___________________________________________________________________________________
template<FELabel T, UInt I>
struct feStaticRealDz
{
  typedef typename feStaticList<T,I>::BASE POLYCARD;
  
  static void eval(sVect<Real> & list, const point3d & Y)
  {
    list(I) = polyStatic<POLYCARD>::evaluateGradientZ(Y);
    feStaticRealDz<T,I-1>::eval(list,Y);
  }
};

template<FELabel T>
struct feStaticRealDz<T, 1>
{
  typedef typename feStaticList<T,1>::BASE POLYCARD;
  
  static void eval(sVect<Real> & list, const point3d & Y)
  {
    list(1) = polyStatic<POLYCARD>::evaluateGradientZ(Y);
  }
};



// REAL TEMPLATE DERIVATIVE CLASS__________________________________________________________________
template<FELabel T, UInt I, UInt dx, UInt dy, UInt dz>
struct feStaticRealDerivative
{
  typedef typename feStaticList<T,I>::BASE POLYCARD;
  
  static void eval(sVect<Real> & list, const point3d & Y)
  {
    list(I) = polyStatic<POLYCARD>::template evaluateDerivative <dx,dy,dz> (Y);
    feStaticRealDerivative<T,I-1, dx,dy,dz>::eval(list,Y);
  }
};

template<FELabel T, UInt dx, UInt dy, UInt dz>
struct feStaticRealDerivative<T, 1, dx, dy, dz>
{
  typedef typename feStaticList<T,1>::BASE POLYCARD;
  
  static void eval(sVect<Real> & list, const point3d & Y)
  {
    list(1) = polyStatic<POLYCARD>::template evaluateDerivative <dx,dy,dz> (Y);
  }
};


// POINT3D EVALUATION CLASS________________________________________________________________________
template<FELabel T, UInt I>
struct feStaticPoint3dEval
{
  typedef typename feStaticListX<T,I>::BASE POLYCARD_X;
  typedef typename feStaticListY<T,I>::BASE POLYCARD_Y;
  typedef typename feStaticListZ<T,I>::BASE POLYCARD_Z;
  
  static void eval(sVect<point3d> & list, const point3d & Y)
  {
    list(I).setX( polyStatic<POLYCARD_X>::evaluateStatic(Y) );
    list(I).setY( polyStatic<POLYCARD_Y>::evaluateStatic(Y) );
    list(I).setZ( polyStatic<POLYCARD_Z>::evaluateStatic(Y) );

    feStaticPoint3dEval<T,I-1>::eval(list,Y);
  }
};

template<FELabel T>
struct feStaticPoint3dEval<T,1>
{
  typedef typename feStaticListX<T,1>::BASE POLYCARD_X;
  typedef typename feStaticListY<T,1>::BASE POLYCARD_Y;
  typedef typename feStaticListZ<T,1>::BASE POLYCARD_Z;
  
  static void eval(sVect<point3d> & list, const point3d & Y)
  {
    list(1).setX( polyStatic<POLYCARD_X>::evaluateStatic(Y) );
    list(1).setY( polyStatic<POLYCARD_Y>::evaluateStatic(Y) );
    list(1).setZ( polyStatic<POLYCARD_Z>::evaluateStatic(Y) );
  }
};


// POINT 3D DX CLASS_______________________________________________________________________________
template<FELabel T, UInt I>
struct feStaticPoint3dDx
{
  typedef typename feStaticListX<T,I>::BASE POLYCARD_X;
  typedef typename feStaticListY<T,I>::BASE POLYCARD_Y;
  typedef typename feStaticListZ<T,I>::BASE POLYCARD_Z;
  
  static void eval(sVect<point3d> & list, const point3d & Y)
  {
    list(I).setX( polyStatic<POLYCARD_X>::evaluateGradientX(Y) );
    list(I).setY( polyStatic<POLYCARD_Y>::evaluateGradientX(Y) );
    list(I).setZ( polyStatic<POLYCARD_Z>::evaluateGradientX(Y) );
    
    feStaticPoint3dDx<T,I-1>::eval(list,Y);
  }
};

template<FELabel T>
struct feStaticPoint3dDx<T, 1>
{
  typedef typename feStaticListX<T,1>::BASE POLYCARD_X;
  typedef typename feStaticListY<T,1>::BASE POLYCARD_Y;
  typedef typename feStaticListZ<T,1>::BASE POLYCARD_Z;
  
  static void eval(sVect<point3d> & list, const point3d & Y)
  {
    list(1).setX( polyStatic<POLYCARD_X>::evaluateGradientX(Y) );
    list(1).setY( polyStatic<POLYCARD_Y>::evaluateGradientX(Y) );
    list(1).setZ( polyStatic<POLYCARD_Z>::evaluateGradientX(Y) );
  }
};



// POINT 3D DY CLASS_______________________________________________________________________________
template<FELabel T, UInt I>
struct feStaticPoint3dDy
{
  typedef typename feStaticListX<T,I>::BASE POLYCARD_X;
  typedef typename feStaticListY<T,I>::BASE POLYCARD_Y;
  typedef typename feStaticListZ<T,I>::BASE POLYCARD_Z;
  
  static void eval(sVect<point3d> & list, const point3d & Y)
  {
    list(I).setX( polyStatic<POLYCARD_X>::evaluateGradientY(Y) );
    list(I).setY( polyStatic<POLYCARD_Y>::evaluateGradientY(Y) );
    list(I).setZ( polyStatic<POLYCARD_Z>::evaluateGradientY(Y) );
    
    feStaticPoint3dDy<T,I-1>::eval(list,Y);
  }
};

template<FELabel T>
struct feStaticPoint3dDy<T, 1>
{
  typedef typename feStaticListX<T,1>::BASE POLYCARD_X;
  typedef typename feStaticListY<T,1>::BASE POLYCARD_Y;
  typedef typename feStaticListZ<T,1>::BASE POLYCARD_Z;
  
  static void eval(sVect<point3d> & list, const point3d & Y)
  {
    list(1).setX( polyStatic<POLYCARD_X>::evaluateGradientY(Y) );
    list(1).setY( polyStatic<POLYCARD_Y>::evaluateGradientY(Y) );
    list(1).setZ( polyStatic<POLYCARD_Z>::evaluateGradientY(Y) );
  }
};



// POINT 3D DZ CLASS_______________________________________________________________________________
template<FELabel T, UInt I>
struct feStaticPoint3dDz
{
  typedef typename feStaticListX<T,I>::BASE POLYCARD_X;
  typedef typename feStaticListY<T,I>::BASE POLYCARD_Y;
  typedef typename feStaticListZ<T,I>::BASE POLYCARD_Z;
  
  static void eval(sVect<point3d> & list, const point3d & Y)
  {
    list(I).setX( polyStatic<POLYCARD_X>::evaluateGradientZ(Y) );
    list(I).setY( polyStatic<POLYCARD_Y>::evaluateGradientZ(Y) );
    list(I).setZ( polyStatic<POLYCARD_Z>::evaluateGradientZ(Y) );
    
    feStaticPoint3dDz<T,I-1>::eval(list,Y);
  }
};

template<FELabel T>
struct feStaticPoint3dDz<T, 1>
{
  typedef typename feStaticListX<T,1>::BASE POLYCARD_X;
  typedef typename feStaticListY<T,1>::BASE POLYCARD_Y;
  typedef typename feStaticListZ<T,1>::BASE POLYCARD_Z;
  
  static void eval(sVect<point3d> & list, const point3d & Y)
  {
    list(1).setX( polyStatic<POLYCARD_X>::evaluateGradientZ(Y) );
    list(1).setY( polyStatic<POLYCARD_Y>::evaluateGradientZ(Y) );
    list(1).setZ( polyStatic<POLYCARD_Z>::evaluateGradientZ(Y) );
  }
};



// POINT3D TEMPLATE DERIVATIVE CLASS_______________________________________________________________
template<FELabel T, UInt I, UInt dx, UInt dy, UInt dz>
struct feStaticPoint3dDerivative
{
  typedef typename feStaticListX<T,I>::BASE POLYCARD_X;
  typedef typename feStaticListY<T,I>::BASE POLYCARD_Y;
  typedef typename feStaticListZ<T,I>::BASE POLYCARD_Z;
  
  static void eval(sVect<point3d> & list, const point3d & Y)
  {
    list(I).setX( polyStatic<POLYCARD_X>::template evaluateDerivative <dx,dy,dz> (Y) );
    list(I).setY( polyStatic<POLYCARD_Y>::template evaluateDerivative <dx,dy,dz> (Y) );
    list(I).setZ( polyStatic<POLYCARD_Z>::template evaluateDerivative <dx,dy,dz> (Y) ); 
    
    feStaticPoint3dDerivative<T,I-1, dx,dy,dz>::eval(list,Y);
  }
};

template<FELabel T, UInt dx, UInt dy, UInt dz>
struct feStaticPoint3dDerivative<T,1,dx,dy,dz>
{
  typedef typename feStaticListX<T,1>::BASE POLYCARD_X;
  typedef typename feStaticListY<T,1>::BASE POLYCARD_Y;
  typedef typename feStaticListZ<T,1>::BASE POLYCARD_Z;
  
  static void eval(sVect<point3d> & list, const point3d & Y)
  {
    list(1).setX( polyStatic<POLYCARD_X>::template evaluateDerivative <dx,dy,dz> (Y) );
    list(1).setY( polyStatic<POLYCARD_Y>::template evaluateDerivative <dx,dy,dz> (Y) );
    list(1).setZ( polyStatic<POLYCARD_Z>::template evaluateDerivative <dx,dy,dz> (Y) ); 
  }
};



#endif
