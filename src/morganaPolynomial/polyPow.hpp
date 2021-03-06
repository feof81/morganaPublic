/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef POLYPOW_HPP
#define POLYPOW_HPP

#include <assert.h>
#include <iostream>
#include <fstream>

#include "typesInterface.hpp"
#include "polyCards.h"


//FACTORIAL CLASS----------------------------------------------------------------------------------
template<UInt D>
struct polyFactorial
{
  static Real eval(const UInt & N)
  {return(polyFactorial<D-1>::eval(N) * (N - D + 1) ); }
};

template<>
struct polyFactorial<0>
{
  static Real eval(const UInt & N)
  {
    assert(N == N);
    return(1.0);
  }
};



//STATIC POWER CLASS-------------------------------------------------------------------------------
template<UInt N> 
struct polyPow
{
    static Real eval(const Real & x)
    { return(polyPow<N-1>::eval(x) * x); }
};

template<> 
struct polyPow<0>
{
    static Real eval(const Real & x)
    {
      assert(x == x);
      return(1.0);
    }
};


//SWITCH CLASS-------------------------------------------------------------------------------------
template<bool,Int N> class polySwitch;

template<Int N> 
struct polySwitch<true,N>
{
  static Real eval(const Real & x)
  { return(polyPow<N>::eval(x)); }
};

template<Int N> 
struct polySwitch<false,N>
{
  static Real eval(const Real & x)
  {
    assert(x == x);
    return(0.0);
  }
};



//STATIC POSITIVE POWER CLASS----------------------------------------------------------------------
template<Int N> 
struct polyPosPow
{
  static Real eval(const Real & x)
  { return(polySwitch< N>=0 , N>::eval(x)); }
};



//STATIC ITERATOR FOR EVALUATION-------------------------------------------------------------------
template<typename POLYCARD, UInt S> class polyEvalIter;

template<typename POLYCARD, UInt S>
struct polyEvalIter
{
  typedef typename POLYCARD::SON SON; 
  
  static Real eval(const point3d & P)
  {
    return(polyEvalIter<SON,S-1>::eval(P) +
    POLYCARD::C *
    polyPosPow<POLYCARD::Px>::eval(P.getX()) *
    polyPosPow<POLYCARD::Py>::eval(P.getY()) *
    polyPosPow<POLYCARD::Pz>::eval(P.getZ()) );
  }
};

template<typename POLYCARD>
struct polyEvalIter<POLYCARD,1>
{
  static Real eval(const point3d & P)
  {
    return(POLYCARD::C *
    polyPosPow<POLYCARD::Px>::eval(P.getX()) *
    polyPosPow<POLYCARD::Py>::eval(P.getY()) *
    polyPosPow<POLYCARD::Pz>::eval(P.getZ()) );
  }
};



//STATIC ITERATOR FOR DERIVATIVE EVALUATION--------------------------------------------------------
template<typename POLYCARD, UInt dx, UInt dy, UInt dz, UInt S> class polyDerIter;

template<typename POLYCARD, UInt dx, UInt dy, UInt dz, UInt S>
struct polyDerIter
{
  typedef typename POLYCARD::SON SON;
  
  static Real eval(const point3d & P)
  {
    return(polyDerIter<SON,dx,dy,dz,S-1>::eval(P) +
    POLYCARD::C *
    polyFactorial<dx>::eval(UInt(POLYCARD::Px)) * polyPosPow<UInt(POLYCARD::Px - dx)>::eval(P.getX()) *
    polyFactorial<dy>::eval(UInt(POLYCARD::Py)) * polyPosPow<UInt(POLYCARD::Py - dy)>::eval(P.getY()) *
    polyFactorial<dz>::eval(UInt(POLYCARD::Pz)) * polyPosPow<UInt(POLYCARD::Pz - dz)>::eval(P.getZ()) );
  }
};

template<typename POLYCARD, UInt dx, UInt dy, UInt dz>
struct polyDerIter<POLYCARD,dx,dy,dz,1>
{
  static Real eval(const point3d & P)
  {   
    return(POLYCARD::C *
    polyFactorial<dx>::eval(UInt(POLYCARD::Px)) * polyPosPow<UInt(POLYCARD::Px - dx)>::eval(P.getX()) *
    polyFactorial<dy>::eval(UInt(POLYCARD::Py)) * polyPosPow<UInt(POLYCARD::Py - dy)>::eval(P.getY()) *
    polyFactorial<dz>::eval(UInt(POLYCARD::Pz)) * polyPosPow<UInt(POLYCARD::Pz - dz)>::eval(P.getZ()) );
  }
};



//STATIC ITERATOR FOR DEGREE DETERMINATION---------------------------------------------------------
template<typename POLYCARD, UInt S> class polyDegree;

template<typename POLYCARD, UInt S>
struct polyDegree
{
  typedef typename POLYCARD::SON SON; 
  
  static UInt maxDegree()
  {
    return(std::max(polyDegree<SON,S-1>::maxDegree() , POLYCARD::Px + POLYCARD::Py + POLYCARD::Pz) );
  }
};

template<typename POLYCARD>
struct polyDegree<POLYCARD,1>
{
  static UInt maxDegree()
  {
    return(POLYCARD::Px + POLYCARD::Py + POLYCARD::Pz);
  }
};


#endif
