/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef DOFTYPETRAITS_HPP
#define DOFTYPETRAITS_HPP

#include "typesInterface.hpp"
#include "EpetraExt_DistArray.h"

using namespace EpetraExt;

/*! Generic traits - empty */
template<typename TYPE>
class doftypeTraits;


//-------------------------------------------------------------------------------------------------
/*! Real trait */
template<> class doftypeTraits<Real>
{
  public:
    doftypeTraits();
    void paraviewString(ofstream & out, const string & value, const string & unit);
    void printDof(ofstream & out, const Real & a);
    void printDof(DistArray<double> & vectorData, const UInt & rowL, const Real & a);
};

doftypeTraits<Real>::
doftypeTraits()
{ }

void
doftypeTraits<Real>::
paraviewString(ofstream & out, const string & value, const string & unit)
{
  out << "1 1" << endl;
  out << value << ", " << unit << endl;
}

void
doftypeTraits<Real>::
printDof(ofstream & out, const Real & a)
{
  out << a << endl;
}

void
doftypeTraits<Real>::
printDof(DistArray<double> & vectorData, const UInt & rowL, const Real & a)
{
  vectorData(rowL,0) = a;
}


//-------------------------------------------------------------------------------------------------
/*! Point2d trait */
template<> class doftypeTraits<point2d>
{
  public:
    doftypeTraits();
    void paraviewString(ofstream & out, const string & value, const string & unit);
    void printDof(ofstream & out, const point2d & a);
    void printDof(DistArray<double> & vectorData, const UInt & rowL, const point2d & a);
};

doftypeTraits<point2d>::
doftypeTraits()
{ }

void
doftypeTraits<point2d>::
paraviewString(ofstream & out, const string & value, const string & unit)
{
  out << "1 3" << endl;
  out << value << ", " << unit << endl;
}

void
doftypeTraits<point2d>::
printDof(ofstream & out, const point2d & a)
{
  out << a.getX() << " " << a.getY() << " " << 0 << endl;
}

void
doftypeTraits<point2d>::
printDof(DistArray<double> & vectorData, const UInt & rowL, const point2d & a)
{
  vectorData(rowL,0) = a.getX();
  vectorData(rowL,1) = a.getY();
}


//-------------------------------------------------------------------------------------------------
/*! Tensor2d trait */
template<> class doftypeTraits<tensor2d>
{
  public:
    doftypeTraits();
    void paraviewString(ofstream & out, const string & value, const string & unit);
    void printDof(ofstream & out, const tensor2d & a);
    void printDof(DistArray<double> & vectorData, const UInt & rowL, const tensor2d & a);
};

doftypeTraits<tensor2d>::
doftypeTraits()
{ }

void
doftypeTraits<tensor2d>::
paraviewString(ofstream & out, const string & value, const string & unit)
{
  out << "4 1 1 1 1"  << endl;
  out << value << "xx" << ", " << unit << endl;
  out << value << "xy" << ", " << unit << endl;
  out << value << "yx" << ", " << unit << endl;
  out << value << "yy" << ", " << unit << endl;
}

void
doftypeTraits<tensor2d>::
printDof(ofstream & out, const tensor2d & a)
{
  out << a.getIJ(1,1) << " " << a.getIJ(1,2) << " " << a.getIJ(2,1) << " " << a.getIJ(2,2) << endl;
}

void
doftypeTraits<tensor2d>::
printDof(DistArray<double> & vectorData, const UInt & rowL, const tensor2d & a)
{
  vectorData(rowL,0) = a.getIJ(1,1);
  vectorData(rowL,1) = a.getIJ(1,2);
  vectorData(rowL,2) = a.getIJ(2,1);
  vectorData(rowL,3) = a.getIJ(2,2);
}


//-------------------------------------------------------------------------------------------------
/*! Point3d trait */
template<> class doftypeTraits<point3d>
{
  public:
    doftypeTraits();
    void paraviewString(ofstream & out, const string & value, const string & unit);
    void printDof(ofstream & out, const point3d & a);
    void printDof(DistArray<double> & vectorData, const UInt & rowL, const point3d & a);
};

doftypeTraits<point3d>::
doftypeTraits()
{ }

void
doftypeTraits<point3d>::
paraviewString(ofstream & out, const string & value, const string & unit)
{
  out << "1 3" << endl;
  out << value << ", " << unit << endl;
}

void
doftypeTraits<point3d>::
printDof(ofstream & out, const point3d & a)
{
  out << a.getX() << " " << a.getY() << " " << a.getZ() << endl;
}

void
doftypeTraits<point3d>::
printDof(DistArray<double> & vectorData, const UInt & rowL, const point3d & a)
{
  vectorData(rowL,0) = a.getX();
  vectorData(rowL,1) = a.getY();
  vectorData(rowL,2) = a.getZ();
}



//-------------------------------------------------------------------------------------------------
/*! Tensor3d trait */
template<> class doftypeTraits<tensor3d>
{
  public:
    doftypeTraits();
    void paraviewString(ofstream & out, const string & value, const string & unit);
    void printDof(ofstream & out, const tensor3d & a);
    void printDof(DistArray<double> & vectorData, const UInt & rowL, const tensor3d & a);
};

doftypeTraits<tensor3d>::
doftypeTraits()
{ }

void
doftypeTraits<tensor3d>::
paraviewString(ofstream & out, const string & value, const string & unit)
{
  out << "9 1 1 1 1 1 1 1 1 1"  << endl;
  out << value << "xx" << ", " << unit << endl;
  out << value << "xy" << ", " << unit << endl;
  out << value << "xz" << ", " << unit << endl;
  out << value << "yx" << ", " << unit << endl;
  out << value << "yy" << ", " << unit << endl;
  out << value << "yz" << ", " << unit << endl;
  out << value << "zx" << ", " << unit << endl;
  out << value << "zy" << ", " << unit << endl;
  out << value << "zz" << ", " << unit << endl;
}

void
doftypeTraits<tensor3d>::
printDof(ofstream & out, const tensor3d & a)
{
  out << a.getIJ(1,1) << " " << a.getIJ(1,2) << " " << a.getIJ(1,3)
      << a.getIJ(2,1) << " " << a.getIJ(2,2) << " " << a.getIJ(2,3)
      << a.getIJ(3,1) << " " << a.getIJ(3,2) << " " << a.getIJ(3,3) << endl;
}

void
doftypeTraits<tensor3d>::
printDof(DistArray<double> & vectorData, const UInt & rowL, const tensor3d & a)
{
  vectorData(rowL,0) = a.getIJ(1,1);
  vectorData(rowL,1) = a.getIJ(1,2);
  vectorData(rowL,2) = a.getIJ(1,3);
  
  vectorData(rowL,3) = a.getIJ(2,1);
  vectorData(rowL,4) = a.getIJ(2,2);
  vectorData(rowL,5) = a.getIJ(2,3);
  
  vectorData(rowL,6) = a.getIJ(3,1);
  vectorData(rowL,7) = a.getIJ(3,2);
  vectorData(rowL,8) = a.getIJ(3,3);
}


//-------------------------------------------------------------------------------------------------
/*! Static vector trait */
template<size_t N> class doftypeTraits<staticVector<N> >
{
  public:
    doftypeTraits();
    void paraviewString(ofstream & out, const string & value, const string & unit);
    void printDof(ofstream & out, const staticVector<N> & a);
    void printDof(DistArray<double> & vectorData, const UInt & rowL, const staticVector<N> & a);
};

template<size_t N>
doftypeTraits<staticVector<N> >::
doftypeTraits()
{ }

template<size_t N>
void
doftypeTraits<staticVector<N> >::
paraviewString(ofstream & out, const string & value, const string & unit)
{
  out << N << " ";
  
  for(UInt i=1; i <= N; ++i)
  {
    out << " " << 1;
  }
  out << endl;
  
  for(UInt i=1; i <= N; ++i)
  {
    out << value << i << ", " << unit << endl;
  }  
}

template<size_t N>
void
doftypeTraits<staticVector<N> >::
printDof(ofstream & out, const staticVector<N> & a)
{
  for(UInt i=1; i <= N; ++i)
  {
    out << a(i) << " ";
  }
  out << endl;
}

template<size_t N>
void
doftypeTraits<staticVector<N> >::
printDof(DistArray<double> & vectorData, const UInt & rowL, const staticVector<N> & a)
{
  for(UInt i=1; i <= N; ++i)
  {
    vectorData(rowL,i-1) = a(i);
  }
}

#endif
