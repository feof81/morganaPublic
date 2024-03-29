/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef DATABASELIN1D_H
#define DATABASELIN1D_H

#include "morganaTypes.hpp"
#include "simpleFormats.hpp"
#include <map>


/*! Data base 1d with linear interpolation  */
template<typename X, typename Y>
class dataBaseLin1d
{
    /*! @name Internal data */ //@{
  public:
    Real scaleX, scaleY;
    std::map<X,Y> data;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    dataBaseLin1d();
    dataBaseLin1d(const sVect<X> & xx, const sVect<Y> & yy);
    dataBaseLin1d(const std::map<X,Y> & Data);
    dataBaseLin1d(const dataBaseLin1d & A);
    dataBaseLin1d & operator=(const dataBaseLin1d & A);
    //@}
    
    /*! @name Data management */ //@{
  public:
    void clear();
    void addPoint(const X & x, const Y & y);
    void addPoint(const sVect<X> & xx, const sVect<Y> & yy);
    bool replacePoint(const X & x, const Y & y);
    void setScaleX(const Real & ScaleX);
    void setScaleY(const Real & ScaleY);
    //@}
    
    /*! @name Interpolation */ //@{ 
  public:
    const X & getMinX() const;
    const X & getMaxX() const;
    bool isInternal(const X & x) const;
    Y interp(const X & x);
    //@}
    
    /*! @name Outstream operators */ //@{
  public:
    template<typename XX, typename YY>
    friend ostream & operator<<(ostream & f, const dataBaseLin1d<XX,YY> & M);
    
    size_t memSize() const;
    //@}
};


//_______________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------------
template<typename X, typename Y>
dataBaseLin1d<X,Y>::
dataBaseLin1d()
{
  scaleX  = 1.0;
  scaleY  = 1.0;
}

template<typename X, typename Y>
dataBaseLin1d<X,Y>::
dataBaseLin1d(const sVect<X> & xx, const sVect<Y> & yy)
{
  assert(xx.size() == yy.size());
  
  scaleX  = 1.0;
  scaleY  = 1.0;
  
  std::pair<X,Y> paio;
  
  for(UInt i=1; i <= xx.size(); ++i)
  {
    paio.first  = xx(i);
    paio.second = yy(i);
    data.insert(paio);
  }
}

template<typename X, typename Y>
dataBaseLin1d<X,Y>::
dataBaseLin1d(const std::map<X,Y> & Data)
{
  scaleX  = 1.0;
  scaleY  = 1.0;
  data    = Data;
}

template<typename X, typename Y>
dataBaseLin1d<X,Y>::
dataBaseLin1d(const dataBaseLin1d & A)
{
  scaleX  = A.scaleX;
  scaleY  = A.scaleY;
  data    = A.data;
}

template<typename X, typename Y>
dataBaseLin1d<X,Y> &
dataBaseLin1d<X,Y>::
operator=(const dataBaseLin1d & A)
{
  scaleX  = A.scaleX;
  scaleY  = A.scaleY;
  data    = A.data;
  
  return(*this);
}


//_______________________________________________________________________________________________________
// DATA MANAGEMENT
//-------------------------------------------------------------------------------------------------------
template<typename X, typename Y>
void
dataBaseLin1d<X,Y>::
clear()
{
  data.clear();
}

template<typename X, typename Y>
void
dataBaseLin1d<X,Y>::
addPoint(const X & x, const Y & y)
{
  std::pair<X,Y> paio;
  paio.first  = x;
  paio.second = y;
  
  data.insert(paio);
}

template<typename X, typename Y>
void
dataBaseLin1d<X,Y>::
addPoint(const sVect<X> & xx, const sVect<Y> & yy)
{
  assert(xx.size() == yy.size());
  
  std::pair<X,Y> paio;
  
  for(UInt i=1; i <= xx.size(); ++i)
  {
    paio.first  = xx(i);
    paio.second = yy(i);
    data.insert(paio);
  }
}

template<typename X, typename Y>
bool
dataBaseLin1d<X,Y>::
replacePoint(const X & x, const Y & y)
{
  typedef typename std::map<X,Y>::iterator ITER;
  
  ITER iter;
  bool isPresent = (data.count(x) == 1);
  
  if(isPresent)
  {
    iter = data.find(x);
    iter->second = y;
  }
  
  return(isPresent);
}

template<typename X, typename Y>
void
dataBaseLin1d<X,Y>::
setScaleX(const Real & ScaleX)
{
  scaleX = ScaleX;
}

template<typename X, typename Y>
void
dataBaseLin1d<X,Y>::
setScaleY(const Real & ScaleY)
{
  scaleY = ScaleY;
}


//_______________________________________________________________________________________________________
// INTERPOLATION
//-------------------------------------------------------------------------------------------------------
template<typename X, typename Y>
const X &
dataBaseLin1d<X,Y>::
getMinX() const
{
  return(data.begin()->first);
}

template<typename X, typename Y>
const X &
dataBaseLin1d<X,Y>::
getMaxX() const
{
  return(data.rbegin()->first);
}

template<typename X, typename Y>
bool
dataBaseLin1d<X,Y>::
isInternal(const X & x) const
{
  return( (x >= data.begin()->first) && (x <= data.rbegin()->first) );
}

template<typename X, typename Y>
Y
dataBaseLin1d<X,Y>::
interp(const X & x)
{
  typedef typename std::map<X,Y>::iterator ITER;
  
  Real realX = x * scaleX;
  
  if(x <=  data.begin()->first)
  { return(data.begin()->second); }
  
  if(x >= data.rbegin()->first)
  { return(data.rbegin()->second); }
  
  ITER itMax = data.upper_bound(realX);
  ITER itMin = itMax;
  itMin--;
  
  X xMin = itMin->first;
  Y yMin = itMin->second;
  
  X xMax = itMax->first;
  Y yMax = itMax->second;
  
  Y outY = yMax * ((realX - xMin) / (xMax - xMin)) +
           yMin * ((realX - xMax) / (xMin - xMax));
    
  return(outY * scaleY);
}


//_______________________________________________________________________________________________________
// OUTSTREAM OPERATORS
//-------------------------------------------------------------------------------------------------------
template<typename XX, typename YY>
ostream &
operator<<(ostream & f, const dataBaseLin1d<XX,YY> & M)
{
  typedef typename std::map<XX,YY>::const_iterator ITER;
  
  f << "Database" << endl;
  
  for(ITER it = M.data.begin(); it != M.data.end(); ++it)
  {
    f << it->first << " " << it->second << endl;
  }
  
  return(f);
}

template<typename X, typename Y>
size_t
dataBaseLin1d<X,Y>::
memSize() const
{
  if(data.size() == 0)
  { return(2 * sizeof(Real)); }
  else 
  {
    return(2 * sizeof(Real) + (dynamicSizeOf(data.begin()->first) 
                             + dynamicSizeOf(data.begin()->second)
                             + 2 * sizeof(Real*) ) * data.size() );
  }
}

#endif

