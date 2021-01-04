/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef LGMSEARCHDATA_HPP
#define LGMSEARCHDATA_HPP

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "typesInterface.hpp"
#include "morganaTypes.hpp"
#include "searchData.hpp"

#include <assert.h>
#include <iostream>

using namespace std;

/*! Contains the data relative to a given point: \c elMap the mapItem of the element, \c locCoord the local coordinates,
 \c the distance (if the point is not internal in the mesh), \c found whether the point has been found. */
template<typename MAP>
class lgmSearchData
{
    /*! @name Typedefs */ //@{
  public:
    typedef searchData<MAP> PDATA;
    //@}
  
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    /*! @name Internal data */ //@{
  public:
    MAP            elMap;
    sVect<point3d> locCoords;
    bool           found;
    
    sVect<PDATA> pData;
    bool isNested;
    //@}
    
    /*! @name Constructors and access functions */ //@{
  public:
    lgmSearchData();
    lgmSearchData(const lgmSearchData & data);
    //@}
    
    /*! @name Set functions */ //@{
  public:
    void setElMap(const MAP & ElMap);
    void setLocCoords(const sVect<point3d> & LocCoords);
    void setFound(const bool & Found);
    void setPData(const sVect<PDATA> & PData);
    void setIsNested(const bool & IsNested);
    void clearLocCoords();
    void clearPData();
    //@}
    
    /*! @name Get functions */ //@{
  public:
    MAP            & getElMap();
    sVect<point3d> & getLocCoords();
    point3d        & getLocCoord(const UInt & i);
    bool           & getFound();
    sVect<PDATA>   & getPData();
    PDATA          & getPData(const UInt & i);
    bool           & getIsNested();
    //@}
    
    /*! @name Const get functions */ //@{
  public:
    const MAP            & getElMap() const;
    const sVect<point3d> & getLocCoords() const;
    const point3d        & getLocCoord(const UInt & i) const;
    const bool           & getFound() const;
    const sVect<PDATA>   & getPData() const;
    const PDATA          & getPData(const UInt & i) const;
    const bool           & getIsNested() const;
    //@}
    
    /*! @name Logic operators */ //@{
  public:
    /*! Equality operator */
    lgmSearchData & operator=(const lgmSearchData & E);
    
    /*! Less operator */
    bool operator<(const lgmSearchData & E) const;
    
    /*! Not equal operator */
    bool operator!=(const lgmSearchData & E) const;
    //@}
    
    /*! @name Outstream operators */ //@{
  public:
    template<typename T>
    friend ostream & operator<<(ostream & f, const lgmSearchData<T> & M);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename MAP>
lgmSearchData<MAP>::
lgmSearchData() : found(false), isNested(true)
{
}

template<typename MAP>
lgmSearchData<MAP>::
lgmSearchData(const lgmSearchData & data)
{
  elMap     = data.elMap;
  locCoords = data.locCoords;
  found     = data.found;
  
  pData    = data.pData;
  isNested = data.isNested;
}
  
template<typename MAP>
void
lgmSearchData<MAP>::
setElMap(const MAP & ElMap)
{
  elMap = ElMap;
}

template<typename MAP>
void
lgmSearchData<MAP>::
setLocCoords(const sVect<point3d> & LocCoords)
{
  locCoords = LocCoords;
}

template<typename MAP>
void
lgmSearchData<MAP>::
setFound(const bool & Found)
{
  found = Found;
}

template<typename MAP>
void
lgmSearchData<MAP>::
setPData(const sVect<PDATA> & PData)
{
  pData = PData;
}

template<typename MAP>
void
lgmSearchData<MAP>::
setIsNested(const bool & IsNested)
{
  isNested = IsNested;
}

template<typename MAP>
void
lgmSearchData<MAP>::
clearLocCoords()
{
  locCoords.clear();
}

template<typename MAP>
void
lgmSearchData<MAP>::
clearPData()
{
  pData.clear();
}

template<typename MAP>
MAP &
lgmSearchData<MAP>::
getElMap()
{
  return(elMap);
}

template<typename MAP>
sVect<point3d> &
lgmSearchData<MAP>::
getLocCoords()
{
  return(locCoords);
}

template<typename MAP>
point3d &
lgmSearchData<MAP>::
getLocCoord(const UInt & i)
{
  assert(i <= locCoords.size());
  return(locCoords(i));
}

template<typename MAP>
bool &
lgmSearchData<MAP>::
getFound()
{
  return(found);
}

template<typename MAP>
sVect<typename lgmSearchData<MAP>::PDATA> &
lgmSearchData<MAP>::
getPData()
{
  return(pData);
}

template<typename MAP>
typename lgmSearchData<MAP>::PDATA &
lgmSearchData<MAP>::
getPData(const UInt & i)
{
  assert(i <= pData.size());
  return(pData(i));
}

template<typename MAP>
bool &
lgmSearchData<MAP>::
getIsNested()
{
  return(isNested);
}

template<typename MAP>
const MAP &
lgmSearchData<MAP>::
getElMap() const
{
  return(elMap);
}

template<typename MAP>
const sVect<point3d> &
lgmSearchData<MAP>::
getLocCoords() const
{
  return(locCoords);
}

template<typename MAP>
const point3d &
lgmSearchData<MAP>::
getLocCoord(const UInt & i) const
{
  assert(i <= locCoords.size());
  return(locCoords(i));
}

template<typename MAP>
const bool &
lgmSearchData<MAP>::
getFound() const
{
  return(found);
}

template<typename MAP>
const sVect<typename lgmSearchData<MAP>::PDATA> &
lgmSearchData<MAP>::
getPData() const
{
  return(pData);
}

template<typename MAP>
const typename lgmSearchData<MAP>::PDATA &
lgmSearchData<MAP>::
getPData(const UInt & i) const
{
  assert(i <= pData.size());
  return(pData(i));
}

template<typename MAP>
const bool &
lgmSearchData<MAP>::
getIsNested() const
{
  return(isNested);
}

template<typename MAP>
lgmSearchData<MAP> &
lgmSearchData<MAP>::
operator=(const lgmSearchData & E)
{
  elMap     = E.elMap;
  locCoords = E.locCoords;
  found     = E.found;
  
  pData    = E.pData;
  isNested = E.isNested;
  
  return(*this);
}

template<typename MAP>
bool
lgmSearchData<MAP>::
operator<(const lgmSearchData & E) const
{
  return(elMap.operator<(E));
}

template<typename MAP>
bool
lgmSearchData<MAP>::
operator!=(const lgmSearchData & E) const
{
  return(elMap.operator!=(E));
}

template<typename MAP>
template<class ARK>
void
lgmSearchData<MAP>::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & elMap;
  ar & locCoords;
  ar & found;
  ar & pData;
  ar & isNested;
}

template<typename T>
ostream & operator<<(ostream & f, const lgmSearchData<T> & M)
{
  f << endl;
  f << "Element           : " << M.elMap;
  f << "Local coordinates : " << M.locCoords << endl;
  f << "Found             : " << M.found     << endl;
  f << "point Data        : " << M.pData     << endl;
  f << "isNested          : " << M.isNested  << endl;
  
  return(f);
}

#endif
