/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef SEARCHDATA_HPP
#define SEARCHDATA_HPP

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "typesInterface.hpp"
#include "morganaTypes.hpp"

#include <assert.h>
#include <iostream>

using namespace std;


/*! Contains the data relative to a given point: \c elMap the mapItem of the element, \c locCoord the local coordinates,
 \c the distance (if the point is not internal in the mesh), \c found whether the point has been found. */
template<typename MAP>
class searchData
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    /*! @name Internal data */ //@{
  public:
    MAP     elMap;
    point3d locCoord;
    Real    distance;
    bool    found;
    //@}
    
    /*! @name Constructors and access functions */ //@{
  public:
    searchData();
    searchData(const searchData & data);
    //@}
    
    /*! @name Set functions */ //@{
  public:
    void setElMap(const MAP & ElMap);
    void setLocCoord(const point3d & Y);
    void setDistance(const Real & d);
    void setFound(const bool & Found);
    //@}
    
    /*! @name Get functions */ //@{
  public:
    MAP     & getElMap();
    point3d & getLocCoord();
    Real    & getDistance();
    bool    & getFound();
    //@}
    
    /*! @name Const get functions */ //@{
  public:
    const MAP     & getElMap() const;
    const point3d & getLocCoord() const;
    const Real    & getDistance() const;
    const bool    & getFound() const;
    //@}
    
    /*! @name Logic operators */ //@{
  public:
    /*! Equality operator */
    searchData & operator=(const searchData & E);
    
    /*! Less operator */
    bool operator<(const searchData & E) const;
    
    /*! Not equal operator */
    bool operator!=(const searchData & E) const;
    //@}
    
    /*! @name Outstream operators */ //@{
  public:
    template<typename T>
    friend ostream & operator<<(ostream & f, const searchData<T> & M);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename MAP>
searchData<MAP>::
searchData() : found(false)
{
}

template<typename MAP>
searchData<MAP>::
searchData(const searchData & data)
{
  elMap    = data.elMap;
  locCoord = data.locCoord;
  distance = data.distance;
  found    = data.found;
}
    
template<typename MAP>
void
searchData<MAP>::
setElMap(const MAP & ElMap)
{
  elMap = ElMap;
}

template<typename MAP>
void
searchData<MAP>::
setLocCoord(const point3d & Y)
{
  locCoord = Y;
}

template<typename MAP>
void
searchData<MAP>::
setDistance(const Real & d)
{
  distance = d;
}

template<typename MAP>
void
searchData<MAP>::
setFound(const bool & Found)
{
  found = Found;
}

template<typename MAP>
MAP &
searchData<MAP>::
getElMap()
{
  return(elMap);
}

template<typename MAP>
point3d &
searchData<MAP>::
getLocCoord()
{
  return(locCoord);
}

template<typename MAP>
Real &
searchData<MAP>::
getDistance()
{
  return(distance);
}

template<typename MAP>
bool &
searchData<MAP>::
getFound()
{
  return(found);
}

template<typename MAP>
const MAP &
searchData<MAP>::
getElMap() const
{
  return(elMap);
}

template<typename MAP>
const point3d &
searchData<MAP>::
getLocCoord() const
{
  return(locCoord);
}

template<typename MAP>
const Real &
searchData<MAP>::
getDistance() const
{
  return(distance);
}

template<typename MAP>
const bool &
searchData<MAP>::
getFound() const
{
  return(found);
}

template<typename MAP>
searchData<MAP> &
searchData<MAP>::
operator=(const searchData & E)
{
  elMap    = E.elMap;
  locCoord = E.locCoord;
  distance = E.distance;
  found    = E.found;
  
  return(*this);
}

template<typename MAP>
bool
searchData<MAP>::
operator<(const searchData & E) const
{
  return(elMap.operator<(E));
}

template<typename MAP>
bool
searchData<MAP>::
operator!=(const searchData & E) const
{
  return(elMap.operator!=(E));
}

template<typename MAP>
template<class ARK>
void
searchData<MAP>::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & elMap;
  ar & locCoord;
  ar & distance; 
  ar & found;
}

template<typename T>
ostream & operator<<(ostream & f, const searchData<T> & M)
{
  f << endl;
  f << "Element          : " << M.elMap;
  f << "Local coordinate : " << M.locCoord;
  f << "Distance         : " << M.distance << endl;
  f << "Found            : " << M.found << endl;
  return(f);
}


#endif
