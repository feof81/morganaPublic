/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef POINTELEMENT_HPP
#define POINTELEMENT_HPP

#include "pGraphItem.h"
#include "geoShapes.h"


/*! Element defined by its points */
template<typename GEOSHAPE>
class pointElement : public GEOSHAPE
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
  
    /*! @name Internal data */ //@{
  public:
    /*! Static data */
    static const UInt numVertices = GEOSHAPE::numVertices;
    static const UInt numPoints   = GEOSHAPE::numPoints;
      
    /*! The geometrical id */
    UInt geoId;
    
    /*! The points list */
    sVect<point3d> points;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    pointElement();
    pointElement(const UInt & GeoId, const sVect<point3d> & Points);
    pointElement(const pointElement & P);
    pointElement & operator=(const pointElement & P);
    //@}
    
    /*! @name Sort functions */ //@{
  public:
    void reorder();
    bool isOrdered() const;
    bool operator<(const pointElement & E) const;
    bool operator!=(const pointElement & E) const;
    //@}
    
    /*! @name Other functions */ //@{
  public:
    void setGeoId(const UInt & GeoId);
    void setPoint(const UInt & i, const point3d & P);
    void setPoints(const sVect<point3d> & Points);
    const UInt    & getGeoId() const;
          UInt    & getGeoId();
    const point3d & getPoint(const UInt & i) const;
          point3d & getPoint(const UInt & i);
    const sVect<point3d> & getPoints() const;
          sVect<point3d> & getPoints();
    //@}
  
    /*! @name Outstream operator */ //@{
  public:
    template<typename T>
    friend ostream & operator<<(ostream & f, const pointElement<T> & P);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE>
pointElement<GEOSHAPE>::
pointElement()
{
  points.resize(GEOSHAPE::numPoints);
}

template<typename GEOSHAPE>
pointElement<GEOSHAPE>::
pointElement(const UInt & GeoId, const sVect<point3d> & Points)
{
  assert(Points.size() == GEOSHAPE::numPoints);
  points.resize(GEOSHAPE::numPoints);
  
  geoId  = GeoId;
  points = Points;
}
    
template<typename GEOSHAPE>
pointElement<GEOSHAPE>::
pointElement(const pointElement & P)
{
  points.resize(GEOSHAPE::numPoints);
  
  geoId  = P.geoId;
  points = P.points;
}
    
template<typename GEOSHAPE>
pointElement<GEOSHAPE> &
pointElement<GEOSHAPE>::
operator=(const pointElement & P)
{
  points.resize(GEOSHAPE::numPoints);
  
  geoId  = P.geoId;
  points = P.points;
  
  return(*this);
}


//_________________________________________________________________________________________________
// SORT FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE>
void
pointElement<GEOSHAPE>::
reorder()
{
  typedef std::set<point3d>::iterator ITER;
  std::set<point3d> orderSet;
  
  for(UInt i=1; i <= points.size(); ++i)
  { orderSet.insert(points(i)); }
  
  UInt k = 1;
  
  for(ITER iter = orderSet.begin(); iter != orderSet.end(); ++iter)
  {
    points(k) = *iter;
    ++k;
  }
}

template<typename GEOSHAPE>
bool
pointElement<GEOSHAPE>::
isOrdered() const
{
  typedef std::set<point3d>::iterator ITER;
  std::set<point3d> orderSet;
  
  for(UInt i=1; i <= points.size(); ++i)
  { orderSet.insert(points(i)); }
  
  UInt k = 1;
  bool ok = true;
  
  for(ITER iter = orderSet.begin(); iter != orderSet.end(); ++iter)
  {
    ok = ok & (!(points(k) != *iter));
    ++k;
  }
  
  return(ok);
}

template<typename GEOSHAPE>
bool
pointElement<GEOSHAPE>::
operator<(const pointElement & E) const
{
  bool less, noteq;
  
  for(UInt i=1; i <= E.points.size(); ++i)
  {
    less  = points(i) <  E.points(i);
    noteq = points(i) != E.points(i);
    
    if(less)
    { return(true); }
    
    if( (!less) && (noteq) )
    { return(false); }
  }
  
  return(false);
}

template<typename GEOSHAPE>
bool
pointElement<GEOSHAPE>::
operator!=(const pointElement & E) const
{
  bool item, flag = true;
  
  for(UInt i=1; i <= E.points.size(); ++i)
  {
    item = (points(i) != E.points(i));
    flag = flag & (!item);
  }
  
  return(!flag);
}


//_________________________________________________________________________________________________
// OTHER FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE>
void
pointElement<GEOSHAPE>::
setGeoId(const UInt & GeoId)
{
  geoId = GeoId;
}

template<typename GEOSHAPE>
void
pointElement<GEOSHAPE>::
setPoint(const UInt & i, const point3d & P)
{
  assert(i > 0);
  assert(i <= GEOSHAPE::numPoints);
  
  points(i) = P;
}

template<typename GEOSHAPE>
void
pointElement<GEOSHAPE>::
setPoints(const sVect<point3d> & Points)
{
  assert(Points.size() <= GEOSHAPE::numPoints);
  
  points = Points;
}

template<typename GEOSHAPE>
const UInt &
pointElement<GEOSHAPE>::
getGeoId() const
{
  return(geoId);
}

template<typename GEOSHAPE>
UInt &
pointElement<GEOSHAPE>::
getGeoId()
{
  return(geoId);
}

template<typename GEOSHAPE>
const point3d &
pointElement<GEOSHAPE>::
getPoint(const UInt & i) const
{
  assert(i > 0);
  assert(i <= GEOSHAPE::numPoints);
  
  return(points(i));
}

template<typename GEOSHAPE>
point3d &
pointElement<GEOSHAPE>::
getPoint(const UInt & i)
{
  assert(i > 0);
  assert(i <= GEOSHAPE::numPoints);
  
  return(points(i));
}

template<typename GEOSHAPE>
const sVect<point3d> &
pointElement<GEOSHAPE>::
getPoints() const
{
  return(points);
}

template<typename GEOSHAPE>
sVect<point3d> &
pointElement<GEOSHAPE>::
getPoints()
{
  return(points);
}


//_________________________________________________________________________________________________
// OUTSTREAM
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE>
template<class ARK>
void
pointElement<GEOSHAPE>::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & geoId;
  ar & points;
}

template<typename GEOSHAPE>
ostream & operator<<(ostream & f, const pointElement<GEOSHAPE> & P)
{
  f << "geoId:  " << P.geoId << endl;
  f << "points: " << P.points;
  return(f);
}


#endif
