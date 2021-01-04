/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef GEOMAPSUPPORT3D_HPP
#define GEOMAPSUPPORT3D_HPP

#include "typesInterface.hpp"
#include "intStaticCards.h"

#include "geoMapInterface.hpp"


/*! Geometry support class for all the one-dimensional \c GEOSHAPE s */
template<typename GEOSHAPE>
class geoMapSupport3d : public geoMapInterface<GEOSHAPE>
{
    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 3;
    //@}

    /*! @name Internal data */ //@{
  public:
    sVect<point3d> points;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    geoMapSupport3d();
    void setPoints(const sVect<point3d> & Points);
    //@}
    
    /*! @name Functions */ //@{
  public:
     /*! Given a point \c P finds the local coordinates \c Y. This algorithm finds
    the \c Y local coordinate such that the mapped point M(Y) has a minimum distance with respect to P.
    A nonlinear projection method is applied and a success flag is returned. */
    bool mapGlobalToLocal(const point3d & P, point3d & Y) const;
    
    /*! First the \c mapGlobalToLocal algorithm is applied and the minimum distance from the 3d object to
    the point \c P is found. Then it is checked whether the \c Y coordinate is or not inside the reference
    geometric element. */
    bool isInternal(const point3d & P, bool & IsInternal) const;
    
    /*! First the \c mapGlobalToLocal algorithm is applied and the minimum distance from the 3d object to
    the point \c P is found. Then the \c Y local coordinate is projected inside the reference element */
    bool projection(const point3d & P, point3d & Y) const;
    
    /*! Compute the normal to a given local face.
    The normal points outwar the 3d element.
    \param Ys face local coordinate
    \param Face face local id */
    point3d normal(const point3d & Ys, const UInt & Face) const;
    
    /*! Estimates, using the numerical integration, the volume of the 3d element */
    template<intTypes TYPE, UInt PRECISION>
    Real volume() const;
    //@}
};

template<typename GEOSHAPE>
geoMapSupport3d<GEOSHAPE>::
geoMapSupport3d()
{
  assert(GEOSHAPE::nDim == 3);
}

template<typename GEOSHAPE>
void
geoMapSupport3d<GEOSHAPE>::
setPoints(const sVect<point3d> & Points)
{
  points = Points;
  assert(GEOSHAPE::numPoints == points.size());
}

template<typename GEOSHAPE>
bool
geoMapSupport3d<GEOSHAPE>::
mapGlobalToLocal(const point3d & P, point3d & Y) const
{
  UInt nIter = 400;
  Real toll  = 1e-12;
  point3d Pk;
  tensor3d T;
  
  for(UInt i=1; i<=nIter; ++i)
  {    
    Pk = geoMapInterface<GEOSHAPE>::getPosition(points,Y);
    T  = geoMapInterface<GEOSHAPE>::getGradient(points,Y);
    
    if( (point3d::norm2(P-Pk) <= toll) && (i >= 2))
    { return(true); }
    
    T.computeInverse();
    
    Y = Y + T.secondIndexSaturation(P - Pk);
  }
   
  return(false);
}

template<typename GEOSHAPE>
bool
geoMapSupport3d<GEOSHAPE>::
isInternal(const point3d & P, bool & IsInternal) const
{
  point3d Y;
  bool flag  = mapGlobalToLocal(P,Y);
  IsInternal = GEOSHAPE::isInside(Y);
  
  return(flag);
}

template<typename GEOSHAPE>
bool
geoMapSupport3d<GEOSHAPE>::
projection(const point3d & P, point3d & Y) const
{
  bool flag = mapGlobalToLocal(P,Y);
  
  if(GEOSHAPE::isInside(Y))
  {
    return(flag);
  }
  else
  {
    point3d Q;
    GEOSHAPE::projectInside(points,P,Q);
    flag = mapGlobalToLocal(Q,Y);
  }
  
  return(flag);
}

template<typename GEOSHAPE>
point3d
geoMapSupport3d<GEOSHAPE>::
normal(const point3d & Ys, const UInt & Face) const
{
  typedef typename GEOSHAPE::GEOBSHAPE BSHAPE;
  
  assert(Face >= 1);
  assert(Face <= geoMapInterface<GEOSHAPE>::getNumFaces());
  assert(GEOSHAPE::numPoints == points.size());
  
  //Surface interface
  geoMapInterface<BSHAPE> refSurfShape;
 
  //Local points
  sVect<point3d> facePoints(geoMapInterface<GEOSHAPE>::getPointsOnFace(Face));
  
  for(UInt i=1; i <= geoMapInterface<GEOSHAPE>::getPointsOnFace(Face); ++i)
  {
    facePoints(i) = points(geoMapInterface<GEOSHAPE>::faceToPoint(Face,i));
  }
  
  //Compute normal
  point3d N = refSurfShape.getDerX(facePoints,Ys) ^ refSurfShape.getDerY(facePoints,Ys);
  
  return(N / N.norm2());
}

template<typename GEOSHAPE>
template<intTypes TYPE, UInt PRECISION>
Real
geoMapSupport3d<GEOSHAPE>::
volume() const
{
  typedef intStaticCard<GEOSHAPE,TYPE,PRECISION> INTTAB;
  
  Real vol = 0.0;
  
  for(UInt i=1; i <= INTTAB::N; ++i)
  {
    vol += geoMapInterface<GEOSHAPE>::getGradientDet(points,INTTAB::getYn(i)) * INTTAB::getWn(i);
  }
  
  return(vol);
}

#endif
