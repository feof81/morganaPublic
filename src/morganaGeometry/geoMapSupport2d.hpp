/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef GEOMAPSUPPORT2D_HPP
#define GEOMAPSUPPORT2D_HPP

#include "typesInterface.hpp"
#include "intStaticCards.h"

#include "geoMapInterface.hpp"


/*! Geometry support class for all the two-dimensional \c GEOSHAPE s */
template<typename GEOSHAPE>
class geoMapSupport2d: public geoMapInterface<GEOSHAPE>
{
    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 2;
    //@}

    /*! @name Internal data */ //@{
  public:
    sVect<point3d> points;
    //@}
  
    /*! @name Constructors and set functions */ //@{
  public:
    geoMapSupport2d();
    void setPoints(const sVect<point3d> & Points);
    //@}
    
    /*! @name Functions */ //@{
  public:
    /*! Given a point \c P finds the local coordinates \c Y. This algorithm finds
    the \c Y local coordinate such that the mapped point M(Y) has a minimum distance with respect to P.
    A nonlinear projection method is applied and a success flag is returned. */
    bool mapGlobalToLocal(const point3d & P, point3d & Y) const;
    
    /*! First the \c mapGlobalToLocal algorithm is applied and the minimum distance from the 2d object to
    the point \c P is found. Then it is checked whether the \c Y coordinate is or not inside the reference
    geometric element. In other terms this algorithm checks whether the projection on the 2d surface stays
    inside the 2d element */
    bool isInternal(const point3d & P, bool & IsInternal) const;
    
    /*! First the \c mapGlobalToLocal algorithm is applied and the minimum distance from the 2d object to
    the point \c P is found. Then the \c Y local coordinate is projected inside the reference element */
    bool projection(const point3d & P, point3d & Y) const;
    
    /*! Compute the normal to a given local edge.
    The normal is both tangent to the 2d element and normal to the 1d boundary element.
    The normal points outward the 2d element.
    \param Ys edge local coordinate
    \param Edge edge local id */
    point3d normal(const point3d & Ys, const UInt & Edge);
    
    /*! Estimates, using the numerical integration, the area of the 2d element */
    template<intTypes TYPE, UInt PRECISION>
    Real volume() const;
    //@}
    
    /*! @name Special functions */ //@{
  public:
    template<intTypes TYPE, UInt PRECISION>
    void surfaceNormal(point3d & N, Real & surface);
    //@}
};

template<typename GEOSHAPE>
geoMapSupport2d<GEOSHAPE>::
geoMapSupport2d()
{
  assert(GEOSHAPE::nDim == 2);
}

template<typename GEOSHAPE>
void
geoMapSupport2d<GEOSHAPE>::
setPoints(const sVect<point3d> & Points)
{
  points = Points;
  assert(GEOSHAPE::numPoints == points.size());
}

template<typename GEOSHAPE>
bool
geoMapSupport2d<GEOSHAPE>::
mapGlobalToLocal(const point3d & P, point3d & Y) const
{
  UInt nIter = 400;
  Real toll  = 1e-12;
  point3d Pk, R, Px, Py;
  tensor3d T;
  
  //2d problem -> set third component to zero
  Y.setZ(0.0);  
  
  for(UInt i=1; i<=nIter; ++i)
  {   
    //Compute Pk
    Pk = geoMapInterface<GEOSHAPE>::getPosition(points,Y);
    
    //Compute first order derivatives
    Px = geoMapInterface<GEOSHAPE>::getDerX(points,Y);
    Py = geoMapInterface<GEOSHAPE>::getDerY(points,Y);
    
    //Compute residual
    R.setX( - point3d::dot(Px,Pk-P) );
    R.setY( - point3d::dot(Py,Pk-P) );
    R.setZ(0.0);
    
    if( (point3d::norm2(R) <= toll) && (i >= 2) )
    { return(true); }
    
    //Compute tensor 
    T(1,1) = point3d::dot(geoMapInterface<GEOSHAPE>::template getDerivative<2,0,0>(points,Y), Pk - P) + point3d::dot(Px,Px);
    T(1,2) = point3d::dot(geoMapInterface<GEOSHAPE>::template getDerivative<1,1,0>(points,Y), Pk - P) + point3d::dot(Px,Py);
    T(2,1) = T(1,2);
    T(2,2) = point3d::dot(geoMapInterface<GEOSHAPE>::template getDerivative<0,2,0>(points,Y), Pk - P) + point3d::dot(Py,Py);
    T(3,3) = 1.0;
    
    T.computeInverse();
    
    //Local coordinate update
    Y = Y + T.secondIndexSaturation(R);
  }
  
  return(false);
}

template<typename GEOSHAPE>
bool
geoMapSupport2d<GEOSHAPE>::
isInternal(const point3d & P, bool & IsInternal) const
{
  point3d Y;
  bool flag  = mapGlobalToLocal(P,Y);
  IsInternal = GEOSHAPE::isInside(Y);
  
  return(flag);
}
    
template<typename GEOSHAPE>
bool
geoMapSupport2d<GEOSHAPE>::
projection(const point3d & P, point3d & Y) const
{
  bool flag  = mapGlobalToLocal(P,Y);
  
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
geoMapSupport2d<GEOSHAPE>::
normal(const point3d & Ys, const UInt & Edge)
{
  typedef typename GEOSHAPE::GEOBSHAPE BSHAPE;
  
  assert(Edge >= 1);
  assert(Edge <= geoMapInterface<GEOSHAPE>::getNumEdges());
  assert(GEOSHAPE::numPoints == points.size());
  
  //Surface interface
  geoMapInterface<BSHAPE> refSurfShape;
  
  //Local points
  sVect<point3d> edgePoints(geoMapInterface<GEOSHAPE>::getPointsOnEdge(Edge));
  sVect<point3d> refEdgePoints(geoMapInterface<GEOSHAPE>::getPointsOnEdge(Edge));
  
  for(UInt i=1; i <= geoMapInterface<GEOSHAPE>::getPointsOnEdge(Edge); ++i)
  {
    edgePoints(i)    = points(geoMapInterface<GEOSHAPE>::edgeToPoint(Edge,i));
    refEdgePoints(i) = geoMapInterface<GEOSHAPE>::getRefNodes(geoMapInterface<GEOSHAPE>::edgeToPoint(Edge,i));
  }
  
  //Compute normal
  point3d Yv = refSurfShape.getPosition(refEdgePoints,Ys);  
  point3d N  = geoMapInterface<GEOSHAPE>::getDerX(points,Yv) ^ geoMapInterface<GEOSHAPE>::getDerY(points,Yv);
  point3d T  = refSurfShape.getDerX(edgePoints,Ys) ^ N;

  return(T / T.norm2());
}

template<typename GEOSHAPE>
template<intTypes TYPE, UInt PRECISION>
Real
geoMapSupport2d<GEOSHAPE>::
volume() const
{
  typedef intStaticCard<GEOSHAPE,TYPE,PRECISION> INTTAB;
  
  Real vol = 0.0;
  point3d Vx, Vy;
  
  for(UInt i=1; i <= INTTAB::N; ++i)
  {
    Vx = getDerX(points,INTTAB::getYn(i));
    Vy = getDerY(points,INTTAB::getYn(i));
    
    vol += point3d::norm2(Vx ^ Vy) * INTTAB::getWn(i);
  }
  
  return(vol);
}

template<typename GEOSHAPE>
template<intTypes TYPE, UInt PRECISION>
void
geoMapSupport2d<GEOSHAPE>::
surfaceNormal(point3d & N, Real & surface)
{
  typedef intStaticCard<GEOSHAPE,TYPE,PRECISION> INTTAB;
  
  //Set to zero the components
  N.setX(0.0); N.setY(0.0); N.setZ(0.0);
  surface = 0.0;
  
  //Allocations
  point3d Vx, Vy;
  
  for(UInt i=1; i <= INTTAB::N; ++i)
  {
    Vx = geoMapInterface<GEOSHAPE>::getDerX(points,INTTAB::getYn(i));
    Vy = geoMapInterface<GEOSHAPE>::getDerY(points,INTTAB::getYn(i));
    
    N += (Vx ^ Vy) * INTTAB::getWn(i);
  }
  
  surface = point3d::norm2(N);
  N       = N / surface;
}

#endif
