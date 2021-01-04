/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef GEOMAPSUPPORT1D_HPP
#define GEOMAPSUPPORT1D_HPP

#include "typesInterface.hpp"
#include "intStaticCards.h"

#include "geoMapInterface.hpp"


/*! Geometry support class for all the one-dimensional \c GEOSHAPE s */
template<typename GEOSHAPE>
class geoMapSupport1d: public geoMapInterface<GEOSHAPE>
{
    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 1;
    //@}

    /*! @name Internal data */ //@{
  public:
    sVect<point3d> points;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    geoMapSupport1d();
    void setPoints(const sVect<point3d> & Points);
    //@}
    
    /*! @name Functions */ //@{
  public:
    /*! Given a point \c P finds the local coordinates \c Y. This algorithm finds
    the \c Y local coordinate such that the mapped point M(Y) has a minimum distance with respect to P.
    A nonlinear projection method is applied and a success flag is returned. */
    bool mapGlobalToLocal(const point3d & P, point3d & Y) const;
    
    /*! First the \c mapGlobalToLocal algorithm is applied and the minimum distance from the 1d object to
    the point \c P is found. Then it is checked whether the \c Y coordinate is or not inside the reference
    geometric element. */
    bool isInternal(const point3d & P, bool & IsInternal) const;
    
    /*! First the \c mapGlobalToLocal algorithm is applied and the minimum distance from the 1d object to
    the point \c P is found. Then the \c Y local coordinate is projected inside the reference element */
    bool projection(const point3d & P, point3d & Y) const;
    
    /*! Computes the tangent vector to the 1d interval in the positions Ys = (0,0,0), Side = 1 and 
    Ys = (1,0,0), Side = 2 */
    point3d normal(const point3d & Ys, const UInt & Side);
    
    /*! Estimates, using the numerical integration, the length of the 1d element */
    template<intTypes TYPE, UInt PRECISION> 
    Real volume() const;
    //@}
};

template<typename GEOSHAPE>
geoMapSupport1d<GEOSHAPE>::
geoMapSupport1d()
{
  assert(GEOSHAPE::nDim == 1);
}

template<typename GEOSHAPE>
void
geoMapSupport1d<GEOSHAPE>::
setPoints(const sVect<point3d> & Points)
{
  points = Points;
  assert(GEOSHAPE::numPoints == points.size());
}

template<typename GEOSHAPE>
bool
geoMapSupport1d<GEOSHAPE>::
mapGlobalToLocal(const point3d & P, point3d & Y) const
{
  UInt nIter = 400;
  Real toll  = 1e-8 * point3d::norm2(P);
  point3d Pk, Px;
  Real R, T;
  
  //1d problem -> set the second and third component to zero
  Y.setY(0.0);
  Y.setZ(0.0);
  
  for(UInt i=1; i<=nIter; ++i)
  {   
    //Compute Pk
    Pk = geoMapInterface<GEOSHAPE>::getPosition(points,Y);
    
    //Compute first order derivative
    Px = geoMapInterface<GEOSHAPE>::getDerX(points,Y);
    
    //Compute residual
    R = - point3d::dot(Px,Pk-P);
    
    if( (abs(R) <= toll) && (i >= 2))
    { return(true); }
    
    //Compute tensor 
    T = point3d::dot(geoMapInterface<GEOSHAPE>::template getDerivative<2,0,0>(points,Y), Pk - P) + point3d::dot(Px,Px);
    
    //Local coordinate update
    Y.setX( Y.getX() + (R / T) );
  }
  
  return(false);
}

template<typename GEOSHAPE>
bool
geoMapSupport1d<GEOSHAPE>::
isInternal(const point3d & P, bool & IsInternal) const
{
  point3d Y;
  bool flag  = mapGlobalToLocal(P,Y);
  IsInternal = GEOSHAPE::isInside(Y);
  
  return(flag);
}

template<typename GEOSHAPE>
bool
geoMapSupport1d<GEOSHAPE>::
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
geoMapSupport1d<GEOSHAPE>::
normal(const point3d & Ys, const UInt & Side)
{
  assert(Side >= 1);
  assert(Side <= 2);
  
  point3d N = geoMapInterface<GEOSHAPE>::getDerX(points,Ys) * ( (Side == 2) - (Side == 1) );
    
  return(N / N.norm2());
}

template<typename GEOSHAPE>
template<intTypes TYPE, UInt PRECISION>
Real
geoMapSupport1d<GEOSHAPE>::
volume() const
{
  typedef intStaticCard<GEOSHAPE,TYPE,PRECISION> INTTAB;
  
  Real vol = 0.0;
  point3d Vx;
  
  for(UInt i=1; i <= INTTAB::N; ++i)
  {
    Vx   = geoMapInterface<GEOSHAPE>::getDerX(points,INTTAB::getYn(i));
    vol += point3d::norm2(Vx) * INTTAB::getWn(i);
  }
  
  return(vol);
}

#endif

