/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef GEOMAPINTERFACE_HPP
#define GEOMAPINTERFACE_HPP

#include "polyStatic.hpp"
#include "geoShapes.h"


//_________________________________________________________________________________________________
// STATIC ITERATION FUNCTION FOR EVALUATION
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, UInt S>
struct baseEvalIter
{
  typedef typename geoMap<GEOSHAPE,S>::BASE POLYCARD;
  
  static point3d eval(const sVect<point3d> & points, const point3d & Y)
  {
    return(points(S) * polyStatic<POLYCARD>::evaluateStatic(Y) +
           baseEvalIter<GEOSHAPE,S-1>::eval(points,Y));
  }
};

template<typename GEOSHAPE>
struct baseEvalIter<GEOSHAPE, 1>
{
  typedef typename geoMap<GEOSHAPE,1>::BASE POLYCARD;
  
  static point3d eval(const sVect<point3d> & points, const point3d & Y)
  {
    return(points(1) * polyStatic<POLYCARD>::evaluateStatic(Y) );
  }
};



//_________________________________________________________________________________________________
// STATIC ITERATION FUNCTION FOR GRADIENT - X
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, UInt S>
struct baseDerXIter
{
  typedef typename geoMap<GEOSHAPE,S>::BASE POLYCARD;
  
  static point3d eval(const sVect<point3d> & points, const point3d & Y)
  {
    return(points(S) * polyStatic<POLYCARD>::evaluateGradientX(Y) +
           baseDerXIter<GEOSHAPE,S-1>::eval(points,Y));
  }
};

template<typename GEOSHAPE>
struct baseDerXIter<GEOSHAPE, 1>
{
  typedef typename geoMap<GEOSHAPE,1>::BASE POLYCARD;
  
  static point3d eval(const sVect<point3d> & points, const point3d & Y)
  {
    return(points(1) * polyStatic<POLYCARD>::evaluateGradientX(Y) );
  }
};



//_________________________________________________________________________________________________
// STATIC ITERATION FUNCTION FOR GRADIENT - Y
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, UInt S>
struct baseDerYIter
{
  typedef typename geoMap<GEOSHAPE,S>::BASE POLYCARD;
  
  static point3d eval(const sVect<point3d> & points, const point3d & Y)
  {
    return(points(S) * polyStatic<POLYCARD>::evaluateGradientY(Y) +
           baseDerYIter<GEOSHAPE,S-1>::eval(points,Y));
  }
};

template<typename GEOSHAPE>
struct baseDerYIter<GEOSHAPE, 1>
{
  typedef typename geoMap<GEOSHAPE,1>::BASE POLYCARD;
  
  static point3d eval(const sVect<point3d> & points, const point3d & Y)
  {
    return(points(1) * polyStatic<POLYCARD>::evaluateGradientY(Y) );
  }
};



//_________________________________________________________________________________________________
// STATIC ITERATION FUNCTION FOR GRADIENT - Z
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, UInt S>
struct baseDerZIter
{
  typedef typename geoMap<GEOSHAPE,S>::BASE POLYCARD;
  
  static point3d eval(const sVect<point3d> & points, const point3d & Y)
  {
    return(points(S) * polyStatic<POLYCARD>::evaluateGradientZ(Y) +
           baseDerZIter<GEOSHAPE,S-1>::eval(points,Y));
  }
};

template<typename GEOSHAPE>
struct baseDerZIter<GEOSHAPE, 1>
{
  typedef typename geoMap<GEOSHAPE,1>::BASE POLYCARD;
  
  static point3d eval(const sVect<point3d> & points, const point3d & Y)
  {
    return(points(1) * polyStatic<POLYCARD>::evaluateGradientZ(Y) );
  }
};


//_________________________________________________________________________________________________
// STATIC ITERATION FUNCTION FOR GENERIC DERIVATIVE
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, UInt S, UInt dx, UInt dy, UInt dz>
struct baseDerIter
{
  typedef typename geoMap<GEOSHAPE,S>::BASE POLYCARD;
  
  static point3d eval(const sVect<point3d> & points, const point3d & Y)
  {
    return(points(S) * polyStatic<POLYCARD>::template evaluateDerivative<dx,dy,dz>(Y) +
           baseDerIter<GEOSHAPE,S-1,dx,dy,dz>::eval(points,Y));
  }
};

template<typename GEOSHAPE, UInt dx, UInt dy, UInt dz>
struct baseDerIter<GEOSHAPE, 1, dx, dy, dz>
{
  typedef typename geoMap<GEOSHAPE,1>::BASE POLYCARD;
  
  static point3d eval(const sVect<point3d> & points, const point3d & Y)
  {
    return(points(1) * polyStatic<POLYCARD>::template evaluateDerivative<dx,dy,dz> (Y) );
  }
};


//_________________________________________________________________________________________________
// MAXIMUM
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, UInt S>
struct maxIter
{
  typedef typename geoMap<GEOSHAPE,S>::BASE POLYCARD;
  
  static UInt getDegree()
  {
    return(max(polyStatic<POLYCARD>::getDegree(),
	      maxIter<GEOSHAPE,S-1>::getDegree()));
  }
};

template<typename GEOSHAPE>
struct maxIter<GEOSHAPE, 1>
{
  typedef typename geoMap<GEOSHAPE,1>::BASE POLYCARD;
  
  static UInt getDegree()
  {
    return(polyStatic<POLYCARD>::getDegree());
  }
};


//_________________________________________________________________________________________________
// CLASS INTERFACE
//-------------------------------------------------------------------------------------------------

/*! Geometric interface class, produces information of the real-reference element mapping. */
template<typename GEOSHAPE>
class geoMapInterface : public GEOSHAPE
{ 
    /*! @name Constructors */ //@{
  public:
    /*! The constructor */
    geoMapInterface();
    //@}
    
    /*! @name Static get functions */ //@{
  public:
    /*! Return the index of the reference shape */
    static UInt getReferenceShapes();
    
    /*! Return the index of the reference geometry */
    static UInt getReferenceGeometry();
    
    /*! Return the element dimension */
    static UInt getNDim();
    
    /*! Return the number of the faces */
    static UInt getNumFaces();
    
    /*! Return the number of the edges */
    static UInt getNumEdges();
    
    /*! Return the number of the vertices */
    static UInt getNumVertices();
    
    /*! Return the number of the points (is equal to the number of the basis) */
    static UInt getNumPoints();
    
    /*! Returns for each \c geotype the number of elements: 
    VERTEX (the number of vertices), EDGE (the number of edges), FACE (the number of faces) */
    static UInt getNumGeoTypes(const ReferenceGeometry & geoType);
    //@}
    
    /*! @name Map interface functions functions */ //@{
  public:
    /*! Returns the position in the real element 
    \param points the list of the points of the real element 
    \param Y the coordinate point in the reference element 
    \return the position in the real element */
    point3d getPosition(const sVect<point3d> & points, const point3d & Y) const;
    
    /*! Computes the Jacobian of the trasformation, the first coloumn represents the derivative
    with respect to x, the second with respect to y and the last is the derivative with respect 
    to z
    \param points the list of the points of the real element 
    \param Y the coordinate point in the reference element 
    \return the jacobian */
    tensor3d getGradient(const sVect<point3d> & points, const point3d & Y) const;
    
    /*! Computes the determinant of the jacobian
    \param points the list of the points of the real element 
    \param Y the coordinate point in the reference element 
    \return the determinant of the jacobian */
    Real getGradientDet(const sVect<point3d> & points, const point3d & Y) const;
    
    /*! Computes the derivative of the map with respect to x
    \param points the list of the points of the real element 
    \param Y the coordinate point in the reference element 
    \return the derivative vector */
    point3d getDerX(const sVect<point3d> & points, const point3d & Y) const;
    
    /*! Computes the derivative of the map with respect to y
    \param points the list of the points of the real element 
    \param Y the coordinate point in the reference element 
    \return the derivative vector */
    point3d getDerY(const sVect<point3d> & points, const point3d & Y) const;
    
    /*! Computes the derivative of the map with respect to z
    \param points the list of the points of the real element 
    \param Y the coordinate point in the reference element 
    \return the derivative vector */
    point3d getDerZ(const sVect<point3d> & points, const point3d & Y) const;
    
    /*! Generic deirvative evaluation
    \param points the list of the points of the real element 
    \param Y the coordinate point in the reference element 
    \return the derivative vector */
    template<UInt dx, UInt dy, UInt dz>
    point3d getDerivative(const sVect<point3d> & points, const point3d & Y) const;
    
    /*! Returns the maximum degree of the map polynomial */
    UInt getDegree() const;
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE>
geoMapInterface<GEOSHAPE>::
geoMapInterface()
{ }

template<typename GEOSHAPE>
UInt
geoMapInterface<GEOSHAPE>::
getReferenceShapes()
{
  return(GEOSHAPE::Shape);
}

template<typename GEOSHAPE>
UInt
geoMapInterface<GEOSHAPE>::
getReferenceGeometry()
{
  return(GEOSHAPE::Geometry);
}

template<typename GEOSHAPE>
UInt
geoMapInterface<GEOSHAPE>::
getNDim()
{
  return(GEOSHAPE::nDim);
}

template<typename GEOSHAPE>
UInt
geoMapInterface<GEOSHAPE>::
getNumFaces()
{
  return(GEOSHAPE::numFaces);
}

template<typename GEOSHAPE>
UInt
geoMapInterface<GEOSHAPE>::
getNumEdges()
{
  return(GEOSHAPE::numEdges);
}

template<typename GEOSHAPE>
UInt
geoMapInterface<GEOSHAPE>::
getNumVertices()
{
  return(GEOSHAPE::numVertices);
}

template<typename GEOSHAPE>
UInt
geoMapInterface<GEOSHAPE>::
getNumPoints()
{
  return(GEOSHAPE::numPoints);
}

template<typename GEOSHAPE>
UInt
geoMapInterface<GEOSHAPE>::
getNumGeoTypes(const ReferenceGeometry & geoType)
{ 
  static const UInt v[3] = {
    GEOSHAPE::numVertices,
    GEOSHAPE::numEdges,
    GEOSHAPE::numFaces };
    
  assert(UInt(geoType) <= 3);
    
  return(v[ UInt(geoType) - 1]);
}

template<typename GEOSHAPE>
point3d
geoMapInterface<GEOSHAPE>::
getPosition(const sVect<point3d> & points, const point3d & Y) const
{
  assert(points.size() == GEOSHAPE::numPoints);
  
  return(baseEvalIter<GEOSHAPE,GEOSHAPE::numPoints>::eval(points,Y));
}

template<typename GEOSHAPE>
tensor3d
geoMapInterface<GEOSHAPE>::
getGradient(const sVect<point3d> & points, const point3d & Y) const
{
  assert(points.size() == GEOSHAPE::numPoints);
  
  return(tensor3d(
  baseDerXIter<GEOSHAPE,GEOSHAPE::numPoints>::eval(points,Y),
  baseDerYIter<GEOSHAPE,GEOSHAPE::numPoints>::eval(points,Y),
  baseDerZIter<GEOSHAPE,GEOSHAPE::numPoints>::eval(points,Y)) );
}

template<typename GEOSHAPE>
Real
geoMapInterface<GEOSHAPE>::
getGradientDet(const sVect<point3d> & points, const point3d & Y) const
{
  assert(points.size() == GEOSHAPE::numPoints);
  
  tensor3d T(
  baseDerXIter<GEOSHAPE,GEOSHAPE::numPoints>::eval(points,Y),
  baseDerYIter<GEOSHAPE,GEOSHAPE::numPoints>::eval(points,Y),
  baseDerZIter<GEOSHAPE,GEOSHAPE::numPoints>::eval(points,Y));
  
  return(T.getThirdInvariant());
}

template<typename GEOSHAPE>
point3d
geoMapInterface<GEOSHAPE>::
getDerX(const sVect<point3d> & points, const point3d & Y) const
{
  assert(points.size() == GEOSHAPE::numPoints);
  
  return(baseDerXIter<GEOSHAPE,GEOSHAPE::numPoints>::eval(points,Y));
}

template<typename GEOSHAPE>
point3d
geoMapInterface<GEOSHAPE>::
getDerY(const sVect<point3d> & points, const point3d & Y) const
{
  assert(points.size() == GEOSHAPE::numPoints);
  
  return(baseDerYIter<GEOSHAPE,GEOSHAPE::numPoints>::eval(points,Y));
}

template<typename GEOSHAPE>
point3d
geoMapInterface<GEOSHAPE>::
getDerZ(const sVect<point3d> & points, const point3d & Y) const
{
  assert(points.size() == GEOSHAPE::numPoints);
  
  return(baseDerZIter<GEOSHAPE,GEOSHAPE::numPoints>::eval(points,Y));
}

template<typename GEOSHAPE>
template<UInt dx, UInt dy, UInt dz>
point3d
geoMapInterface<GEOSHAPE>::
getDerivative(const sVect<point3d> & points, const point3d & Y) const
{
  assert(points.size() == GEOSHAPE::numPoints);
  
  return(baseDerIter<GEOSHAPE, GEOSHAPE::numPoints, dx,dy,dz>::eval(points,Y));
}

template<typename GEOSHAPE>
UInt
geoMapInterface<GEOSHAPE>::
getDegree() const
{
  return(maxIter<GEOSHAPE,GEOSHAPE::numPoints>::getDegree());
}

#endif
