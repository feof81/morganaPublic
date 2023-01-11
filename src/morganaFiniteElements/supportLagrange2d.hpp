/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef SUPPORTLAGRANGE2D_HPP
#define SUPPORTLAGRANGE2D_HPP

#include "mesh1d.hpp"
#include "mesh2d.hpp"
#include "connect1d.hpp"
#include "connect2d.hpp"

#include "geoMapSupport1d.hpp"
#include "geoMapSupport2d.hpp"


/*! Support for Lagrange integration 2d */
template<typename GEOSHAPE2D, typename PMAPTYPE>
class supportLagrange2d
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename GEOSHAPE2D::GEOBSHAPE          GEOSHAPE1D;
    typedef mesh1d<GEOSHAPE1D,PMAPTYPE,PMAPTYPE>    MESH1D;
    typedef mesh2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>    MESH2D;
    typedef connect1d<GEOSHAPE1D,PMAPTYPE,PMAPTYPE> CONNECT1D;
    typedef connect2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> CONNECT2D;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    geoMapSupport2d<GEOSHAPE2D> faceInterface;
    geoMapSupport1d<GEOSHAPE1D> edgeInterface;
    UInt el2d, el1d, locEd;
    sVect<point3d> refEdgePoints;     /*! Local  coordinates of the local edge \c locEd of the reference element */
    sVect<point3d> globEdgePoints;    /*! Global coordinates of the local edge \c locEd of the element \c el2d */
    sVect<point3d> globElementPoints; /*! Global coordinates of the element \c el2d */
    sVect<point3d> globSurfacePoints; /*! Global coordinates ot the surface element 1d */
    //@}
    
     /*! @name Links */ //@{
  public:
    bool geometryLoaded;
    Teuchos::RCP<MESH2D>     grid2d;
    Teuchos::RCP<MESH1D>     grid1d;
    Teuchos::RCP<CONNECT2D>  connectGrid2d;
    Teuchos::RCP<CONNECT1D>  connectGrid1d;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    supportLagrange2d();
    supportLagrange2d(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnectGrid1d);
    supportLagrange2d(MESH2D & Grid2d, CONNECT2D & ConnectGrid2d, MESH1D & Grid1d, CONNECT1D & ConnectGrid1d);
    void setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnectGrid1d);
    void setGeometry(MESH2D & Grid2d, CONNECT2D & ConnectGrid2d, MESH1D & Grid1d, CONNECT1D & ConnectGrid1d);
    //@}
    
    /*! @name Operation functions */ //@{
  public:
    void setElement2d(const UInt & El);
    void setLocalEdge(const UInt & Edge);
    bool isBoundary() const;
    UInt getElement1d() const;
    
    point3d mapVolumeY(const point3d & Yd) const;     /*! Volume  -> element2d */
    point3d mapSurfaceY(const point3d & Yd) const;    /*! Surface -> element1d */
    point3d computeNormal(const point3d & Yd) const;  /*! Normal computed with the edge nodes */
    const sVect<point3d> & getGlobEdgePoints() const; /*! Edge -> edge */
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE2D, typename PMAPTYPE>
supportLagrange2d<GEOSHAPE2D,PMAPTYPE>::
supportLagrange2d()
{
  geometryLoaded = false;
  
  el1d  = 0;
  el2d  = 0;
  locEd = 0;
  
  refEdgePoints.resize(GEOSHAPE1D::numPoints);
  globEdgePoints.resize(GEOSHAPE1D::numPoints);
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
supportLagrange2d<GEOSHAPE2D,PMAPTYPE>::
supportLagrange2d(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnectGrid1d)
{
  geometryLoaded = true;
  
  el1d  = 0;
  el2d  = 0;
  locEd = 0;
  
  refEdgePoints.resize(GEOSHAPE1D::numPoints);
  globEdgePoints.resize(GEOSHAPE1D::numPoints);
  
  grid1d        = Grid1d;
  grid2d        = Grid2d;
  connectGrid1d = ConnectGrid1d;
  connectGrid2d = ConnectGrid2d;
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
supportLagrange2d<GEOSHAPE2D,PMAPTYPE>::
supportLagrange2d(MESH2D & Grid2d, CONNECT2D & ConnectGrid2d, MESH1D & Grid1d, CONNECT1D & ConnectGrid1d)
{
  geometryLoaded = true;
  
  el1d  = 0;
  el2d  = 0;
  locEd = 0;
  
  refEdgePoints.resize(GEOSHAPE1D::numPoints);
  globEdgePoints.resize(GEOSHAPE1D::numPoints);
  
  grid1d        = Teuchos::rcp(new MESH1D(Grid1d));
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  connectGrid1d = Teuchos::rcp(new CONNECT1D(ConnectGrid1d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
void
supportLagrange2d<GEOSHAPE2D,PMAPTYPE>::
setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnectGrid1d)
{
  geometryLoaded = true;
  
  grid1d        = Grid1d;
  grid2d        = Grid2d;
  connectGrid1d = ConnectGrid1d;
  connectGrid2d = ConnectGrid2d;
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
void
supportLagrange2d<GEOSHAPE2D,PMAPTYPE>::
setGeometry(MESH2D & Grid2d, CONNECT2D & ConnectGrid2d, MESH1D & Grid1d, CONNECT1D & ConnectGrid1d)
{
  geometryLoaded = true;
  
  grid1d        = Teuchos::rcp(new MESH1D(Grid1d));
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  connectGrid1d = Teuchos::rcp(new CONNECT1D(ConnectGrid1d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
}



//_________________________________________________________________________________________________
// OPERATIONS FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE2D, typename PMAPTYPE>
void
supportLagrange2d<GEOSHAPE2D,PMAPTYPE>::
setElement2d(const UInt & El)
{
  assert(geometryLoaded);
  assert(El >= 1);
  assert(El <= grid2d->getNumElements());
  
  el2d = El;
  globElementPoints = grid2d->getElementNodesL(el2d);
  
  faceInterface.setPoints(globElementPoints);
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
void
supportLagrange2d<GEOSHAPE2D,PMAPTYPE>::
setLocalEdge(const UInt & Edge)
{
  assert(geometryLoaded);
  assert(Edge >= 1);
  assert(Edge <= GEOSHAPE2D::numEdges);
  
  locEd = Edge;
  
  if( connectGrid2d->getEdgeIsBoundary(connectGrid2d->getElementToEdge(el2d,locEd)) )
  {
    el1d = connectGrid2d->getEdgeBEdge(connectGrid2d->getElementToEdge(el2d,locEd));
    globSurfacePoints = grid1d->getElementNodesL(el1d);
  }
  
  for(UInt i = 1; i <= edgeInterface.getNumPoints(); ++i )
  {
    globEdgePoints(i) = globElementPoints( GEOSHAPE2D::edgeToPoint(locEd,i) );
    refEdgePoints(i)  = GEOSHAPE2D::getRefNodes( GEOSHAPE2D::edgeToPoint(locEd,i) );
  }
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
bool
supportLagrange2d<GEOSHAPE2D,PMAPTYPE>::
isBoundary() const
{
  assert(geometryLoaded);
  return( connectGrid2d->getEdgeIsBoundary(connectGrid2d->getElementToEdge(el2d,locEd)) );
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
UInt
supportLagrange2d<GEOSHAPE2D,PMAPTYPE>::
getElement1d() const
{
  assert(geometryLoaded);
  return(el1d);
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
point3d
supportLagrange2d<GEOSHAPE2D,PMAPTYPE>::
mapVolumeY(const point3d & Yd) const
{
  assert(geometryLoaded);
  return(edgeInterface.getPosition(refEdgePoints,Yd));
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
point3d
supportLagrange2d<GEOSHAPE2D,PMAPTYPE>::
mapSurfaceY(const point3d & Yd) const
{
  assert(geometryLoaded);
  
  //Global mapping
  point3d P = edgeInterface.getPosition(globEdgePoints,Yd);
  
  //Local mapping  
  point3d Y;
  edgeInterface.setPoints(globSurfacePoints);
  bool ok = edgeInterface.mapGlobalToLocal(P,Y);
  
  //Return
  assert(ok);
  assert( !(P != edgeInterface.getPosition(globSurfacePoints,Y)) );
  return(Y);
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
point3d
supportLagrange2d<GEOSHAPE2D,PMAPTYPE>::
computeNormal(const point3d & Yd) const
{
  assert(geometryLoaded);
  return(faceInterface.normal(Yd,locEd));
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
const sVect<point3d> &
supportLagrange2d<GEOSHAPE2D,PMAPTYPE>::
getGlobEdgePoints() const
{
  assert(geometryLoaded);
  return(globEdgePoints);
}

#endif
