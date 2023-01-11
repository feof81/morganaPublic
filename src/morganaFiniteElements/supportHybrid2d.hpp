/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef SUPPORTHYBRID2D_HPP
#define SUPPORTHYBRID2D_HPP

#include "mesh1d.hpp"
#include "mesh2d.hpp"
#include "connect1d.hpp"
#include "connect2d.hpp"

#include "geoMapSupport2d.hpp"


/*! Support for hybrid fe builder */
template<typename GEOSHAPE2D, typename PMAPTYPE>
class supportHybrid2d
{
  /*! @name Typedefs */ //@{
  public:
    typedef typename GEOSHAPE2D::GEOBSHAPE          GEOSHAPE1D;
    typedef mesh2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>    MESH2D;
    typedef connect2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> CONNECT2D;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    geoMapSupport2d<GEOSHAPE2D> faceInterface;
    geoMapInterface<GEOSHAPE1D> edgeInterface;
    UInt el2d, locEd;
    sVect<point3d> refEdgePoints;     /*! Local  coordinates of the local edge \c locEd of the reference element */
    sVect<point3d> globEdgePoints;    /*! Global coordinates of the local edge \c locEd of the element \c el2d */
    sVect<point3d> globElementPoints; /*! Global coordinates of the element \c el2d */
    //@}
    
    /*! @name Links */ //@{
  public:
    bool geometryLoaded;
    Teuchos::RCP<MESH2D>     grid2d;
    Teuchos::RCP<CONNECT2D>  connectGrid2d;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    supportHybrid2d();
    supportHybrid2d(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    supportHybrid2d(MESH2D & Grid2d, CONNECT2D & ConnectGrid2d);
    void setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    void setGeometry(MESH2D & Grid2d, CONNECT2D & ConnectGrid2d);
    //@}
    
    /*! @name Operation functions */ //@{
  public:
    void setElement2d(const UInt & El);
    void setLocalEdge(const UInt & Edge);
    point3d mapVolumeY(const point3d & Ys) const;
    point3d computeNormal(const point3d & Ys) const;
    const sVect<point3d> & getGlobEdgePoints() const;
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE2D, typename PMAPTYPE>
supportHybrid2d<GEOSHAPE2D,PMAPTYPE>::
supportHybrid2d()
{
  geometryLoaded = false;
  
  el2d  = 0;
  locEd = 0;
  
  refEdgePoints.resize(GEOSHAPE1D::numPoints);
  globEdgePoints.resize(GEOSHAPE1D::numPoints);
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
supportHybrid2d<GEOSHAPE2D,PMAPTYPE>::
supportHybrid2d(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  
  el2d  = 0;
  locEd = 0;
  
  refEdgePoints.resize(GEOSHAPE1D::numPoints);
  globEdgePoints.resize(GEOSHAPE1D::numPoints);
  
  grid2d        = Grid2d;
  connectGrid2d = ConnectGrid2d;
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
supportHybrid2d<GEOSHAPE2D,PMAPTYPE>::
supportHybrid2d(MESH2D & Grid2d, CONNECT2D & ConnectGrid2d)
{
  geometryLoaded = true;
  
  el2d  = 0;
  locEd = 0;
  
  refEdgePoints.resize(GEOSHAPE1D::numPoints);
  globEdgePoints.resize(GEOSHAPE1D::numPoints);
  
  grid2d        = Grid2d;
  connectGrid2d = ConnectGrid2d;
  
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
void
supportHybrid2d<GEOSHAPE2D,PMAPTYPE>::
setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  
  grid2d        = Grid2d;
  connectGrid2d = ConnectGrid2d;
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
void
supportHybrid2d<GEOSHAPE2D,PMAPTYPE>::
setGeometry(MESH2D & Grid2d, CONNECT2D & ConnectGrid2d)
{
  geometryLoaded = true;
  
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
}



//_________________________________________________________________________________________________
// OPERATIONS FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE2D, typename PMAPTYPE>
void
supportHybrid2d<GEOSHAPE2D,PMAPTYPE>::
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
supportHybrid2d<GEOSHAPE2D,PMAPTYPE>::
setLocalEdge(const UInt & Edge)
{
  assert(geometryLoaded);
  assert(Edge >= 1);
  assert(Edge <= GEOSHAPE2D::numEdges);
  
  locEd = Edge;
  
  for(UInt i = 1; i <= edgeInterface.getNumPoints(); ++i )
  {
    globEdgePoints(i) = globElementPoints( GEOSHAPE2D::edgeToPoint(locEd,i) );
    refEdgePoints(i)  = GEOSHAPE2D::getRefNodes( GEOSHAPE2D::edgeToPoint(locEd,i) );
  }
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
point3d
supportHybrid2d<GEOSHAPE2D,PMAPTYPE>::
mapVolumeY(const point3d & Ys) const
{
  assert(geometryLoaded);
  return(edgeInterface.getPosition(refEdgePoints,Ys));
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
point3d
supportHybrid2d<GEOSHAPE2D,PMAPTYPE>::
computeNormal(const point3d & Ys) const
{
  assert(geometryLoaded);
  return(faceInterface.normal(Ys,locEd));
}

template<typename GEOSHAPE2D, typename PMAPTYPE>
const sVect<point3d> &
supportHybrid2d<GEOSHAPE2D,PMAPTYPE>::
getGlobEdgePoints() const
{
  assert(geometryLoaded);
  return(globEdgePoints);
}


#endif
