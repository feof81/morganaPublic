/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef SUPPORTLAGRANGE3D_HPP
#define SUPPORTLAGRANGE3D_HPP

#include "morganaTypes.hpp"
#include "mesh2d.hpp"
#include "mesh3d.hpp"
#include "connect2d.hpp"
#include "connect3d.hpp"
#include "geoMapSupport2d.hpp"


/*! Support for Lagrange integration 3d */
template<typename GEOSHAPE3D, typename PMAPTYPE>
class supportLagrange3d
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename GEOSHAPE3D::GEOBSHAPE          GEOSHAPE2D;
    typedef mesh2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>    MESH2D;
    typedef mesh3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>    MESH3D;
    typedef connect2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> CONNECT2D;
    typedef connect3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> CONNECT3D;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    geoMapSupport2d<GEOSHAPE2D> surfInterface;
    UInt el3d, el2d, locFc;
    sVect<point3d> refFacePoints;       /*! Local  coordinates of the local face \c locFc of the reference element */
    sVect<point3d> globFacePoints;      /*! Global coordinates of the local face \c locFc of the element \c el3d */
    sVect<point3d> globElement3dPoints; /*! Global coordinates of the element \c el3d */
    sVect<point3d> globElement2dPoints; /*! Global coordinates ot the surface element 2d */
    //@}
    
    /*! @name Links */ //@{
  public:
    bool geometryLoaded;
    Teuchos::RCP<MESH3D>     grid3d;
    Teuchos::RCP<MESH2D>     grid2d;
    Teuchos::RCP<CONNECT3D>  connectGrid3d;
    Teuchos::RCP<CONNECT2D>  connectGrid2d;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    supportLagrange3d();
    supportLagrange3d(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    supportLagrange3d(MESH3D & Grid3d, CONNECT3D & ConnectGrid3d, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d);
    void setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    void setGeometry(MESH3D & Grid3d, CONNECT3D & ConnectGrid3d, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d);
    //@}
    
    /*! @name Operation functions */ //@{
  public:
    void setElement3d(const UInt & El);
    void setLocalFace(const UInt & Face);
    bool isBoundary() const;
    UInt getElement2d() const;
    point3d mapVolumeY(const point3d & Yf) const;     /*! Volume  -> element3d */
    point3d mapSurfaceY(const point3d & Yf);          /*! Surface -> element2d */
    point3d computeNormal(const point3d & Yf) const;  /*! Normal computed with the face nodes */
    const sVect<point3d> & getGlobFacePoints() const; /*! Face -> face */
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE3D, typename PMAPTYPE>
supportLagrange3d<GEOSHAPE3D,PMAPTYPE>::
supportLagrange3d()
{
  geometryLoaded = false;
  
  el3d  = 0;
  el2d  = 0;
  locFc = 0;
  
  refFacePoints.resize(GEOSHAPE2D::numPoints);
  globFacePoints.resize(GEOSHAPE2D::numPoints);
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
supportLagrange3d<GEOSHAPE3D,PMAPTYPE>::
supportLagrange3d(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  
  el3d  = 0;
  el2d  = 0;
  locFc = 0;
  
  refFacePoints.resize(GEOSHAPE2D::numPoints);
  globFacePoints.resize(GEOSHAPE2D::numPoints);
  
  grid2d        = Grid2d;
  grid3d        = Grid3d;
  connectGrid2d = ConnectGrid2d;
  connectGrid3d = ConnectGrid3d;
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
supportLagrange3d<GEOSHAPE3D,PMAPTYPE>::
supportLagrange3d(MESH3D & Grid3d, CONNECT3D & ConnectGrid3d, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d)
{
  geometryLoaded = true;
  
  el3d  = 0;
  el2d  = 0;
  locFc = 0;
  
  refFacePoints.resize(GEOSHAPE2D::numPoints);
  globFacePoints.resize(GEOSHAPE2D::numPoints);
  
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  grid3d        = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
  connectGrid3d = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
void
supportLagrange3d<GEOSHAPE3D,PMAPTYPE>::
setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d)
{
  geometryLoaded = true;
  
  grid2d        = Grid2d;
  grid3d        = Grid3d;
  connectGrid2d = ConnectGrid2d;
  connectGrid3d = ConnectGrid3d;
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
void
supportLagrange3d<GEOSHAPE3D,PMAPTYPE>::
setGeometry(MESH3D & Grid3d, CONNECT3D & ConnectGrid3d, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d)
{
  geometryLoaded = true;
  
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  grid3d        = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
  connectGrid3d = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
}



//_________________________________________________________________________________________________
// OPERATIONS FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE3D, typename PMAPTYPE>
void
supportLagrange3d<GEOSHAPE3D,PMAPTYPE>::
setElement3d(const UInt & El)
{
  assert(geometryLoaded);
  assert(El >= 1);
  assert(El <= grid3d->getNumElements());
  
  el3d = El;
  globElement3dPoints = grid3d->getElementNodesL(el3d);
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
void
supportLagrange3d<GEOSHAPE3D,PMAPTYPE>::
setLocalFace(const UInt & Face)
{
  assert(geometryLoaded);
  assert(Face >= 1);
  assert(Face <= GEOSHAPE3D::numFaces);
  
  locFc = Face;
  
  if( connectGrid3d->getFaceIsBoundary(connectGrid3d->getElementToFace(el3d,locFc)) )
  {
    el2d = connectGrid3d->getFaceBFace(connectGrid3d->getElementToFace(el3d,locFc));
    globElement2dPoints = grid2d->getElementNodesL(el2d);
  }
  
  for(UInt i = 1; i <= surfInterface.getNumPoints(); ++i )
  {
    globFacePoints(i) = globElement3dPoints( GEOSHAPE3D::faceToPoint(locFc,i) );
    refFacePoints(i)  = GEOSHAPE3D::getRefNodes( GEOSHAPE3D::faceToPoint(locFc,i) );
  }
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
bool
supportLagrange3d<GEOSHAPE3D,PMAPTYPE>::
isBoundary() const
{
  assert(geometryLoaded);
  return( connectGrid3d->getFaceIsBoundary(connectGrid3d->getElementToFace(el3d,locFc)) );
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
UInt
supportLagrange3d<GEOSHAPE3D,PMAPTYPE>::
getElement2d() const
{
  assert(geometryLoaded);
  return(el2d);
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
point3d
supportLagrange3d<GEOSHAPE3D,PMAPTYPE>::
mapVolumeY(const point3d & Yf) const
{
  assert(geometryLoaded);
  return(surfInterface.getPosition(refFacePoints,Yf));
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
point3d
supportLagrange3d<GEOSHAPE3D,PMAPTYPE>::
mapSurfaceY(const point3d & Yf)
{
  assert(geometryLoaded);
  
  //Global mapping
  point3d P = surfInterface.getPosition(globFacePoints,Yf);
  
  //Local mapping  
  point3d Y;
  surfInterface.setPoints(globElement2dPoints);
  bool ok = surfInterface.mapGlobalToLocal(P,Y);
  
  //Return
  assert(ok);
  assert( point3d::norm2(P - surfInterface.getPosition(globElement2dPoints,Y)) <= geoToll );
  
  return(Y);
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
point3d
supportLagrange3d<GEOSHAPE3D,PMAPTYPE>::
computeNormal(const point3d & Yf) const
{
  assert(geometryLoaded);
  point3d N = surfInterface.getDerX(globFacePoints,Yf) ^ surfInterface.getDerY(globFacePoints,Yf);
  
  return(N / N.norm2());
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
const sVect<point3d> &
supportLagrange3d<GEOSHAPE3D,PMAPTYPE>::
getGlobFacePoints() const 
{
  assert(geometryLoaded);
  return(globFacePoints);
}


#endif
