/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef SUPPORTHYBRID3D_HPP
#define SUPPORTHYBRID3D_HPP

#include "mesh2d.hpp"
#include "mesh3d.hpp"
#include "connect2d.hpp"
#include "connect3d.hpp"


/*! Support for hybrid fe builder */
template<typename GEOSHAPE3D, typename PMAPTYPE>
class supportHybrid3d
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename GEOSHAPE3D::GEOBSHAPE          GEOSHAPE2D;
    typedef mesh3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE>    MESH3D;
    typedef connect3d<GEOSHAPE3D,PMAPTYPE,PMAPTYPE> CONNECT3D;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    geoMapInterface<GEOSHAPE2D> surfInterface;
    UInt el3d, locFc;
    sVect<point3d> refFacePoints;     /*! Local  coordinates of the local face \c locFc of the reference element */
    sVect<point3d> globFacePoints;    /*! Global coordinates of the local face \c locFc of the element \c el3d */
    sVect<point3d> globElementPoints; /*! Global coordinates of the element \c el3d */
    //@}
    
    /*! @name Links */ //@{
  public:
    bool geometryLoaded;
    Teuchos::RCP<MESH3D>     grid3d;
    Teuchos::RCP<CONNECT3D>  connectGrid3d;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    supportHybrid3d();
    supportHybrid3d(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d);
    supportHybrid3d(MESH3D & Grid3d, CONNECT3D & ConnectGrid3d);
    void setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d);
    void setGeometry(MESH3D & Grid3d, CONNECT3D & ConnectGrid3d);
    //@}
    
    /*! @name Operation functions */ //@{
  public:
    void setElement3d(const UInt & El);
    void setLocalFace(const UInt & Face);
    point3d mapVolumeY(const point3d & Ys) const;
    point3d computeNormal(const point3d & Ys) const;
    const sVect<point3d> & getGlobFacePoints() const;
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE3D, typename PMAPTYPE>
supportHybrid3d<GEOSHAPE3D,PMAPTYPE>::
supportHybrid3d()
{
  geometryLoaded = false;
  
  el3d  = 0;
  locFc = 0;
  
  refFacePoints.resize(GEOSHAPE2D::numPoints);
  globFacePoints.resize(GEOSHAPE2D::numPoints);
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
supportHybrid3d<GEOSHAPE3D,PMAPTYPE>::
supportHybrid3d(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d)
{
  geometryLoaded = true;
  
  el3d  = 0;
  locFc = 0;
  
  refFacePoints.resize(GEOSHAPE2D::numPoints);
  globFacePoints.resize(GEOSHAPE2D::numPoints);
  
  grid3d        = Grid3d;
  connectGrid3d = ConnectGrid3d;
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
supportHybrid3d<GEOSHAPE3D,PMAPTYPE>::
supportHybrid3d(MESH3D & Grid3d, CONNECT3D & ConnectGrid3d)
{
  geometryLoaded = true;
  
  el3d  = 0;
  locFc = 0;
  
  refFacePoints.resize(GEOSHAPE2D::numPoints);
  globFacePoints.resize(GEOSHAPE2D::numPoints);
  
  grid3d        = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid3d = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
void
supportHybrid3d<GEOSHAPE3D,PMAPTYPE>::
setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d)
{
  geometryLoaded = true;
  
  grid3d        = Grid3d;
  connectGrid3d = ConnectGrid3d;
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
void
supportHybrid3d<GEOSHAPE3D,PMAPTYPE>::
setGeometry(MESH3D & Grid3d, CONNECT3D & ConnectGrid3d)
{
  geometryLoaded = true;
  
  grid3d        = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid3d = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
}



//_________________________________________________________________________________________________
// OPERATIONS FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE3D, typename PMAPTYPE>
void
supportHybrid3d<GEOSHAPE3D,PMAPTYPE>::
setElement3d(const UInt & El)
{
  assert(geometryLoaded);
  assert(El >= 1);
  assert(El <= grid3d->getNumElements());
  
  el3d = El;
  globElementPoints = grid3d->getElementNodesL(el3d);
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
void
supportHybrid3d<GEOSHAPE3D,PMAPTYPE>::
setLocalFace(const UInt & Face)
{
  assert(geometryLoaded);
  assert(Face >= 1);
  assert(Face <= GEOSHAPE3D::numFaces);
  
  locFc = Face;
  
  for(UInt i = 1; i <= surfInterface.getNumPoints(); ++i )
  {
    globFacePoints(i) = globElementPoints( GEOSHAPE3D::faceToPoint(locFc,i) );
    refFacePoints(i)  = GEOSHAPE3D::getRefNodes( GEOSHAPE3D::faceToPoint(locFc,i) );
  }
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
point3d
supportHybrid3d<GEOSHAPE3D,PMAPTYPE>::
mapVolumeY(const point3d & Ys) const
{
  assert(geometryLoaded);
  return(surfInterface.getPosition(refFacePoints,Ys));
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
point3d
supportHybrid3d<GEOSHAPE3D,PMAPTYPE>::
computeNormal(const point3d & Ys) const
{
  assert(geometryLoaded);
  point3d N = surfInterface.getDerX(globFacePoints,Ys) ^ surfInterface.getDerY(globFacePoints,Ys);
  
  return(N / N.norm2());
}

template<typename GEOSHAPE3D, typename PMAPTYPE>
const sVect<point3d> &
supportHybrid3d<GEOSHAPE3D,PMAPTYPE>::
getGlobFacePoints() const 
{
  assert(geometryLoaded);
  return(globFacePoints);
}


#endif
