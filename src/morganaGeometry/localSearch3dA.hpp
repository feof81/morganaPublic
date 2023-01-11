/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef LOCALSEARCH3DA_HPP
#define LOCALSEARCH3DA_HPP

#include "mesh2d.hpp"
#include "mesh3d.hpp"
#include "geoMapSupport2d.hpp"
#include "geoMapSupport3d.hpp"
#include "morganaTypes.hpp"
#include "searchData.hpp"
#include "geoOctTree.hpp"
#include "searchBoundingBox.h"

using namespace std;
namespace mpi = boost::mpi;


template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class localSearch3dA
{
    /*! @name Typedefs */ //@{
  public:
    typedef GEOSHAPE                       GEOSHAPE3D;
    typedef typename GEOSHAPE::GEOBSHAPE   GEOSHAPE2D;
    typedef geoElement<GEOSHAPE2D>       GEOELEMENT2D;
    typedef geoElement<GEOSHAPE3D>       GEOELEMENT3D;
    
    typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>     MESH2D;
    typedef mesh3d<GEOSHAPE3D,ELMAP,NODEMAP>     MESH3D;
    typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>  CONNECT2D;
    typedef connect3d<GEOSHAPE3D,ELMAP,NODEMAP>  CONNECT3D;
    
    typedef geoOctTree<MESH2D,CONNECT2D> OCTTREE2D;
    typedef geoOctTree<MESH3D,CONNECT3D> OCTTREE3D;
    
    typedef searchData<ELMAP>           SEARCHDATA;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 3;
    //@}
    
    /*! @name Internal flags */ //@{
  public:
    bool startupLocal;
    bool meshLoaded;
    //@}
    
    /*! @name Internal structures */ //@{
  public:
    geoMapSupport2d<GEOSHAPE2D> geoSupport2d;
    geoMapSupport3d<GEOSHAPE3D> geoSupport3d;
    UInt localQuery, numFeasibElements;
    //@}
    
    /*! @name Local internal data */ //@{
  public:
    OCTTREE2D octTree2d;                       //! octTree2d for 2d geometry namagement
    OCTTREE3D octTree3d;                       //! octTree3d for 3d geometry namagement
    Real Xmax, Xmin, Ymax, Ymin, Zmax, Zmin;   //! coordinate of the bounding box 
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<const MESH2D>    grid2d;
    Teuchos::RCP<const MESH3D>    grid3d;
    Teuchos::RCP<const CONNECT2D> connectGrid2d;
    Teuchos::RCP<const CONNECT3D> connectGrid3d;
    //@}
    
    /*! @name Constructor and setting functions */ //@{
  public:
    localSearch3dA();
    void setMesh(const Teuchos::RCP<const MESH2D> & Grid2d,
                 const Teuchos::RCP<const MESH3D> & Grid3d,
                 const Teuchos::RCP<const CONNECT2D> & ConnectGrid2d,
                 const Teuchos::RCP<const CONNECT3D> & ConnectGrid3d);
    
    void setMesh(const MESH2D & Grid2d,
                 const MESH3D & Grid3d,
                 const CONNECT2D & ConnectGrid2d,
                 const CONNECT3D & ConnectGrid3d);
    //@}
    
     /*! @name Startup and search function */ //@{
  public:
    void       localInit(const Real & OctToll = 2.0);
    point3d    projection(const point3d & P);
    SEARCHDATA findLocal(const sVect<UInt> & matchingElements, const point3d & P, const point3d & Pj);
    SEARCHDATA findLocal(const point3d & P);
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
localSearch3dA<GEOSHAPE,ELMAP,NODEMAP>::
localSearch3dA()
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  startupLocal  = false;
  meshLoaded    = false;
  
  localQuery        = 0;
  numFeasibElements = 0;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
localSearch3dA<GEOSHAPE,ELMAP,NODEMAP>::
setMesh(const Teuchos::RCP<const MESH2D> & Grid2d,
        const Teuchos::RCP<const MESH3D> & Grid3d,
        const Teuchos::RCP<const CONNECT2D> & ConnectGrid2d,
        const Teuchos::RCP<const CONNECT3D> & ConnectGrid3d)
{
  grid2d        = Grid2d;
  grid3d        = Grid3d;
  connectGrid2d = ConnectGrid2d;
  connectGrid3d = ConnectGrid3d;
  
  octTree2d.setMesh(Grid2d,ConnectGrid2d);
  octTree3d.setMesh(Grid3d,ConnectGrid3d);
  
  meshLoaded    = true;
  startupLocal  = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
localSearch3dA<GEOSHAPE,ELMAP,NODEMAP>::
setMesh(const MESH2D & Grid2d,
        const MESH3D & Grid3d,
        const CONNECT2D & ConnectGrid2d,
        const CONNECT3D & ConnectGrid3d)
{
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  grid3d        = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
  connectGrid3d = Teuchos::rcp(new CONNECT3D(ConnectGrid3d));
  
  octTree2d.setMesh(Grid2d,ConnectGrid2d);
  octTree3d.setMesh(Grid3d,ConnectGrid3d);
  
  meshLoaded    = true;
  startupLocal  = false;
}


//_________________________________________________________________________________________________
// STARTUP FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
localSearch3dA<GEOSHAPE,ELMAP,NODEMAP>::
localInit(const Real & OctToll)
{
  //Assert and controls--------------------------------------------------------
  assert(meshLoaded);
  startupLocal = true;
  
  //Startup tree---------------------------------------------------------------
  octTree2d.localInit(OctToll);
  octTree3d.localInit(OctToll);
  
  //Compute boundaing box------------------------------------------------------
  point3d P;
  
  if(grid3d->getNumNodes() != 0)
  {
    Xmax = grid3d->getNodeL(1).getX(); Xmin = grid3d->getNodeL(1).getX();
    Ymax = grid3d->getNodeL(1).getY(); Ymin = grid3d->getNodeL(1).getY();
    Zmax = grid3d->getNodeL(1).getZ(); Zmin = grid3d->getNodeL(1).getZ();
  }
  
  for(UInt i=1; i <= grid3d->getNumNodes(); ++i)
  {
    P = grid3d->getNodeL(i);
      
    Xmax = max(P.getX(),Xmax);
    Ymax = max(P.getY(),Ymax);
    Zmax = max(P.getZ(),Zmax);
      
    Xmin = min(P.getX(),Xmin);
    Ymin = min(P.getY(),Ymin);
    Zmin = min(P.getZ(),Zmin);
  }
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
point3d
localSearch3dA<GEOSHAPE,ELMAP,NODEMAP>::
projection(const point3d & P)
{
  UInt el;
  bool prjFlag;
  point3d Pg, Pfin, Y;
  Real d, dmin = numeric_limits<Real>::infinity();
  sVect<UInt> matchingElements = octTree2d.getMatchingElementsExt(P);
  
  for(UInt i=1; i <= matchingElements.size(); ++i)
  {
    el = matchingElements(i);
    
    geoSupport2d.setPoints( grid2d->getElementNodesL(el) );
    prjFlag = geoSupport2d.projection(P,Y);
    assert(prjFlag);
    
    Pg = geoSupport2d.getPosition(grid2d->getElementNodesL(el), Y);
    d  = point3d::norm2(Pg - P);
    
    if(d < dmin)
    {
      dmin = d;
      Pfin = Pg;
    }
  }
  
  return(Pfin);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
typename localSearch3dA<GEOSHAPE,ELMAP,NODEMAP>::SEARCHDATA
localSearch3dA<GEOSHAPE,ELMAP,NODEMAP>::
findLocal(const sVect<UInt> & matchingElements, const point3d & P, const point3d & Pj)
{
  //Alloc------------------------------------------------------------
  bool prjFlag, isInternal;
  UInt el;
  point3d Pg, Y, Pmin, Pmax;
  sVect<point3d> points;
  Real d = numeric_limits<Real>::infinity();
  
  SEARCHDATA outData;
  outData.setDistance(d);
  
  //Exit-------------------------------------------------------------
  if(matchingElements.size() == 0)
  {
     outData.setLocCoord(Y);
     outData.setFound(false);
  }
  
  //Search-----------------------------------------------------------
  for(UInt i=1; i <= matchingElements.size(); ++i)
  {
    el     = matchingElements(i);
    points = grid3d->getElementNodesL(el);
    
    GEOSHAPE::boundingBox(points,Pmin,Pmax);
    isInternal = (Pj.getX() >= (Pmin.getX() - geoToll)) && (Pj.getX() <= (Pmax.getX() + geoToll)) &&
                 (Pj.getY() >= (Pmin.getY() - geoToll)) && (Pj.getY() <= (Pmax.getY() + geoToll)) &&
                 (Pj.getZ() >= (Pmin.getZ() - geoToll)) && (Pj.getZ() <= (Pmax.getZ() + geoToll));
                 
    if(isInternal)
    {
      geoSupport3d.setPoints(points);
      prjFlag = geoSupport3d.projection(Pj,Y);
      assert(prjFlag);
      
      Pg = geoSupport3d.getPosition(points, Y);
      d  = point3d::norm2(Pg - P);

      if(d < outData.getDistance())
      {
        outData.setElMap(grid3d->getElements().getRowMapL(el)); 
        outData.setLocCoord(Y);
        outData.setDistance(d);
        outData.setFound(true);
      }
    }
  } 
  return(outData);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
typename localSearch3dA<GEOSHAPE,ELMAP,NODEMAP>::SEARCHDATA
localSearch3dA<GEOSHAPE,ELMAP,NODEMAP>::
findLocal(const point3d & P)
{
  assert(startupLocal);
  localQuery++;
  
  //Find the matching elements
  point3d Pj;
  sVect<UInt> matchingElements = octTree3d.getMatchingElements(P);
  
  if(matchingElements.size() == 0)  //If the point is outside this is projected in the domain                                
  {    
                    Pj = projection(P);
      matchingElements = octTree3d.getMatchingElements(Pj);
    numFeasibElements += matchingElements.size();    
    
    SEARCHDATA outData = findLocal(matchingElements,P,Pj);
    
    assert( !( (!outData.getFound()) && (grid2d->getNumElements() != 0) )  );
    
    return(outData);
  }
  else                             //If the search cell is not empy
  {
    SEARCHDATA outData = findLocal(matchingElements,P,P);
    
    if(outData.getDistance() <= geoToll)
    {     
      numFeasibElements += matchingElements.size();
      return(outData);
    }
    else                          //If the distance is not zero a pojection attempt is carried out
    {
                      Pj = projection(P);
        matchingElements = octTree3d.getMatchingElements(Pj);
      numFeasibElements += matchingElements.size();
    
      return( findLocal(matchingElements,P,Pj) );
    }
  }
}

#endif
