/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef LOCALSEARCH2DA_HPP
#define LOCALSEARCH2DA_HPP

#include "mesh1d.hpp"
#include "mesh2d.hpp"
#include "geoMapSupport1d.hpp"
#include "geoMapSupport2d.hpp"
#include "morganaTypes.hpp"
#include "searchData.hpp"
#include "geoOctTree.hpp"

using namespace std;
namespace mpi = boost::mpi;


//_________________________________________________________________________________________________
// SEARCHING CLASS
//-------------------------------------------------------------------------------------------------

/*! Local-serial search algorithm, 2d geometries */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class localSearch2dA
{
    /*! @name Typedefs */ //@{
  public:
    typedef GEOSHAPE               GEOSHAPE2D;
    typedef geoElement<GEOSHAPE2D> GEOELEMENT2D;
    
    typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>     MESH2D;
    typedef connect2d<GEOSHAPE2D,ELMAP,NODEMAP>  CONNECT2D;
    
    typedef geoOctTree<MESH2D,CONNECT2D> OCTTREE2D;
    
    typedef searchData<ELMAP>           SEARCHDATA;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 2;
    //@}
    
     /*! @name Internal flags */ //@{
  public:
    bool startupLocal;
    bool meshLoaded;
    //@}
    
    /*! @name Internal structures */ //@{
  public:
    geoMapSupport2d<GEOSHAPE2D> geoSupport2d;
    UInt localQuery, numFeasibElements;
    //@}
    
    /*! @name Local internal data */ //@{
  public:
    OCTTREE2D octTree2d;                       //! octTree2d for 2d geometry namagement
    Real Xmax, Xmin, Ymax, Ymin, Zmax, Zmin;   //! coordinate of the bounding box 
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<const MESH2D>    grid2d;
    Teuchos::RCP<const CONNECT2D> connectGrid2d;
    //@}
    
    /*! @name Constructor and setting functions */ //@{
  public:
    localSearch2dA();
    void setMesh(const Teuchos::RCP<const MESH2D>    & Grid2d,
                 const Teuchos::RCP<const CONNECT2D> & ConnectGrid2d);
    
    void setMesh(const MESH2D    & Grid2d,
                 const CONNECT2D & ConnectGrid2d);
    //@}
    
     /*! @name Startup and search function */ //@{
  public:
    void localInit(const Real & OctToll = 2.0);
    SEARCHDATA findLocal(const sVect<UInt> & matchingElements, const point3d & P);
    SEARCHDATA findLocal(const point3d & P);
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
localSearch2dA<GEOSHAPE,ELMAP,NODEMAP>::
localSearch2dA()
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  
  startupLocal  = false;
  meshLoaded    = false;
  
  localQuery        = 0;
  numFeasibElements = 0;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
localSearch2dA<GEOSHAPE,ELMAP,NODEMAP>::
setMesh(const Teuchos::RCP<const MESH2D>    & Grid2d,
        const Teuchos::RCP<const CONNECT2D> & ConnectGrid2d)
{
  grid2d        = Grid2d;
  connectGrid2d = ConnectGrid2d;
  
  octTree2d.setMesh(Grid2d,ConnectGrid2d);
  
  meshLoaded    = true;
  startupLocal  = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
localSearch2dA<GEOSHAPE,ELMAP,NODEMAP>::
setMesh(const MESH2D    & Grid2d,
        const CONNECT2D & ConnectGrid2d)
{
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnectGrid2d));
  
  octTree2d.setMesh(Grid2d,ConnectGrid2d);
  
  meshLoaded    = true;
  startupLocal  = false;
}


//_________________________________________________________________________________________________
// STARTUP FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
localSearch2dA<GEOSHAPE,ELMAP,NODEMAP>::
localInit(const Real & OctToll)
{
  //Assert and controls--------------------------------------------------------
  assert(meshLoaded);
  startupLocal = true;
  
  //Startup tree---------------------------------------------------------------
  octTree2d.localInit(OctToll);
  
  //Compute boundaing box------------------------------------------------------
  point3d P;
  
  if(grid2d->getNumNodes() != 0)
  {
    Xmax = grid2d->getNodeL(1).getX(); Xmin = grid2d->getNodeL(1).getX();
    Ymax = grid2d->getNodeL(1).getY(); Ymin = grid2d->getNodeL(1).getY();
    Zmax = grid2d->getNodeL(1).getZ(); Zmin = grid2d->getNodeL(1).getZ();
  }
  
  for(UInt i=1; i <= grid2d->getNumNodes(); ++i)
  {
    P = grid2d->getNodeL(i);
      
    Xmax = max(P.getX(),Xmax);
    Ymax = max(P.getY(),Ymax);
    Zmax = max(P.getZ(),Zmax);
      
    Xmin = min(P.getX(),Xmin);
    Ymin = min(P.getY(),Ymin);
    Zmin = min(P.getZ(),Zmin);
  }
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
typename localSearch2dA<GEOSHAPE,ELMAP,NODEMAP>::SEARCHDATA
localSearch2dA<GEOSHAPE,ELMAP,NODEMAP>::
findLocal(const sVect<UInt> & matchingElements, const point3d & P)
{
  bool prjFlag, isInternal;
  UInt el;
  point3d Pg, Y, Pmin, Pmax;
  sVect<point3d> points;
  Real d = numeric_limits<Real>::infinity();
  
  SEARCHDATA outData;
  outData.setDistance(d);
  
  for(UInt i=1; i <= matchingElements.size(); ++i)
  {
    el     = matchingElements(i);
    points = grid2d->getElementNodesL(el);

    geoSupport2d.setPoints(points);
    prjFlag = geoSupport2d.projection(P,Y);
    assert(prjFlag);
      
    Pg = geoSupport2d.getPosition(grid2d->getElementNodesL(el), Y);
    d  = point3d::norm2(Pg - P);

    if(d < outData.getDistance())
    {
      outData.setElMap(grid2d->getElements().getRowMapL(el)); 
      outData.setLocCoord(Y);
      outData.setDistance(d);
      outData.setFound(true);
    }
  }
  
  return(outData);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
typename localSearch2dA<GEOSHAPE,ELMAP,NODEMAP>::SEARCHDATA
localSearch2dA<GEOSHAPE,ELMAP,NODEMAP>::
findLocal(const point3d & P)
{
  assert(startupLocal);
  localQuery++;
  
  sVect<UInt> matchingElements = octTree2d.getMatchingElementsExt(P);
  numFeasibElements += matchingElements.size();
  
  return( findLocal(matchingElements,P) );
}

#endif
