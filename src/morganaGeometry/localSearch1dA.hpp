/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef LOCALSEARCH1DA_HPP
#define LOCALSEARCH1DA_HPP

#include "mesh1d.hpp"
#include "connect1d.hpp"
#include "geoMapSupport1d.hpp"
#include "morganaTypes.hpp"
#include "searchData.hpp"
#include "geoOctTree.hpp"

using namespace std;
namespace mpi = boost::mpi;


//_________________________________________________________________________________________________
// SEARCHING CLASS
//-------------------------------------------------------------------------------------------------

/*! Local-serial search algorithm, 1d geometries */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class localSearch1dA
{
    /*! @name Typedefs */ //@{
  public:
    typedef GEOSHAPE               GEOSHAPE1D;
    typedef geoElement<GEOSHAPE1D> GEOELEMENT1D;
    
    typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP>     MESH1D;
    typedef connect1d<GEOSHAPE1D,ELMAP,NODEMAP>  CONNECT1D;
    
    typedef geoOctTree<MESH1D,CONNECT1D> OCTTREE1D;
    
    typedef searchData<ELMAP>           SEARCHDATA;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 1;
    //@}
    
     /*! @name Internal flags */ //@{
  public:
    bool startupLocal;
    bool meshLoaded;
    //@}
    
    /*! @name Internal structures */ //@{
  public:
    geoMapSupport1d<GEOSHAPE1D> geoSupport1d;
    UInt localQuery, numFeasibElements;
    //@}
    
    /*! @name Local internal data */ //@{
  public:
    OCTTREE1D octTree1d;                       //! octTree2d for 1d geometry namagement
    Real Xmax, Xmin, Ymax, Ymin, Zmax, Zmin;   //! coordinate of the bounding box 
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<const MESH1D>    grid1d;
    Teuchos::RCP<const CONNECT1D> connectGrid1d;
    //@}
    
    /*! @name Constructor and setting functions */ //@{
  public:
    localSearch1dA();
    void setMesh(const Teuchos::RCP<const MESH1D>    & Grid1d,
                 const Teuchos::RCP<const CONNECT1D> & ConnectGrid1d);
    
    void setMesh(const MESH1D    & Grid1d,
                 const CONNECT1D & ConnectGrid1d);
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
localSearch1dA<GEOSHAPE,ELMAP,NODEMAP>::
localSearch1dA()
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  startupLocal  = false;
  meshLoaded    = false;
  
  localQuery        = 0;
  numFeasibElements = 0;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
localSearch1dA<GEOSHAPE,ELMAP,NODEMAP>::
setMesh(const Teuchos::RCP<const MESH1D> & Grid1d, const Teuchos::RCP<const CONNECT1D> & ConnectGrid1d)
{ 
  grid1d        = Grid1d;
  connectGrid1d = ConnectGrid1d;
  
  octTree1d.setMesh(Grid1d,ConnectGrid1d);
  
  meshLoaded    = true;
  startupLocal  = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
localSearch1dA<GEOSHAPE,ELMAP,NODEMAP>::
setMesh(const MESH1D & Grid1d, const CONNECT1D & ConnectGrid1d)
{
  grid1d        = Teuchos::rcp(new MESH1D(Grid1d));
  connectGrid1d = Teuchos::rcp(new CONNECT1D(ConnectGrid1d));
  
  octTree1d.setMesh(Grid1d,ConnectGrid1d);
  
  meshLoaded    = true;
  startupLocal  = false;
}


//_________________________________________________________________________________________________
// STARTUP FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
localSearch1dA<GEOSHAPE,ELMAP,NODEMAP>::
localInit(const Real & OctToll)
{
  //Assert and controls--------------------------------------------------------
  assert(meshLoaded);
  startupLocal = true;
  
  //Startup tree---------------------------------------------------------------
  octTree1d.localInit(OctToll);
  
  //Compute boundaing box------------------------------------------------------
  point3d P;
  
  if(grid1d->getNumNodes() != 0)
  {
    Xmax = grid1d->getNodeL(1).getX(); Xmin = grid1d->getNodeL(1).getX();
    Ymax = grid1d->getNodeL(1).getY(); Ymin = grid1d->getNodeL(1).getY();
    Zmax = grid1d->getNodeL(1).getZ(); Zmin = grid1d->getNodeL(1).getZ();
  }
  
  for(UInt i=1; i <= grid1d->getNumNodes(); ++i)
  {
    P = grid1d->getNodeL(i);
      
    Xmax = max(P.getX(),Xmax);
    Ymax = max(P.getY(),Ymax);
    Zmax = max(P.getZ(),Zmax);
      
    Xmin = min(P.getX(),Xmin);
    Ymin = min(P.getY(),Ymin);
    Zmin = min(P.getZ(),Zmin);
  }
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
typename localSearch1dA<GEOSHAPE,ELMAP,NODEMAP>::SEARCHDATA
localSearch1dA<GEOSHAPE,ELMAP,NODEMAP>::
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
    points = grid1d->getElementNodesL(el);

    geoSupport1d.setPoints(points);
    prjFlag = geoSupport1d.projection(P,Y);
    assert(prjFlag);
      
    Pg = geoSupport1d.getPosition(grid1d->getElementNodesL(el), Y);
    d  = point3d::norm2(Pg - P);

    if(d < outData.getDistance())
    {
      outData.setElMap(grid1d->getElements().getRowMapL(el)); 
      outData.setLocCoord(Y);
      outData.setDistance(d);
      outData.setFound(true);
    }
  }
  
  return(outData);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
typename localSearch1dA<GEOSHAPE,ELMAP,NODEMAP>::SEARCHDATA
localSearch1dA<GEOSHAPE,ELMAP,NODEMAP>::
findLocal(const point3d & P)
{
  assert(startupLocal);
  localQuery++;
  
  sVect<UInt> matchingElements = octTree1d.getMatchingElementsExt(P);
  numFeasibElements += matchingElements.size();
  
  return( findLocal(matchingElements,P) );
}

#endif
