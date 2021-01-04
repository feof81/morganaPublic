/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef REFINE3D_HPP
#define REFINE3D_HPP

#include "sRefine3d.hpp"
#include "../morganaContainer/pMapItem.h"
#include "pVectGlobalManip.hpp"

/*! Parallel refinement class */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class refine3d : public sRefine3d<GEOSHAPE, ELMAP, NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:    
    typedef GEOSHAPE                       GEOSHAPE3D;
    typedef typename GEOSHAPE3D::GEOBSHAPE GEOSHAPE2D;
    typedef typename GEOSHAPE2D::GEOBSHAPE GEOSHAPE1D;
    
    typedef geoElement<GEOSHAPE3D>       GEOELEMENT3D;
    typedef geoElement<GEOSHAPE2D>       GEOELEMENT2D;
    typedef geoElement<GEOSHAPE1D>       GEOELEMENT1D;
    
    typedef pointElement<GEOSHAPE1D> POINT_EDGE;
    
    typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP>  MESH2D;
    typedef mesh3d<GEOSHAPE3D,ELMAP,NODEMAP>  MESH3D;
    
    typedef sRefine3d<GEOSHAPE, ELMAP, NODEMAP> SREFINE;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    Teuchos::RCP<MESH2D> grid2d;
    Teuchos::RCP<MESH3D> grid3d;
    //@}
    
    /*! @name Constructors and io functions */ //@{
  public:
    refine3d();
    refine3d(const communicator & CommDev);
    refine3d(const Teuchos::RCP<communicator> & CommDev);
    
    void setCommDev(const communicator & CommDev);
    void setCommDev(const Teuchos::RCP<communicator> & CommDev);
    
    void upload(const MESH3D & grid3d,
                const MESH2D & grid2d);
    
    void upload(const Teuchos::RCP<MESH3D> & grid3d,
                const Teuchos::RCP<MESH2D> & grid2d);
    
    void download(MESH3D & grid3d,
                  MESH2D & grid2d);
    
    void download(Teuchos::RCP<MESH3D> & grid3d,
                  Teuchos::RCP<MESH2D> & grid2d);
    //@}
    
    /*! @name Refine functions - parallel leb */ //@{
  public:
    void setRefinementParams(const Real & TollH,
                             const UInt & RFactor,
                             const sVect<UInt> & ElList);
    
    sVect<POINT_EDGE> getEdgesTbr(const UInt & maxLength);
    
    UInt refineEdges(const sVect<POINT_EDGE> & edges) const;
    
    bool refineLeb(const UInt & edgesBlock,
                   const UInt & numSteps);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
refine3d<GEOSHAPE,ELMAP,NODEMAP>::
refine3d()
{
  commDevLoaded = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
refine3d<GEOSHAPE,ELMAP,NODEMAP>::
refine3d(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
refine3d<GEOSHAPE,ELMAP,NODEMAP>::
refine3d(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
refine3d<GEOSHAPE,ELMAP,NODEMAP>::
setCommDev(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
refine3d<GEOSHAPE,ELMAP,NODEMAP>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}
    
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
refine3d<GEOSHAPE,ELMAP,NODEMAP>::
upload(const MESH3D & grid3d,
       const MESH2D & grid2d)
{
  //Assert-------------------------------------------------
  assert(commDevLoaded);
  
  //Typedef------------------------------------------------
  typedef sRefine3d<GEOSHAPE, ELMAP, NODEMAP> SREFINE;
  
  //Checking-----------------------------------------------
  mesh2dGlobalManip<GEOSHAPE2D,ELMAP,NODEMAP> manip2d(commDev);
  mesh3dGlobalManip<GEOSHAPE3D,ELMAP,NODEMAP> manip3d(commDev);
  
  assert(manip2d.check(grid2d));
  assert(manip3d.check(grid3d));
  
  //Startup------------------------------------------------
  SREFINE::upload(grid3d, grid2d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
refine3d<GEOSHAPE,ELMAP,NODEMAP>::
upload(const Teuchos::RCP<MESH3D> & grid3d,
       const Teuchos::RCP<MESH2D> & grid2d)
{
  //Assert-------------------------------------------------
  assert(commDevLoaded);
  
  //Typedef------------------------------------------------
  typedef sRefine3d<GEOSHAPE, ELMAP, NODEMAP> SREFINE;
  
  //Checking-----------------------------------------------
  mesh2dGlobalManip<GEOSHAPE2D,ELMAP,NODEMAP> manip2d(commDev);
  mesh3dGlobalManip<GEOSHAPE3D,ELMAP,NODEMAP> manip3d(commDev);
  
  assert(manip2d.check(grid2d));
  assert(manip3d.check(grid3d));
  
  //Startup------------------------------------------------
  SREFINE::upload(grid3d, grid2d);
}

/*template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
refine3d<GEOSHAPE,ELMAP,NODEMAP>::
download(MESH3D & grid3d,
         MESH2D & grid2d)
{
  //Assert-------------------------------------------------
  assert(commDevLoaded);
  
  //Typedef------------------------------------------------
  typedef sRefine3d<GEOSHAPE, ELMAP, NODEMAP> SREFINE;
  
  //Download-----------------------------------------------
  SREFINE::download(grid3d, grid2d);
}*/
    
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
refine3d<GEOSHAPE,ELMAP,NODEMAP>::
download(Teuchos::RCP<MESH3D> & grid3d,
         Teuchos::RCP<MESH2D> & grid2d)
{
  //Assert-------------------------------------------------
  assert(commDevLoaded);
  
  //Typedef------------------------------------------------
  typedef sRefine3d<GEOSHAPE, ELMAP, NODEMAP> SREFINE;
  typedef pGraph<GEOELEMENT3D,ELMAP,NODEMAP>  GRAPH3D;
  typedef pGraph<GEOELEMENT2D,ELMAP,NODEMAP>  GRAPH2D;
  
  //Download-----------------------------------------------
  SREFINE::download(grid3d, grid2d);
  
  grid2d->clearEdges();
  grid3d->clearFaces();
  grid3d->clearEdges();
  
  //Fix the nodes map 3d-----------------------------------
  pVectGlobalManip<point3d,NODEMAP> nodeManip(commDev);
  pVect<point3d,NODEMAP> tempNodes = grid3d->getNodes();
  
  nodeManip.buildGlobalNumbering(tempNodes);
  
  grid3d->setNodes(tempNodes);
  grid3d->transferMap();
  
  //Fix the nodes map 2d-----------------------------------
  tempNodes = grid2d->getNodes();
  
  nodeManip.buildGlobalNumbering(tempNodes);
  
  grid2d->setNodes(tempNodes);
  grid2d->transferMap();
  
  //Fix the global ids of the 3d elements------------------
  pGraphGlobalManip<GEOELEMENT3D,ELMAP,NODEMAP> elManip3d(commDev);
  GRAPH3D tempElements3d = grid3d->getElements();
  
  elManip3d.buildGlobalNumbering(tempElements3d);
  
  grid3d->setElements(tempElements3d);
  
  //Fix the global ids of the 2d elements------------------
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> elManip2d(commDev);
  GRAPH2D tempElements2d = grid2d->getElements();
  
  elManip2d.buildGlobalNumbering(tempElements2d);
  
  grid2d->setElements(tempElements2d);
}


//_________________________________________________________________________________________________
// REFINE FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
refine3d<GEOSHAPE,ELMAP,NODEMAP>::
setRefinementParams(const Real & TollH,
                    const UInt & RFactor,
                    const sVect<UInt> & ElList)
{
  typedef sRefine3d<GEOSHAPE, ELMAP, NODEMAP> SREFINE;
  
  SREFINE::setRefinementParams(TollH, RFactor, ElList);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
sVect<typename refine3d<GEOSHAPE,ELMAP,NODEMAP>::POINT_EDGE>
refine3d<GEOSHAPE,ELMAP,NODEMAP>::
getEdgesTbr(const UInt & maxLength)
{ 
  //Assert-----------------------------------------------------------
  assert(commDevLoaded);
  
  //Typedef----------------------------------------------------------
  typedef sRefine3d<GEOSHAPE, ELMAP, NODEMAP> SREFINE;
  
  //Serial edges-----------------------------------------------------
  sVect<POINT_EDGE> sEdges = SREFINE::getEdgesTbr(maxLength);
  
  //Parallel edges---------------------------------------------------
  sVect<POINT_EDGE> pEdges;
  std::vector<sVect<POINT_EDGE> > globalStack;
  
  all_gather(*commDev, sEdges, globalStack);
  
  for(UInt i=0; i<globalStack.size(); ++i)
  {
    for(UInt j=1; j <= globalStack[i].size(); ++j)
    { pEdges.push_back(globalStack[i](j)); }
  }
  
  return(pEdges);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
refine3d<GEOSHAPE,ELMAP,NODEMAP>::
refineEdges(const sVect<POINT_EDGE> & edges) const
{
  //Assert-------------------------------------------------
  assert(commDevLoaded);
  
  //Typedef----------------------------------------------------------
  typedef sRefine3d<GEOSHAPE, ELMAP, NODEMAP> SREFINE;
  
  return(SREFINE::refineEdges(edges));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
bool
refine3d<GEOSHAPE,ELMAP,NODEMAP>::
refineLeb(const UInt & edgesBlock,
          const UInt & numSteps)
{
  //Assert-------------------------------------------------
  assert(commDevLoaded);
  
  //Typedef----------------------------------------------------------
  typedef sRefine3d<GEOSHAPE, ELMAP, NODEMAP> SREFINE;
  
  //Alloc--------------------------------------------------
  UInt numCut, totalCut;
  bool success = false;
  sVect<POINT_EDGE> edges;
  boost::mpi::maximum<Real> op;
  
  //Refinement cycle---------------------------------------
  for(UInt i=1; i <= numSteps; ++i)
  {
    edges  = getEdgesTbr(edgesBlock);
    numCut = refineEdges(edges);
    
    all_reduce(*commDev, numCut, totalCut, op);
    
    if(totalCut == 0)
    {
      success = true;
      break;
    }
  }
  
  return(success);
}

#endif
