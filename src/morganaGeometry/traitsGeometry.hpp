/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef TRAITSGEOMETRY_HPP
#define TRAITSGEOMETRY_HPP

#include "pMapItem.h"
#include "pMapItemShare.h"
#include "pMap.hpp"
#include "pMapGlobalManip.h"

#include "pGraph.hpp"

#include "geoElement.hpp"

using namespace std;


//_________________________________________________________________________________________________
// AUXILIARY CHECKING FUNCTIONS FOR MESHGLOBAL MANIPULATORS -> 3D
//-------------------------------------------------------------------------------------------------

/*! Auxiliary graph checks 3d - unspecialized, empty */
template<typename GEOSHAPE, typename ROWMAP, typename COLMAP>
class auxiliaryGraphCheck3d
{ };


/*! Auxiliary graph checks 3d - specialization for \c pMapItem */
template<typename GEOSHAPE>
class auxiliaryGraphCheck3d<GEOSHAPE,pMapItemShare,pMapItem>
{
  public:
    typedef geoElement<GEOSHAPE>                      GEOELEMENT;
    typedef pGraph<GEOELEMENT,pMapItemShare,pMapItem> GRAPH;
    
  public:
    auxiliaryGraphCheck3d();
    bool check(const GRAPH & Graph);
};

template<typename GEOSHAPE>
auxiliaryGraphCheck3d<GEOSHAPE,pMapItemShare,pMapItem>::
auxiliaryGraphCheck3d()
{ }

template<typename GEOSHAPE>
bool
auxiliaryGraphCheck3d<GEOSHAPE,pMapItemShare,pMapItem>::
check(const GRAPH & Graph)
{
  return(true);
}


/*! Auxiliary graph checks 3d - specialization for \c pMapItemSharing */
template<typename GEOSHAPE>
class auxiliaryGraphCheck3d<GEOSHAPE,pMapItemShare,pMapItemShare>
{
  public:
    typedef geoElement<GEOSHAPE>                           GEOELEMENT;
    typedef pGraph<GEOELEMENT,pMapItemShare,pMapItemShare> GRAPH;
    
  public:
    auxiliaryGraphCheck3d();
    bool check(const GRAPH & Graph);
};

template<typename GEOSHAPE>
auxiliaryGraphCheck3d<GEOSHAPE,pMapItemShare,pMapItemShare>::
auxiliaryGraphCheck3d()
{
}

template<typename GEOSHAPE>
bool
auxiliaryGraphCheck3d<GEOSHAPE,pMapItemShare,pMapItemShare>::
check(const GRAPH & Graph)
{
  bool flag = true;
  bool nodesSharing;
  UInt cid;
  
  for(UInt i=1; i <= Graph.rowSize(); ++i)
  {
    //Control that all the nodes conncected are shared 
    nodesSharing = true;
      
    for(UInt j=1 ; j <= Graph.rowSizeL(i); ++j)
    {
      cid  = Graph.getCid_LL(i,j);
      nodesSharing = nodesSharing && Graph.getColMap().get(cid).getShared();
    }
    
    //Is shared
    if(Graph.getRowMapL(i).getShared())
    {      
      flag = flag && nodesSharing;
    }
  }
  
  if(!flag)
  {cout << "A shared element should have all the nodes shared" << endl;}
  
  return(flag);
}


//_________________________________________________________________________________________________
// AUXILIARY CHECKING FUNCTIONS 2D
//-------------------------------------------------------------------------------------------------

/*! Auxiliary graph checks 2d - unspecialized, empty */
template<typename GEOSHAPE, typename ROWMAP, typename COLMAP>
class auxiliaryGraphCheck2d
{ };


/*! Auxiliary graph checks 2d - specialization for \c pMapItem */
template<typename GEOSHAPE>
class auxiliaryGraphCheck2d<GEOSHAPE,pMapItemShare,pMapItem>
{
  public:
    typedef geoElement<GEOSHAPE>                      GEOELEMENT;
    typedef pGraph<GEOELEMENT,pMapItemShare,pMapItem> GRAPH;
    
  public:
    auxiliaryGraphCheck2d();
    bool check(const GRAPH & Graph);
};

template<typename GEOSHAPE>
auxiliaryGraphCheck2d<GEOSHAPE,pMapItemShare,pMapItem>::
auxiliaryGraphCheck2d()
{ }

template<typename GEOSHAPE>
bool
auxiliaryGraphCheck2d<GEOSHAPE,pMapItemShare,pMapItem>::
check(const GRAPH & Graph)
{
  assert(Graph.size() == Graph.size());
  return(true);
}


/*! Auxiliary graph checks 2d - specialization for \c pMapItemShare */
template<typename GEOSHAPE>
class auxiliaryGraphCheck2d<GEOSHAPE,pMapItemShare,pMapItemShare>
{
  public:
    typedef geoElement<GEOSHAPE>                           GEOELEMENT;
    typedef pGraph<GEOELEMENT,pMapItemShare,pMapItemShare> GRAPH;
    
  public:
    auxiliaryGraphCheck2d();
    bool check(const GRAPH & Graph);
};

template<typename GEOSHAPE>
auxiliaryGraphCheck2d<GEOSHAPE,pMapItemShare,pMapItemShare>::
auxiliaryGraphCheck2d()
{
}

template<typename GEOSHAPE>
bool
auxiliaryGraphCheck2d<GEOSHAPE,pMapItemShare,pMapItemShare>::
check(const GRAPH & Graph)
{
  bool flag = true;
  bool nodesSharing;
  UInt cid;
  
  for(UInt i=1; i <= Graph.rowSize(); ++i)
  {
    //Control that all the nodes conncected are shared 
    nodesSharing = true;
      
    for(UInt j=1 ; j <= Graph.rowSizeL(i); ++j)
    {
      cid  = Graph.getCid_LL(i,j);
      nodesSharing = nodesSharing && Graph.getColMap().get(cid).getShared();
    }
    
    //Is shared
    if(Graph.getRowMapL(i).getShared())
    {      
      flag = flag && nodesSharing;
    }
  }
  
  if(!flag)
  {cout << "A shared element should have all the nodes shared" << endl;}
  
  return(flag);
}



//_________________________________________________________________________________________________
// AUXILIARY CHECKING FUNCTIONS 1D
//-------------------------------------------------------------------------------------------------

/*! Auxiliary graph checks 1d - unspecialized, empty */
template<typename GEOSHAPE, typename ROWMAP, typename COLMAP>
class auxiliaryGraphCheck1d
{ };


/*! Auxiliary graph checks 1d - specialization for \c pMapItem */
template<typename GEOSHAPE>
class auxiliaryGraphCheck1d<GEOSHAPE,pMapItemShare,pMapItem>
{
  public:
    typedef geoElement<GEOSHAPE>                      GEOELEMENT;
    typedef pGraph<GEOELEMENT,pMapItemShare,pMapItem> GRAPH;
    
  public:
    auxiliaryGraphCheck1d();
    bool check(const GRAPH & Graph);
};

template<typename GEOSHAPE>
auxiliaryGraphCheck1d<GEOSHAPE,pMapItemShare,pMapItem>::
auxiliaryGraphCheck1d()
{ }

template<typename GEOSHAPE>
bool
auxiliaryGraphCheck1d<GEOSHAPE,pMapItemShare,pMapItem>::
check(const GRAPH & Graph)
{
  assert(Graph.size() == Graph.size());
  return(true);
}


/*! Auxiliary graph checks 1d - specialization for \c pMapItemShare */
template<typename GEOSHAPE>
class auxiliaryGraphCheck1d<GEOSHAPE,pMapItemShare,pMapItemShare>
{
  public:
    typedef geoElement<GEOSHAPE>                           GEOELEMENT;
    typedef pGraph<GEOELEMENT,pMapItemShare,pMapItemShare> GRAPH;
    
  public:
    auxiliaryGraphCheck1d();
    bool check(const GRAPH & Graph);
};

template<typename GEOSHAPE>
auxiliaryGraphCheck1d<GEOSHAPE,pMapItemShare,pMapItemShare>::
auxiliaryGraphCheck1d()
{
}

template<typename GEOSHAPE>
bool
auxiliaryGraphCheck1d<GEOSHAPE,pMapItemShare,pMapItemShare>::
check(const GRAPH & Graph)
{
  bool flag = true;
  bool nodesSharing;
  UInt cid;
  
  for(UInt i=1; i <= Graph.rowSize(); ++i)
  {
    //Control that all the nodes conncected are shared 
    nodesSharing = true;
      
    for(UInt j=1 ; j <= Graph.rowSizeL(i); ++j)
    {
      cid  = Graph.getCid_LL(i,j);
      nodesSharing = nodesSharing && Graph.getColMap().get(cid).getShared();
    }
    
    //Is shared
    if(Graph.getRowMapL(i).getShared())
    {      
      flag = flag && nodesSharing;
    }
  }
  
  if(!flag)
  {cout << "A shared element should have all the nodes shared" << endl;}
  
  return(flag);
}



//_________________________________________________________________________________________________
// NODES FIXER FOR MESH DOCTORS
//-------------------------------------------------------------------------------------------------

/*! Nodes map fixer - \c pMapItem */
inline void nodesMapFixer(pMap<pMapItem> & map, const communicator & commDev)
{
  assert(map.size() == map.size());
  assert(commDev.rank() == commDev.rank());
}

/*! Nodes map fixer - \c pMapItemShare */
inline void nodesMapFixer(pMap<pMapItemShare> & map, const communicator & commDev)
{
  pMapGlobalManip<pMapItemShare> colMaxFixer(commDev);
  colMaxFixer.updateOwningSharing(map);
}


#endif
