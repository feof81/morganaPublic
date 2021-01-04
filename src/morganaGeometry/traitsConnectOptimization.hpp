/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef TRAITSCONNECTOPTIMIZATION_HPP
#define TRAITSCONNECTOPTIMIZATION_HPP

#include "pMapItem.h"
#include "pMapItemShare.h"
#include "pMap.hpp"
#include "pMapGlobalManip.h"

#include "pGraph.hpp"

using namespace std;


//General - not implemented class__________________________________________________________________

/*! Class for the optimization of connections.
Determine which items are far away from the boundary among sub-meshes.
Unspecialized, empty.*/
template<typename ITEM, typename ROWMAP, typename COLMAP>
class connectOptimization
{ };


//Class for pMapItem_______________________________________________________________________________

/*! Class for the optimization of connections.
Determine which items are far away from the boundary among sub-meshes.
Specialization for \c pMapItem */
template<typename ITEM, typename ROWMAP>
class connectOptimization<ITEM,ROWMAP,pMapItem>
{
  public:
    typedef pGraph<ITEM,ROWMAP,pMapItem> PGRAPH;
  
  public:
    connectOptimization();
    sVect<bool> getIsLocal(const PGRAPH & Graph);
};

template<typename ITEM, typename ROWMAP>
connectOptimization<ITEM,ROWMAP,pMapItem>::
connectOptimization()
{ }

template<typename ITEM, typename ROWMAP>
sVect<bool>
connectOptimization<ITEM,ROWMAP,pMapItem>::
getIsLocal(const PGRAPH & Graph)
{
  sVect<bool> outVect(Graph.rowSize());
  
  for(UInt i=1; i <= outVect.size(); ++i)
  {
    outVect(i) = false;
  }
  
  return(outVect);
}


//Class for pMapItemShare__________________________________________________________________________


/*! Class for the optimization of connections.
Determine which items are far away from the boundary among sub-meshes.
Specialization for \c pMapItemShare */
template<typename ITEM, typename ROWMAP>
class connectOptimization<ITEM,ROWMAP,pMapItemShare>
{
  public:
    typedef pGraph<ITEM,ROWMAP,pMapItemShare> PGRAPH;
  
  public:
    connectOptimization();
    sVect<bool> getIsLocal(const PGRAPH & Graph);
};

template<typename ITEM, typename ROWMAP>
connectOptimization<ITEM,ROWMAP,pMapItemShare>::
connectOptimization()
{ }

template<typename ITEM, typename ROWMAP>
sVect<bool>
connectOptimization<ITEM,ROWMAP,pMapItemShare>::
getIsLocal(const PGRAPH & Graph)
{
  sVect<bool> outVect(Graph.rowSize());
  bool shared;
  UInt cid;
 
  for(UInt i=1; i <= outVect.size(); ++i)
  {
    shared = true;
    
    for(UInt j=1; j <= Graph.rowSizeL(i); ++j)
    {
      cid    = Graph.getCid_LL(i,j);
      shared = shared & Graph.getColMapL(cid).getShared();
    }
    
    outVect(i) = !shared;
  }
  
  return(outVect);
}

#endif
