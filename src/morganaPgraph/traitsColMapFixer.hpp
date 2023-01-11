/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef TRAITSCOLMAPFIXER_HPP
#define TRAITSCOLMAPFIXER_HPP

#include "Teuchos_RCP.hpp"
#include <Epetra_Map.h>

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "pVectGlobalManip.hpp"
#include "pGraphComm.hpp"
#include "pGraph.hpp"

using namespace std;
using namespace boost::mpi;


inline void colMapFixer_overlap(pMap<pMapItem> & map,
                                const communicator & commDev)
{
  assert(map.size() == map.size());
  assert(commDev.rank() == commDev.rank());
}

inline void colMapFixer_overlap(pMap<pMapItemShare> & map,
                                const communicator & commDev)
{
  pMapGlobalManip<pMapItemShare> colMaxFixer(commDev);
  colMaxFixer.updateSharing(map);
}

inline void colMapFixer_changeMap(pMap<pMapItem> & map,
                                  const communicator & commDev)
{
  assert(map.size() == map.size());
  assert(commDev.rank() == commDev.rank());
}

inline void colMapFixer_changeMap(pMap<pMapItemShare> & map,
                                  const communicator & commDev)
{
  pMapGlobalManip<pMapItemShare> colMaxFixer(commDev);
  colMaxFixer.updateOwningSharing(map);
}

inline void upgradeMap_share(const pMap<pMapItem>      & inMap,
                                   pMap<pMapItemShare> & outMap,
                             const communicator        & commDev)
{
  assert(inMap.size() == inMap.size());
  assert(commDev.rank() == commDev.rank());
  
  pMapItemShare item;
  outMap.resize(inMap.size());
  
  for(UInt i=1; i <= inMap.size(); ++i)
  {
    item.setLid(inMap(i).getLid());
    item.setGid(inMap(i).getGid());
    item.setPid(inMap(i).getPid());
    
    outMap(i) = item;
  }
  
  pMapGlobalManip<pMapItemShare> colMaxFixer(commDev);
  colMaxFixer.updateOwningSharing(outMap);
}

inline void upgradeMap_share(const pMap<pMapItemShare> & inMap,
                                   pMap<pMapItemShare> & outMap,
                             const communicator        & commDev)
{
  outMap = inMap;
}

#endif
