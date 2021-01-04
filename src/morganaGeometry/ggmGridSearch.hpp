/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef GGMGRIDSEARCH_HPP
#define GGMGRIDSEARCH_HPP

#include "ggmGridSearchTraits.hpp"


template<typename TGT_MESH, typename SRC_MESH>
class ggmGridSearch : public ggmGridSearchTraits<TGT_MESH,SRC_MESH,SRC_MESH::nDim>
{
    /*! @name Typedefs */ //@{
  public:
    typedef ggmGridSearchTraits<TGT_MESH,SRC_MESH,SRC_MESH::nDim> TRAITSGS;
    typedef typename TRAITSGS::LOCALSEARCH    LOCALSEARCH;
    typedef typename LOCALSEARCH::OUTVECT     OUTVECT;
    typedef typename LOCALSEARCH::INVECT      INVECT;
    typedef typename LOCALSEARCH::INDATA      INDATA;
    typedef typename LOCALSEARCH::OUTDATA     OUTDATA;
    
    typedef typename TGT_MESH::GRID_GEOSHAPE  TGT_GEOSHAPE;
    typedef typename SRC_MESH::GRID_GEOSHAPE  SRC_GEOSHAPE;
    typedef typename LOCALSEARCH::TGT_NODEMAP TGT_NODEMAP;
    typedef typename LOCALSEARCH::SRC_NODEMAP SRC_NODEMAP;
    typedef typename LOCALSEARCH::TGT_ELMAP   TGT_ELMAP;
    typedef typename LOCALSEARCH::SRC_ELMAP   SRC_ELMAP;   
    //@}
  
    /*! @name Functions */ //@{
  public:
    ggmGridSearch();
    ggmGridSearch(const Teuchos::RCP<const communicator> & CommDev);
    ggmGridSearch(const communicator & CommDev);
    
    sVect<UInt> globalMatchingPids(const sVect<point3d> & points) const;
    OUTVECT search(const sVect<point3d> & locCoords);
    //@}
};


//_________________________________________________________________________________________________
// FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename TGT_MESH, typename SRC_MESH>
ggmGridSearch<TGT_MESH,SRC_MESH>::
ggmGridSearch() : ggmGridSearchTraits<TGT_MESH,SRC_MESH,SRC_MESH::nDim>::ggmGridSearchTraits()
{
}

template<typename TGT_MESH, typename SRC_MESH>
ggmGridSearch<TGT_MESH,SRC_MESH>::
ggmGridSearch(const Teuchos::RCP<const communicator> & CommDev) : ggmGridSearchTraits<TGT_MESH,SRC_MESH,SRC_MESH::nDim>::ggmGridSearchTraits(CommDev)
{
}

template<typename TGT_MESH, typename SRC_MESH>
ggmGridSearch<TGT_MESH,SRC_MESH>::
ggmGridSearch(const communicator & CommDev) : ggmGridSearchTraits<TGT_MESH,SRC_MESH,SRC_MESH::nDim>::ggmGridSearchTraits(CommDev)
{
}

template<typename TGT_MESH, typename SRC_MESH>
sVect<UInt>
ggmGridSearch<TGT_MESH,SRC_MESH>::
globalMatchingPids(const sVect<point3d> & points) const
{
  //Assert-----------------------------------------------------------
  assert(TRAITSGS::commDevOk);
  
  //Typedef----------------------------------------------------------
  typedef std::set<UInt>::iterator ITER;
  
  //Compute----------------------------------------------------------
  sVect<UInt> tempVect, outVect;
  std::set<UInt> outSet;
  
  for(UInt i=1; i <= points.size(); ++i)
  {
    tempVect = TRAITSGS::pointSearch.globalMatchingPids(points(i));
    
    for(UInt k=1; k <= tempVect.size(); ++k)
    { outSet.insert(tempVect(k)); }
  }
  
  for(ITER iter = outSet.begin(); iter != outSet.end(); ++iter)
  { outVect.push_back(*iter); }
  
  return(outVect);
}

template<typename TGT_MESH, typename SRC_MESH>
typename ggmGridSearch<TGT_MESH,SRC_MESH>::OUTVECT
ggmGridSearch<TGT_MESH,SRC_MESH>::
search(const sVect<point3d> & locCoords)
{
  //Typedefs---------------------------------------------------------------------------------------
  typedef searchData<SRC_ELMAP>  SEARCHPDATA;
  
  //Assert-----------------------------------------------------------------------------------------
  assert(TRAITSGS::commDevOk);
  
  //Local searching--------------------------------------------------------------------------------
  UInt pid = TRAITSGS::commDev->rank();
  
  INVECT inVect(TRAITSGS::tgtGrid->getNumElements());
  inVect.setMap(TRAITSGS::tgtGrid->getElements().getMapRef());
  
  for(UInt el=1; el <= TRAITSGS::tgtGrid->getNumElements(); ++el)
  { inVect(el) = INDATA(0,TRAITSGS::tgtGrid->getElementNodesL(el)); }
  
  inVect.updateFinder();
  
  OUTVECT outVect = LOCALSEARCH::localSearch(locCoords,inVect);
  
  //Global comm map creation-----------------------------------------------------------------------
  sVect<UInt> matchingPids;
  pMapItemSendRecv sendItem;
  pMap<pMapItemSendRecv> commMap;
  
  for(UInt i=1; i <= outVect.size(); ++i)
  {
    if(!outVect(i).getFound())
    {
      outVect(i).setIsNested(false);
      matchingPids = globalMatchingPids(TRAITSGS::tgtGrid->getElementNodesL(i));
      
      for(UInt j=1; j <= matchingPids.size(); ++j)
      {
        if(matchingPids(j) != pid)
        {
          assert(int(matchingPids(j)) < TRAITSGS::commDev->size());
  
          sendItem.setPid(pid);
          sendItem.setLid(outVect.getMapL(i).getLid());
          sendItem.setGid(outVect.getMapL(i).getGid());
          sendItem.setSid(pid);
          sendItem.setRid(matchingPids(j));
  
          commMap.push_back(sendItem);
        }
      }
    }
  }
  
  //Communication----------------------------------------------------------------------------------
  INVECT vectOther;
  
  pVectComm<INDATA,TGT_ELMAP> otherComm(TRAITSGS::commDev);
  otherComm.sendRecv(commMap,inVect,vectOther);
  
  //Searching the othed pids nodes-----------------------------------------------------------------
  OUTVECT dataOther = LOCALSEARCH::localSearch(locCoords,vectOther);
  
  //Transmit back----------------------------------------------------------------------------------
  pVectComm<OUTDATA,TGT_ELMAP> otherDataComm(TRAITSGS::commDev);
  otherDataComm.vectorPid(dataOther);
  
  //Correct unfound elements-----------------------------------------------------------------------
  UInt gid;
  
  for(UInt i=1; i <= dataOther.size(); ++i)
  {
    gid = dataOther.getMapL(i).getGid();
    assert(outVect.isG(gid));

    if(!outVect.getG(gid).getFound())
    { outVect.getG(gid) = dataOther(i); }
  }
  
  //Search points if previous algorithm has failed -- Alloc----------------------------------------
  point3d P;
  TGT_ELMAP tgtElmap;
  
  pVect<point3d,TGT_ELMAP>     Pvect;
  pVect<SEARCHPDATA,TGT_ELMAP> pDataVect;
  
  UInt locId    = 1;
  UInt locNum   = 0;
  UInt locStray = 0;
  sVect<UInt> allLocNum(TRAITSGS::commDev->size());
  
  //Search points if previous algorithm has failed -- Count ---------------------------------------
  for(UInt i=1; i <= outVect.size(); ++i)
  {
    if(!outVect(i).getFound())
    { locNum += 1 + locCoords.size(); }
  }
  
  boost::mpi::all_gather(*TRAITSGS::commDev,locNum,allLocNum);
  
  for(UInt i=0; i < TRAITSGS::commDev->rank(); ++i)
  { locStray += allLocNum[i]; }
  
  
  //Search points if previous algorithm has failed -- Alloc ---------------------------------------
  for(UInt el=1; el <= outVect.size(); ++el)
  {
    if(!outVect(el).getFound())
    {
      tgtElmap = outVect.getMapL(el);
      
      //Barycenter 
      P  = LOCALSEARCH::tgtProjector.getPosition(LOCALSEARCH::tgtGrid->getElementNodesL(el),
                                                 TGT_GEOSHAPE::getBarycenter());
      
      tgtElmap.setLid(locId);
      tgtElmap.setGid(locId + locStray);
      locId++;
      
      Pvect.push_back(P,tgtElmap);
      
      //LocCoords
      for(UInt j=1; j <= locCoords.size(); ++j)
      {
        P  = LOCALSEARCH::tgtProjector.getPosition(LOCALSEARCH::tgtGrid->getElementNodesL(el),
                                                   locCoords(j));
      
        tgtElmap.setLid(locId);
        tgtElmap.setGid(locId + locStray);
        locId++;
      
        Pvect.push_back(P,tgtElmap);
      }
    }
  }
  Pvect.updateFinder();
  
  //Search points if previous algorithm has failed -- Find ----------------------------------------
  TRAITSGS::commDev->barrier();
  TRAITSGS::pointSearch.findGlobal(Pvect,pDataVect);
  
  //Search points if previous algorithm has failed -- Unroll --------------------------------------
  locId = 1;
  
  for(UInt el=1; el <= outVect.size(); ++el)
  {
    if(!outVect(el).getFound())
    {
      //Barycenter search
      outVect(el).setElMap(pDataVect(locId).getElMap());
      outVect(el).setFound(pDataVect(locId).getFound());
      locId++;
      
      //Other nodes search
      outVect(el).setFound(true);
      outVect(el).setIsNested(false);
      outVect(el).getPData().resize(locCoords.size());
      
      for(UInt j=1; j <= locCoords.size(); ++j)
      {
        outVect(el).getPData(j) = pDataVect(locId);
        locId++;
      }
    }
  }
  
  //Return-----------------------------------------------------------------------------------------
  return(outVect);
}

#endif
