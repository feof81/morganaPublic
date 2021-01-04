/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MESHDOCTOR1D_HPP
#define MESHDOCTOR1D_HPP

#include <set>
#include "pMapGlobalManip.h"
#include "mesh1d.hpp"

using namespace std;


/*! Mesh Doctor 1d */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class meshDoctor1d
{
  public:
    typedef mesh1d<GEOSHAPE,ELMAP,NODEMAP> MESH1D;

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 1;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    meshDoctor1d();
    meshDoctor1d(const Teuchos::RCP<communicator> & CommDev);
    meshDoctor1d(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Checks */ //@{
  public:
    bool checkGeoIds(const Teuchos::RCP<MESH1D> & Mesh1d) const;
    bool checkGeoIds(const MESH1D & Mesh1d) const;
    //@}
    
    /*! @name Info */ //@{
  public:
    UInt countGeoIds(const Teuchos::RCP<MESH1D> & Mesh1d) const;
    UInt countGeoIds(const MESH1D & Mesh1d) const;
    //@}
    
    /*! @name Operation functions */ //@{
  public:
    void removeUnusedPoints(const Teuchos::RCP<MESH1D> & Mesh1d) const;
    void removeUnusedPoints(MESH1D & Mesh1d) const;
    void fixGeoIds(Teuchos::RCP<MESH1D> & Mesh1d) const;
    void fixGeoIds(MESH1D & Mesh1d) const;
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP>::
meshDoctor1d()
{
  assert(GEOSHAPE::nDim == 1);
  commDevLoaded = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP>::
meshDoctor1d(const Teuchos::RCP<communicator> & CommDev)
{
  assert(GEOSHAPE::nDim == 1);
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP>::
meshDoctor1d(communicator & CommDev)
{
  assert(GEOSHAPE::nDim == 1);
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
bool
meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP>::
checkGeoIds(const Teuchos::RCP<MESH1D> & Mesh1d) const
{
  assert(commDevLoaded);
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> geoIdsSet;
  
  for(UInt i=1; i <= Mesh1d->getNumElements(); ++i)
  {
    geoIdsSet.insert(Mesh1d->getElementL(i).getGeoId());
  }
  
  //Map translation
  UInt k = 1;
  UInt pid = commDev->rank();
  
  ITERATOR iter;
  pMapItem item;
  pMap<pMapItem> idsMap(geoIdsSet.size());
  
  for(iter = geoIdsSet.begin(); iter != geoIdsSet.end(); iter++)
  {
    item.setLid(k);
    item.setGid(*iter);
    item.setPid(pid);
    
    idsMap(k) = item;
    ++k;
  }
  
  //Checking
  pMapGlobalManip<pMapItem> checker(commDev);
  
  return(checker.check(idsMap));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
bool
meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP>::
checkGeoIds(const MESH1D & Mesh1d) const
{
  assert(commDevLoaded);
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> geoIdsSet;
  
  for(UInt i=1; i <= Mesh1d.getNumElements(); ++i)
  {
    geoIdsSet.insert(Mesh1d.getElementL(i).getGeoId());
  }
  
  //Map translation
  UInt k = 1;
  UInt pid = commDev->rank();
  
  ITERATOR iter;
  pMapItem item;
  pMap<pMapItem> idsMap(geoIdsSet.size());
  
  for(iter = geoIdsSet.begin(); iter != geoIdsSet.end(); iter++)
  {
    item.setLid(k);
    item.setGid(*iter);
    item.setPid(pid);
    
    idsMap(k) = item;
    ++k;
  }
  
  //Checking
  pMapGlobalManip<pMapItem> checker(commDev);
  
  return(checker.check(idsMap));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP>::
countGeoIds(const Teuchos::RCP<MESH1D> & Mesh1d) const
{
  assert(commDevLoaded);
  assert(checkGeoIds(Mesh1d));
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> geoIdsSet;
  
  for(UInt i=1; i <= Mesh1d->getNumElements(); ++i)
  {
    geoIdsSet.insert(Mesh1d->getElementL(i).getGeoId());
  }
  
  //Map translation
  UInt k = 1;
  UInt pid = commDev->rank();
  
  ITERATOR iter;
  pMapItem item;
  pMap<pMapItem> idsMap(geoIdsSet.size());
  
  for(iter = geoIdsSet.begin(); iter != geoIdsSet.end(); iter++)
  {
    item.setLid(k);
    item.setGid(*iter);
    item.setPid(pid);
    
    idsMap(k) = item;
    ++k;
  }
  
  //Checking
  pMapGlobalManip<pMapItem> checker(commDev);
  
  return(checker.sizeG(idsMap));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP>::
countGeoIds(const MESH1D & Mesh1d) const
{
  assert(commDevLoaded);
  assert(checkGeoIds(Mesh1d));
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> geoIdsSet;
  
  for(UInt i=1; i <= Mesh1d.getNumElements(); ++i)
  {
    geoIdsSet.insert(Mesh1d.getElementL(i).getGeoId());
  }
  
  //Map translation
  UInt k = 1;
  UInt pid = commDev->rank();
  
  ITERATOR iter;
  pMapItem item;
  pMap<pMapItem> idsMap(geoIdsSet.size());
  
  for(iter = geoIdsSet.begin(); iter != geoIdsSet.end(); iter++)
  {
    item.setLid(k);
    item.setGid(*iter);
    item.setPid(pid);
    
    idsMap(k) = item;
    ++k;
  }
  
  //Checking
  pMapGlobalManip<pMapItem> checker(commDev);
  
  return(checker.sizeG(idsMap));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP>::
removeUnusedPoints(const Teuchos::RCP<MESH1D> & Mesh1d) const
{
  //Checks
  assert(commDevLoaded);
  assert(Mesh1d->getMapTransferred());
  
  //Allocations
  UInt pid     = commDev->rank();
  UInt numPids = commDev->size();
  
  geoElement<GEOSHAPE> geoElement1d(true);
  UInt id;
  
  sVect<bool> active(Mesh1d->getNumNodes());
  sVect<UInt> newOrder(Mesh1d->getNumNodes());
  
  for(UInt i=1; i <= active.size(); ++i)
  { active(i) = false; }
  
  geoMapInterface<GEOSHAPE> interface1d;
  
  //Find active nodes
  for(UInt i=1; i <= Mesh1d->getNumElements(); ++i)
  {
    assert(interface1d.getNumPoints() == Mesh1d->getElementL(i).size());
    
    for(UInt j=1; j <= interface1d.getNumPoints(); ++j)
    {
      id         = Mesh1d->getElementL(i).getCid(j);
      active(id) = true;
    }
  }
  
  //Nodes renumbering and reduction
  typedef pVect<point3d,NODEMAP>  NODESVECT;
  NODESVECT tempNodes;
  NODEMAP   nodeMap;
  
  UInt tot = 1;
  
  for(UInt i=1; i <= active.size(); ++i)
  {
    if(active(i))
    {
      newOrder(i) = tot;
      
      nodeMap = Mesh1d->getNodes().getMapL(i);
      nodeMap.setLid(tot);
      
      tempNodes.push_back(Mesh1d->getNodes().getDataL(i), nodeMap);
      tot++;
    }
  }
  
  tempNodes.updateFinder();
  Mesh1d->setNodes(tempNodes);
   
  //Elements correction
  for(UInt i=1; i <= Mesh1d->getNumElements(); ++i)
  {
    geoElement1d = Mesh1d->getElementL(i);
    
    for(UInt j=1; j <= interface1d.getNumPoints(); ++j)
    {
      id = Mesh1d->getElementL(i).getCid(j);
      geoElement1d.setCid(j,newOrder(id));
    }
    
    Mesh1d->setElementL(i,geoElement1d);
  }
  
  
  //Nodes map global fix---------------------------------------------------------------------------
  pMap<NODEMAP> tempMap = Mesh1d->getNodes().getMapRef();
  tempMap.bufferLids();
  
  for(UInt i=1; i <= tempMap.size(); ++i)
  { tempMap(i).setPid(pid); }
  
  pMapGlobalManip<NODEMAP> nodesMapManip(commDev);
  UInt maxGid = nodesMapManip.sizeG(tempMap);
  
  pMapComm<NODEMAP> nodesComm(commDev);
  nodesComm.vectorNormal(tempMap,maxGid);
  
  //Temp map ordering
  typedef typename multiset<NODEMAP>::iterator ITERATOR;
  
  ITERATOR          iter;
  multiset<NODEMAP> tempSet;
  
  for(UInt i=1; i <= tempMap.size(); ++i)
  { tempSet.insert(tempMap(i)); }
  
  tot = 1;
  for(iter = tempSet.begin(); iter != tempSet.end(); iter++)
  {
    tempMap(tot) = *iter;
    tot++;
  }
  
  //Temp contiguous numbering
  tot = 1;
  newOrder.resize(tempMap.size());
  
  for(UInt i=1; i <= tempMap.size(); ++i)
  {
    if(i != 1)
    {
      if(tempMap(i).getGid() != tempMap(i-1).getGid())
      { tot++; }
    }
    
    newOrder(i) = tot;
  }
  
  for(UInt i=1; i <= tempMap.size(); ++i)
  {
    tempMap(i).setGid(newOrder(i));
  }
  
  //Offset determination 
  sVect<UInt> localSizes(numPids);
  all_gather(*commDev,tot,localSizes);
  
  UInt offset = 0;
  for(UInt i=0; i < pid; ++i)
  {
    offset += localSizes[i];
  }
  
  //Global numbering
  for(UInt i=1; i <= tempMap.size(); ++i)
  {
    id = tempMap(i).getGid();
    tempMap(i).setGid(id + offset);
  }
  
  //Pid map communication
  nodesComm.vectorPid(tempMap);
  
  //Local pid reordering
  pMapManip<NODEMAP> localMapManip;
  
  tempMap.restoreLids();
  localMapManip.setIndexing(tempMap);
  
  //Ownership fixing
  nodesMapFixer(tempMap,*commDev);
  
  //Map substitution
  tempNodes = Mesh1d->getNodes();
  tempNodes.setMap(tempMap);
  tempNodes.updateFinder();
  
  Mesh1d->setNodes(tempNodes);
  Mesh1d->computeNumVertices();
  Mesh1d->transferMap();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP>::
removeUnusedPoints(MESH1D & Mesh1d) const
{
  //Checks
  assert(commDevLoaded);
  assert(Mesh1d.getMapTransferred());
  
  //Allocations
  UInt pid     = commDev->rank();
  UInt numPids = commDev->size();
  
  geoElement<GEOSHAPE> geoElement1d(true);
  UInt id;
  
  sVect<bool> active(Mesh1d.getNumNodes());
  sVect<UInt> newOrder(Mesh1d.getNumNodes());
  
  for(UInt i=1; i <= active.size(); ++i)
  { active(i) = false; }
  
  geoMapInterface<GEOSHAPE> interface1d;
  
  //Find active nodes
  for(UInt i=1; i <= Mesh1d.getNumElements(); ++i)
  {
    for(UInt j=1; j <= interface1d.getNumPoints(); ++j)
    {
      assert(interface1d.getNumPoints() == Mesh1d.getElementL(i).size());
      
      id         = Mesh1d.getElementL(i).getCid(j);
      active(id) = true;
    }
  }
  
  //Nodes renumbering and reduction
  typedef pVect<point3d,NODEMAP>  NODESVECT;
  NODESVECT tempNodes;
  NODEMAP   nodeMap;
  
  UInt tot = 1;
  
  for(UInt i=1; i <= active.size(); ++i)
  {
    if(active(i))
    {
      newOrder(i) = tot;
      
      nodeMap = Mesh1d.getNodes().getMapL(i);
      nodeMap.setLid(tot);
      
      tempNodes.push_back(Mesh1d.getNodes().getDataL(i), nodeMap);
      tot++;
    }
  }
  
  tempNodes.updateFinder();
  Mesh1d.setNodes(tempNodes);
   
  //Elements correction
  for(UInt i=1; i <= Mesh1d.getNumElements(); ++i)
  {
    geoElement1d = Mesh1d.getElementL(i);
    
    for(UInt j=1; j <= interface1d.getNumPoints(); ++j)
    {
      id = Mesh1d.getElementL(i).getCid(j);
      geoElement1d.setCid(j,newOrder(id));
    }
    
    Mesh1d.setElementL(i,geoElement1d);
  }  
  
  
  //Nodes map global fix---------------------------------------------------------------------------
  pMap<NODEMAP> tempMap = Mesh1d.getNodes().getMapRef();
  tempMap.bufferLids();
  
  for(UInt i=1; i <= tempMap.size(); ++i)
  { tempMap(i).setPid(pid); }
  
  pMapGlobalManip<NODEMAP> nodesMapManip(commDev);
  UInt maxGid = nodesMapManip.sizeG(tempMap);
  
  pMapComm<NODEMAP> nodesComm(commDev);
  nodesComm.vectorNormal(tempMap,maxGid);
  
  //Temp map ordering
  typedef typename multiset<NODEMAP>::iterator ITERATOR;
  
  ITERATOR          iter;
  multiset<NODEMAP> tempSet;
  
  for(UInt i=1; i <= tempMap.size(); ++i)
  { tempSet.insert(tempMap(i)); }
  
  tot = 1;
  for(iter = tempSet.begin(); iter != tempSet.end(); iter++)
  {
    tempMap(tot) = *iter;
    tot++;
  }
  
  //Temp contiguous numbering
  tot = 1;
  newOrder.resize(tempMap.size());
  
  for(UInt i=1; i <= tempMap.size(); ++i)
  {
    if(i != 1)
    {
      if(tempMap(i).getGid() != tempMap(i-1).getGid())
      { tot++; }
    }
    
    newOrder(i) = tot;
  }
  
  for(UInt i=1; i <= tempMap.size(); ++i)
  {
    tempMap(i).setGid(newOrder(i));
  }
  
  //Offset determination 
  sVect<UInt> localSizes(numPids);
  all_gather(*commDev,tot,localSizes);
  
  UInt offset = 0;
  for(UInt i=0; i < pid; ++i)
  {
    offset += localSizes[i];
  }
  
  //Global numbering
  for(UInt i=1; i <= tempMap.size(); ++i)
  {
    id = tempMap(i).getGid();
    tempMap(i).setGid(id + offset);
  }
  
  //Pid map communication
  nodesComm.vectorPid(tempMap);
  
  //Local pid reordering
  pMapManip<NODEMAP> localMapManip;
  
  tempMap.restoreLids();
  localMapManip.setIndexing(tempMap);
  
  //Ownership fixing
  nodesMapFixer(tempMap,*commDev);
  
  //Map substitution
  tempNodes = Mesh1d.getNodes();
  tempNodes.setMap(tempMap);
  tempNodes.updateFinder();
  
  Mesh1d.setNodes(tempNodes);
  Mesh1d.computeNumVertices();
  Mesh1d.transferMap();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP>::
fixGeoIds(Teuchos::RCP<MESH1D> & Mesh1d) const
{
  fixGeoIds(*Mesh1d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor1d<GEOSHAPE,ELMAP,NODEMAP>::
fixGeoIds(MESH1D & Mesh1d) const
{
  assert(commDevLoaded);
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> geoIdsSet;
  
  for(UInt i=1; i <= Mesh1d.getNumElements(); ++i)
  { geoIdsSet.insert(Mesh1d.getElementL(i).getGeoId()); }
  
  ITERATOR iter;
  sVect<UInt> serialIds;
  
  for(iter = geoIdsSet.begin(); iter != geoIdsSet.end(); iter++)
  { serialIds.push_back(*iter); }
  
  //All to all communication
  sVect< sVect<UInt> > inValues;
  all_gather(*commDev, serialIds, inValues);
  
  //Final list
  for(UInt i=1; i <= inValues.size(); ++i)
  {
    for(UInt j=1; j <= inValues(i).size(); ++j)
    {
      geoIdsSet.insert(inValues(i)(j));
    }
  }
  
  //Map conversion
  typedef std::map<char,int>::iterator  MAPITER;
  typedef typename MESH1D::GEOELEMENT1D GEOELEMENT1D;
  
  map<UInt,UInt> geoIdMap;
  pair<UInt,UInt> paio;
  UInt k=1;
  
  for(iter = geoIdsSet.begin(); iter != geoIdsSet.end(); iter++)
  {
    paio.first  = *iter;
    paio.second = k;
    
    geoIdMap.insert(paio);
    ++k;
  }
  
  //Change geoIds
  GEOELEMENT1D element;
  
  for(UInt i=1; i <= Mesh1d.getNumElements(); ++i)
  {
    k = Mesh1d.getElementL(i).getGeoId();
    k = geoIdMap.find(k)->second;
    
    element = Mesh1d.getElementL(i);
    element.setGeoId(k);
    Mesh1d.setElementL(i,element);
  }
}


#endif
