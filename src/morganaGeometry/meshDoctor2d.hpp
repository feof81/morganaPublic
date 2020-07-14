/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MESHDOCTOR2D_HPP
#define MESHDOCTOR2D_HPP

#include <set>
#include "pMapGlobalManip.h"
#include "mesh2d.hpp"
#include "mesh3d.hpp"
#include "meshRefineUniform2d.hpp"

using namespace std;


/*! Mesh Doctor 2d */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class meshDoctor2d : public meshRefineUniform2d<GEOSHAPE,ELMAP,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef mesh2d<GEOSHAPE,ELMAP,NODEMAP>              MESH2D;
    typedef connect2d<GEOSHAPE,ELMAP,NODEMAP>           CONNECT3D;
    typedef meshRefineUniform2d<GEOSHAPE,ELMAP,NODEMAP> MESHREFINE;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 2;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    meshDoctor2d();
    meshDoctor2d(const Teuchos::RCP<communicator> & CommDev);
    meshDoctor2d(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Checks */ //@{
  public:
    bool checkGeoIds(const Teuchos::RCP<MESH2D> & Mesh2d) const;
    bool checkGeoIds(const MESH2D & Mesh2d) const;
    //@}
    
    /*! @name Info */ //@{
  public:
    UInt countGeoIds(const Teuchos::RCP<MESH2D> & Mesh2d) const;
    UInt countGeoIds(const MESH2D & Mesh2d) const;
    std::set<UInt> getGeoIds(const Teuchos::RCP<const MESH2D> & Mesh2d) const;
    std::set<UInt> getGeoIds(const MESH2D & Mesh2d) const;
    //@}
    
    /*! @name Operation functions */ //@{
  public:
    void removeUnusedPoints(const Teuchos::RCP<MESH2D> & Mesh2d) const;
    void removeUnusedPoints(MESH2D & Mesh2d) const;
    void fixGeoIds(Teuchos::RCP<MESH2D> & Mesh2d) const;
    void fixGeoIds(MESH2D & Mesh2d) const;
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
meshDoctor2d() : meshRefineUniform2d<GEOSHAPE,ELMAP,NODEMAP>()
{
  assert(GEOSHAPE::nDim == 2);
  commDevLoaded = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
meshDoctor2d(const Teuchos::RCP<communicator> & CommDev) : meshRefineUniform2d<GEOSHAPE,ELMAP,NODEMAP>(CommDev)
{
  assert(GEOSHAPE::nDim == 2);
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
meshDoctor2d(communicator & CommDev) : meshRefineUniform2d<GEOSHAPE,ELMAP,NODEMAP>(CommDev)
{
  assert(GEOSHAPE::nDim == 2);
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
  MESHREFINE::setCommunicator(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
  MESHREFINE::setCommunicator(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
bool
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
checkGeoIds(const Teuchos::RCP<MESH2D> & Mesh2d) const
{
  assert(commDevLoaded);
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> geoIdsSet;
  
  for(UInt i=1; i <= Mesh2d->getNumElements(); ++i)
  {
    geoIdsSet.insert(Mesh2d->getElementL(i).getGeoId());
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
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
checkGeoIds(const MESH2D & Mesh2d) const
{
  assert(commDevLoaded);
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> geoIdsSet;
  
  for(UInt i=1; i <= Mesh2d.getNumElements(); ++i)
  {
    geoIdsSet.insert(Mesh2d.getElementL(i).getGeoId());
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
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
countGeoIds(const Teuchos::RCP<MESH2D> & Mesh2d) const
{
  assert(commDevLoaded);
  assert(checkGeoIds(Mesh2d));
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> geoIdsSet;
  
  for(UInt i=1; i <= Mesh2d->getNumElements(); ++i)
  {
    geoIdsSet.insert(Mesh2d->getElementL(i).getGeoId());
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
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
countGeoIds(const MESH2D & Mesh2d) const
{
  assert(commDevLoaded);
  assert(checkGeoIds(Mesh2d));
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> geoIdsSet;
  
  for(UInt i=1; i <= Mesh2d.getNumElements(); ++i)
  {
    geoIdsSet.insert(Mesh2d.getElementL(i).getGeoId());
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
std::set<UInt>
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
getGeoIds(const Teuchos::RCP<const MESH2D> & Mesh2d) const
{
  assert(commDevLoaded);
  return(getGeoIds(*Mesh2d));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
std::set<UInt>
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
getGeoIds(const MESH2D & Mesh2d) const
{
  //Asserts
  assert(commDevLoaded);
  assert(checkGeoIds(Mesh2d));
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> localGeoIds;
  
  for(UInt i=1; i <= Mesh2d.getNumElements(); ++i)
  { localGeoIds.insert(Mesh2d.getElementL(i).getGeoId()); }
  
  //GlobalGeoIds
  sVect<UInt> inVect, outVect;
  set<UInt> globalGeoIds;
  
  for(ITERATOR iter = localGeoIds.begin(); iter != localGeoIds.end(); iter++)
  { inVect.push_back(*iter); }
  
  UInt * pointer = &inVect[0];
  UInt n = inVect.size();
  
  boost::mpi::all_gather(*commDev,pointer,n,outVect);
  
  for(UInt i=1; i <= outVect.size(); ++i)
  { globalGeoIds.insert(outVect(i)); }
  
  return(globalGeoIds);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
removeUnusedPoints(const Teuchos::RCP<MESH2D> & Mesh2d) const
{
  //Checks
  assert(commDevLoaded);
  assert(Mesh2d->getMapTransferred());
  
  //Clearing
  Mesh2d->clearEdges();
  
  //Allocations
  UInt pid     = commDev->rank();
  UInt numPids = commDev->size();
  
  geoElement<GEOSHAPE> geoElement3d(true);
  UInt id;
  
  sVect<bool> active(Mesh2d->getNumNodes());
  sVect<UInt> newOrder(Mesh2d->getNumNodes());
  
  for(UInt i=1; i <= active.size(); ++i)
  { active(i) = false; }
  
  geoMapInterface<GEOSHAPE> interface2d;
  
  //Find active nodes
  for(UInt i=1; i <= Mesh2d->getNumElements(); ++i)
  {
    assert(interface2d.getNumPoints() == Mesh2d->getElementL(i).size());
    
    for(UInt j=1; j <= interface2d.getNumPoints(); ++j)
    {
      id         = Mesh2d->getElementL(i).getCid(j);
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
      
      nodeMap = Mesh2d->getNodes().getMapL(i);
      nodeMap.setLid(tot);
      
      tempNodes.push_back(Mesh2d->getNodes().getDataL(i), nodeMap);
      tot++;
    }
  }
  
  tempNodes.updateFinder();
  Mesh2d->setNodes(tempNodes);
   
  //Elements correction
  for(UInt i=1; i <= Mesh2d->getNumElements(); ++i)
  {
    geoElement3d = Mesh2d->getElementL(i);
    
    for(UInt j=1; j <= interface2d.getNumPoints(); ++j)
    {
      id = Mesh2d->getElementL(i).getCid(j);
      geoElement3d.setCid(j,newOrder(id));
    }
    
    Mesh2d->setElementL(i,geoElement3d);
  }
  
  
  //Nodes map global fix---------------------------------------------------------------------------
  pMap<NODEMAP> tempMap = Mesh2d->getNodes().getMapRef();
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
  tempNodes = Mesh2d->getNodes();
  tempNodes.setMap(tempMap);
  tempNodes.updateFinder();
  
  Mesh2d->setNodes(tempNodes);
  Mesh2d->computeNumVertices();
  Mesh2d->transferMap();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
removeUnusedPoints(MESH2D & Mesh2d) const
{
  //Checks
  assert(commDevLoaded);
  assert(Mesh2d.getMapTransferred());
  
  //Clearing
  Mesh2d.clearEdges();
  
  //Allocations
  UInt pid     = commDev->rank();
  UInt numPids = commDev->size();
  
  geoElement<GEOSHAPE> geoElement3d(true);
  UInt id;
  
  sVect<bool> active(Mesh2d.getNumNodes());
  sVect<UInt> newOrder(Mesh2d.getNumNodes());
  
  for(UInt i=1; i <= active.size(); ++i)
  { active(i) = false; }
  
  geoMapInterface<GEOSHAPE> interface2d;
  
  //Find active nodes
  for(UInt i=1; i <= Mesh2d.getNumElements(); ++i)
  {
    for(UInt j=1; j <= interface2d.getNumPoints(); ++j)
    {
      assert(interface2d.getNumPoints() == Mesh2d.getElementL(i).size());
      
      id         = Mesh2d.getElementL(i).getCid(j);
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
      
      nodeMap = Mesh2d.getNodes().getMapL(i);
      nodeMap.setLid(tot);
      
      tempNodes.push_back(Mesh2d.getNodes().getDataL(i), nodeMap);
      tot++;
    }
  }
  
  tempNodes.updateFinder();
  Mesh2d.setNodes(tempNodes);
   
  //Elements correction
  for(UInt i=1; i <= Mesh2d.getNumElements(); ++i)
  {
    geoElement3d = Mesh2d.getElementL(i);
    
    for(UInt j=1; j <= interface2d.getNumPoints(); ++j)
    {
      id = Mesh2d.getElementL(i).getCid(j);
      geoElement3d.setCid(j,newOrder(id));
    }
    
    Mesh2d.setElementL(i,geoElement3d);
  }  
  
  
  //Nodes map global fix---------------------------------------------------------------------------
  pMap<NODEMAP> tempMap = Mesh2d.getNodes().getMapRef();
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
  tempNodes = Mesh2d.getNodes();
  tempNodes.setMap(tempMap);
  tempNodes.updateFinder();
  
  Mesh2d.setNodes(tempNodes);
  Mesh2d.computeNumVertices();
  Mesh2d.transferMap();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
fixGeoIds(Teuchos::RCP<MESH2D> & Mesh2d) const
{
  fixGeoIds(*Mesh2d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor2d<GEOSHAPE,ELMAP,NODEMAP>::
fixGeoIds(MESH2D & Mesh2d) const
{
  assert(commDevLoaded);
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> geoIdsSet;
  
  for(UInt i=1; i <= Mesh2d.getNumElements(); ++i)
  { geoIdsSet.insert(Mesh2d.getElementL(i).getGeoId()); }
  
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
  typedef std::map<char,int>::iterator MAPITER;
  typedef typename MESH2D::GEOELEMENT2D GEOELEMENT2D;
  
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
  GEOELEMENT2D element;
  
  for(UInt i=1; i <= Mesh2d.getNumElements(); ++i)
  {
    k = Mesh2d.getElementL(i).getGeoId();
    k = geoIdMap.find(k)->second;
    
    element = Mesh2d.getElementL(i);
    element.setGeoId(k);
    Mesh2d.setElementL(i,element);
  }
}


#endif

