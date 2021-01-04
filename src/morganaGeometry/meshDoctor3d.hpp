/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MESHDOCTOR3D_HPP
#define MESHDOCTOR3D_HPP

#include <set>
#include <boost/iterator/iterator_concepts.hpp>
#include "pMapGlobalManip.h"
#include "traitsGeometry.hpp"
#include "mesh3d.hpp"
#include "meshRefineUniform3d.hpp"

using namespace std;


/*! Mesh Doctor 3d */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class meshDoctor3d : public meshRefineUniform3d<GEOSHAPE,ELMAP,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef mesh3d<GEOSHAPE,ELMAP,NODEMAP>              MESH3D;
    typedef connect3d<GEOSHAPE,ELMAP,NODEMAP>           CONNECT3D;
    typedef meshRefineUniform3d<GEOSHAPE,ELMAP,NODEMAP> MESHREFINE;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 3;
    //@}
  
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    meshDoctor3d();
    meshDoctor3d(const Teuchos::RCP<const communicator> & CommDev);
    meshDoctor3d(const communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    //@}
    
    /*! @name Checks */ //@{
  public:
    bool checkGeoIds(const Teuchos::RCP<const MESH3D> & Mesh3d) const;
    bool checkGeoIds(const MESH3D & Mesh3d) const;
    bool checkJacobian(const Teuchos::RCP<const MESH3D> & Mesh3d) const;
    bool checkJacobian(const MESH3D & Mesh3d) const;
    //@}
    
    /*! @name Info */ //@{
  public:
    UInt countGeoIds(const Teuchos::RCP<const MESH3D> & Mesh3d) const;
    UInt countGeoIds(const MESH3D & Mesh3d) const;
    std::set<UInt> getGeoIds(const Teuchos::RCP<const MESH3D> & Mesh3d) const;
    std::set<UInt> getGeoIds(const MESH3D & Mesh3d) const;
    //@}
    
    /*! @name Operation functions */ //@{
  public:
    void removeUnusedPoints(const Teuchos::RCP<MESH3D> & Mesh3d) const;
    void removeUnusedPoints(MESH3D & Mesh3d) const;
    void fixGeoIds(Teuchos::RCP<MESH3D> & Mesh3d) const;
    void fixGeoIds(MESH3D & Mesh3d) const;
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
meshDoctor3d() : meshRefineUniform3d<GEOSHAPE,ELMAP,NODEMAP>()
{
  assert(GEOSHAPE::nDim == 3);
  commDevLoaded = false;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
meshDoctor3d(const Teuchos::RCP<const communicator> & CommDev) : meshRefineUniform3d<GEOSHAPE,ELMAP,NODEMAP>(CommDev)
{
  assert(GEOSHAPE::nDim == 3);
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
meshDoctor3d(const communicator & CommDev) : meshRefineUniform3d<GEOSHAPE,ELMAP,NODEMAP>(CommDev)
{
  assert(GEOSHAPE::nDim == 3);
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
  MESHREFINE::setCommunicator(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
  MESHREFINE::setCommunicator(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
bool
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
checkGeoIds(const Teuchos::RCP<const MESH3D> & Mesh3d) const
{
  assert(commDevLoaded);
  return(checkGeoIds(*Mesh3d));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
bool
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
checkGeoIds(const MESH3D & Mesh3d) const
{
  assert(commDevLoaded);
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> geoIdsSet;
  
  for(UInt i=1; i <= Mesh3d.getNumElements(); ++i)
  {
    geoIdsSet.insert(Mesh3d.getElementL(i).getGeoId());
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
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
checkJacobian(const Teuchos::RCP<const MESH3D> & Mesh3d) const
{
  assert(commDevLoaded);
  return(checkJacobian(*Mesh3d));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
bool
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
checkJacobian(const MESH3D & Mesh3d) const
{
  assert(commDevLoaded);
  
  //Local checking
  bool flag = true;
  sVect<point3d> points;
  point3d        refNode;
  
  for(UInt el=1; el <= Mesh3d.getNumElements(); ++el)
  {
    for(UInt ln=1; ln <= Mesh3d.getNumPoints(); ++ln)
    {
      points  = Mesh3d.getElementNodesL(el);
      refNode = Mesh3d.getRefNodes(ln);
      
      flag = flag & (Mesh3d.getGradientDet(points,refNode) > 0.0);
    }
  }
  
  //Global checking
  UInt flagTemp = flag;
  sVect<UInt> listTemp;
  
  all_gather(*commDev,flagTemp,listTemp);
  
  for(UInt i=1; i <= listTemp.size(); ++i)
  {
    flag = flag & bool(listTemp(i));
  }
  
  return(flag);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
countGeoIds(const Teuchos::RCP<const MESH3D> & Mesh3d) const
{
  assert(commDevLoaded);
  return(countGeoIds(*Mesh3d));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
UInt
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
countGeoIds(const MESH3D & Mesh3d) const
{
  assert(commDevLoaded);
  assert(checkGeoIds(Mesh3d));
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> geoIdsSet;
  
  for(UInt i=1; i <= Mesh3d.getNumElements(); ++i)
  {
    geoIdsSet.insert(Mesh3d.getElementL(i).getGeoId());
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
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
getGeoIds(const Teuchos::RCP<const MESH3D> & Mesh3d) const
{
  assert(commDevLoaded);
  return(getGeoIds(*Mesh3d));
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
std::set<UInt>
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
getGeoIds(const MESH3D & Mesh3d) const
{
  //Asserts
  assert(commDevLoaded);
  assert(checkGeoIds(Mesh3d));
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> localGeoIds;
  
  for(UInt i=1; i <= Mesh3d.getNumElements(); ++i)
  { localGeoIds.insert(Mesh3d.getElementL(i).getGeoId()); }
  
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
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
removeUnusedPoints(const Teuchos::RCP<MESH3D> & Mesh3d) const
{
  assert(commDevLoaded);
  removeUnusedPoints(*Mesh3d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
removeUnusedPoints(MESH3D & Mesh3d) const
{
  //Checks
  assert(commDevLoaded);
  assert(Mesh3d.getMapTransferred());
  
  //Clearing
  Mesh3d.clearFaces();
  Mesh3d.clearEdges();
  
  //Allocations
  UInt pid     = commDev->rank();
  UInt numPids = commDev->size();
  
  geoElement<GEOSHAPE> geoElement3d(true);
  UInt id, tot=1;
  
  sVect<bool> active(Mesh3d.getNumNodes());
  sVect<UInt> newOrder(Mesh3d.getNumNodes());
  
  for(UInt i=1; i <= active.size(); ++i)
  { active(i) = false; }
  
  geoMapInterface<GEOSHAPE> interface3d;
  
  //Find active nodes
  for(UInt i=1; i <= Mesh3d.getNumElements(); ++i)
  {
    for(UInt j=1; j <= interface3d.getNumPoints(); ++j)
    {
      id         = Mesh3d.getElementL(i).getCid(j);
      active(id) = true;
    }
  }
  
  //Nodes renumbering and reduction
  typedef pVect<point3d,NODEMAP>  NODESVECT;
  NODESVECT tempNodes;
  NODEMAP   nodeMap;
  
  for(UInt i=1; i <= active.size(); ++i)
  {
    if(active(i))
    {
      newOrder(i) = tot;
      
      nodeMap = Mesh3d.getNodes().getMapL(i);
      nodeMap.setLid(tot);
      
      tempNodes.push_back(Mesh3d.getNodes().getDataL(i), nodeMap);
      tot++;
    }
  }
  
  tempNodes.updateFinder();
  Mesh3d.setNodes(tempNodes);
   
  //Elements correction
  for(UInt i=1; i <= Mesh3d.getNumElements(); ++i)
  {
    geoElement3d = Mesh3d.getElementL(i);
    
    for(UInt j=1; j <= interface3d.getNumPoints(); ++j)
    {
      id = Mesh3d.getElementL(i).getCid(j);
      geoElement3d.setCid(j,newOrder(id));
    }
    
    Mesh3d.setElementL(i,geoElement3d);
  }  
  
  //Nodes map global fix
  pMap<NODEMAP> tempMap = Mesh3d.getNodes().getMapRef();
  
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
  localMapManip.setIndexing(tempMap);
  
  //Ownership fixing
  nodesMapFixer(tempMap,*commDev);
  
  //Map substitution
  tempNodes = Mesh3d.getNodes();
  tempNodes.setMap(tempMap);
  tempNodes.updateFinder();
  
  Mesh3d.setNodes(tempNodes);
  Mesh3d.computeNumVertices();
  Mesh3d.transferMap();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
fixGeoIds(Teuchos::RCP<MESH3D> & Mesh3d) const
{
  fixGeoIds(*Mesh3d);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
meshDoctor3d<GEOSHAPE,ELMAP,NODEMAP>::
fixGeoIds(MESH3D & Mesh3d) const
{
  assert(commDevLoaded);
  
  //Local geoIds 
  typedef set<UInt>::iterator ITERATOR;
  set<UInt> geoIdsSet;
  
  for(UInt i=1; i <= Mesh3d.getNumElements(); ++i)
  { geoIdsSet.insert(Mesh3d.getElementL(i).getGeoId()); }
  
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
  typedef typename MESH3D::GEOELEMENT3D GEOELEMENT3D;
  
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
  GEOELEMENT3D element;
  
  for(UInt i=1; i <= Mesh3d.getNumElements(); ++i)
  {
    k = Mesh3d.getElementL(i).getGeoId();
    k = geoIdMap.find(k)->second;
    
    element = Mesh3d.getElementL(i);
    element.setGeoId(k);
    Mesh3d.setElementL(i,element);
  }
}

#endif
