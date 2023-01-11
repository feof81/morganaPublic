/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PGRAPHGLOBALMANIP_HPP
#define PGRAPHGLOBALMANIP_HPP

#include "Teuchos_RCP.hpp"
#include <Epetra_Map.h>

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "pVectGlobalManip.hpp"
#include "pGraphComm.hpp"
#include "pGraph.hpp"
#include "traitsColMapFixer.hpp"

using namespace std;
using namespace boost::mpi;



//_________________________________________________________________________________________________
// UNSPECIALIZED
//-------------------------------------------------------------------------------------------------
/*! Parallel graph manipulator - unspecialized, empty */
template<typename ITEM, typename ROWMAP, typename COLMAP> class pGraphGlobalManip
{ };



//_________________________________________________________________________________________________
// PMAPITEM
//-------------------------------------------------------------------------------------------------
/*! Parallel graph manipulator \c pMapItem specialization */
template<typename ITEM, typename COLMAP> class pGraphGlobalManip<ITEM,pMapItem,COLMAP>
{
  /*! @name Typedefs */ //@{
  public:
    typedef pGraph<ITEM,pMapItem,COLMAP>        PGRAPH;
    typedef pGraph<ITEM,COLMAP,pMapItem>        PINVGRAPH;
    typedef pVect<ITEM,pMapItem>                PVECT;
    typedef pMap<pMapItemSendRecv>              SENDRECV;
    typedef pMap<pMapItem>                      CONTAINER_ROWMAP;
    typedef pMap<COLMAP>                        CONTAINER_COLMAP;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
    
    
    /*! @name Constructors and set functions */ //@{
  public:
    pGraphGlobalManip();
    pGraphGlobalManip(const Teuchos::RCP<const communicator> & CommDev);
    pGraphGlobalManip(const communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    //@}
    
    
    /*! @name Data retriving functions */ //@{
  public:
    /*! Returns the maximum row \c gid on all the processors */
    UInt sizeRowG(const Teuchos::RCP<const PGRAPH> & Graph) const;
    
    /*! Returns the maximum row \c gid on all the processors */
    UInt sizeRowG(const PGRAPH & Graph) const;
    
    /*! Returns the maximum col \c gid on all the processors */
    UInt sizeColG(const Teuchos::RCP<const PGRAPH> & Graph) const;
    
    /*! Returns the maximum col \c gid on all the processors */
    UInt sizeColG(const PGRAPH & Graph) const;
    
    /*! Perform some checks on both the coloumns and the rows */
    bool check(const Teuchos::RCP<const PGRAPH> & Graph) const;
    
    /*! Perform some checks on both the coloumns and the rows */
    bool check(const PGRAPH & Graph) const;
    //@}
    
    
    /*! @name Global manipulations */ //@{
  public:
    /*! Builds the global row numbering of the graph */
    void buildGlobalNumbering(Teuchos::RCP<PGRAPH> & Graph) const;
    
    /*! Builds the global row numbering of the graph */
    void buildGlobalNumbering(PGRAPH & Graph) const;
    
    /*! A higher performance method of \c buildGlobalNumbering. The local elements are marked with 
     the \c isLocal vector and they are not broadcasted*/
    void buildGlobalNumbering(Teuchos::RCP<PGRAPH> & Graph, const Teuchos::RCP<const sVect<bool> > & isRowLocal) const;
    
    /*! A higher performance method of \c buildGlobalNumbering. The local elements are marked with 
     the \c isLocal vector and they are not broadcasted*/
    void buildGlobalNumbering(PGRAPH & Graph, const sVect<bool> & isRowLocal) const;
    
    /*! Global inversion of the graph */
    void inversion(const Teuchos::RCP<const PGRAPH> & originGraph, Teuchos::RCP<PINVGRAPH> & targetGraph) const;
    
    /*! Global inversion of the graph */
    void inversion(const PGRAPH & originGraph, PINVGRAPH & targetGraph) const;
    
    /*! Enlarges the graph adding the subdomain neighbouring elements. The row and col numbering is retained while the neighbouring elements being added 
    are renumbered. */
    void overlap(Teuchos::RCP<PGRAPH> & Graph) const;
    
    /*! Enlarges the graph adding the subdomain neighbouring elements. The row and col numbering is retained while the neighbouring elements being added 
    are renumbered. */
    void overlap(PGRAPH & Graph) const;
    
    /*! Given a new map \c NewMap moves accordingly the elements across the communicator. The \c colMap s are updated accordingly */
    void changeRowMap(Teuchos::RCP<PGRAPH> & Graph, const Teuchos::RCP<const CONTAINER_ROWMAP> & NewMap) const;
    
    /*! Given a new map \c NewMap moves accordingly the elements across the communicator. The \c colMap s are updated accordingly */
    void changeRowMap(PGRAPH & Graph, const CONTAINER_ROWMAP & NewMap) const;
    //@}
    
    
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const PGRAPH>       & OldPgraph,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<PGRAPH>             & NewPgraph);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const PGRAPH       & OldPgraph,
                            const communicator & NewCommDev,
                                  PGRAPH       & NewPgraph);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const PGRAPH>       & OldPgraph,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<PGRAPH>             & NewPgraph);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const PGRAPH       & OldPgraph,
                            const communicator & NewCommDev,
                                  PGRAPH       & NewPgraph);
    //@}
    
    /*! @name Memory functions */ //@{
  public:
    size_t memSize() const;
    //@}
};


template<typename ITEM, typename COLMAP>
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
pGraphGlobalManip()
{
  commDevLoaded = false;
}

template<typename ITEM, typename COLMAP>
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
pGraphGlobalManip(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}
    
template<typename ITEM, typename COLMAP>
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
pGraphGlobalManip(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}
    
template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}
    
template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename ITEM, typename COLMAP>
UInt
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
sizeRowG(const Teuchos::RCP<const PGRAPH> & Graph) const
{
  pVectGlobalManip<ITEM,pMapItem> tempManip(commDev);
  return(tempManip.sizeG(Graph));
}

template<typename ITEM, typename COLMAP>
UInt
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
sizeRowG(const PGRAPH & Graph) const
{
  pVectGlobalManip<ITEM,pMapItem> tempManip(commDev);
  return(tempManip.sizeG(Graph));
}

template<typename ITEM, typename COLMAP>
UInt
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
sizeColG(const Teuchos::RCP<const PGRAPH> & Graph) const
{
  pMapManip<COLMAP> tempManip;
  return(tempManip.sizeG(Graph));
}

template<typename ITEM, typename COLMAP>
UInt
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
sizeColG(const PGRAPH & Graph) const
{
  pMapManip<COLMAP> tempManip;
  return(tempManip.sizeG(Graph));
}

template<typename ITEM, typename COLMAP>
bool
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
check(const Teuchos::RCP<const PGRAPH> & Graph) const
{
  assert(commDevLoaded);
  return(check(*Graph));
}

template<typename ITEM, typename COLMAP>
bool
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
check(const PGRAPH & Graph) const
{
  assert(commDevLoaded);
  
  //Pid
  UInt pid = commDev->rank();
  
  //Graph copy
  PGRAPH  graphCopy(Graph);
  
  //Row checking
  pVectGlobalManip<ITEM,pMapItem> rowChecker(commDev);
  bool flagRL = rowChecker.check(Graph);
  
  if(!flagRL)
  { cout << "ERROR: Graph row local check failed. Pid: " << pid << endl; }
  
  //Col checking
  pMapGlobalManip<COLMAP> colChecker(commDev);
  bool flagCG = colChecker.check(Graph.getColMap());
  
  if(!flagCG)
  { cout << "ERROR: Graph col global check failed. Pid: " << pid << endl; }
  
  //Local colChecking
  pGraphManip<ITEM,pMapItem,COLMAP> localCheker;
  bool flagCL = localCheker.checkCol(graphCopy);
  
  if(!flagCL)
  { cout << "ERROR: Graph column local check failed. Pid: " << pid << endl; }
  
  //Check that elements with the same gid have the same pattern
  pVectComm<ITEM,pMapItem> vectComm(commDev);
  
  graphCopy.pushToGlobal();
  vectComm.vectorNormal(graphCopy);
  
  pVectManip<ITEM,pMapItem> graphOrder;
  graphOrder.orderData(graphCopy);
  
  UInt gid = 0;
  bool flagSM = true;
  
  if(graphCopy.rowSize() > 0)
  { gid = graphCopy.getRowMapL(1).getGid(); }
  
  for(UInt i=2; i <= graphCopy.rowSize(); ++i)
  {
    if(gid == graphCopy.getRowMapL(i).getGid())
    {
      flagSM = flagSM & (!(graphCopy.getItemL(i) !=  graphCopy.getItemL(i-1)));
    }
    
     gid = graphCopy.getRowMapL(i).getGid();
  }
  
  if(!flagSM)
  { cout << "ERROR: Graph shared items not consistent. Pid: " << pid << endl; }
  
  //Flag bathering
  bool flag = flagRL & flagCG & flagCL & flagSM;
  UInt size = commDev->size();
  
  sVect<UInt> flags;
  all_gather(*commDev,UInt(flag),flags);
  
  for(UInt i=1; i <= size; ++i)
  { flag = flag & bool(flags(i));}
  
  return(flag);
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
buildGlobalNumbering(Teuchos::RCP<PGRAPH> & Graph) const
{
  assert(commDevLoaded);
  buildGlobalNumbering(*Graph);
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
buildGlobalNumbering(PGRAPH & Graph) const
{ 
  assert(commDevLoaded);
  
  //Graph copy
  PGRAPH  graphCopy(Graph);
  
  //Going global
  graphCopy.pushToGlobal();
  
  //Numbering
  pVectGlobalManip<ITEM,pMapItem> vectManip(commDev);
  vectManip.buildGlobalNumbering(graphCopy);
  
  //New map setting
  Graph.setRowMap(graphCopy.getRowMap());
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
buildGlobalNumbering(Teuchos::RCP<PGRAPH> & Graph, const Teuchos::RCP<const sVect<bool> > & isRowLocal) const
{ 
  assert(commDevLoaded);
  buildGlobalNumbering(*Graph,*isRowLocal);
}
    
template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
buildGlobalNumbering(PGRAPH & Graph, const sVect<bool> & isRowLocal) const
{ 
   assert(commDevLoaded);
  
  //Graph copy
  PGRAPH  graphCopy(Graph);
  
  //Going global
  graphCopy.pushToGlobal();
  
  //Numbering
  pVectGlobalManip<ITEM,pMapItem> vectManip(commDev);
  vectManip.buildGlobalNumbering(graphCopy,isRowLocal);
  
  //New map setting
  Graph.setRowMap(graphCopy.getRowMap());
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
inversion(const Teuchos::RCP<const PGRAPH> & originGraph, Teuchos::RCP<PINVGRAPH> & targetGraph) const
{
  assert(commDevLoaded);
  assert(check(originGraph));
  inversion(*originGraph,*targetGraph);
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
inversion(const PGRAPH & originGraph, PINVGRAPH & targetGraph) const
{
  assert(commDevLoaded);
  assert(check(originGraph));
  
  //Typedef
  typedef pair<pMapItem,ITEM>                   PAIR_ROW;
  typedef pair<COLMAP,ITEM>                     PAIR_COL;
  typedef typename map<pMapItem,ITEM>::iterator ITERATOR_ROW;
  typedef typename map<COLMAP,ITEM>::iterator   ITERATOR_COL;
  typedef pair<ITERATOR_ROW,bool>               BOOLPAIR_ROW;
  typedef pair<ITERATOR_COL,bool>               BOOLPAIR_COL;
  typedef typename set<pMapItem>::iterator      SETITER_ROW;
  typedef typename set<COLMAP>::iterator        SETITER_COL;
  typedef pair<SETITER_ROW,bool>                SETPAIR_ROW;
  typedef pair<SETITER_COL,bool>                SETPAIR_COL;
  
  //Clearing
  targetGraph.clear();
  
  //Allocations
  UInt proc = commDev->rank();
  
  //Local inversion
  PINVGRAPH graphInverse;
  pGraphManip<ITEM,pMapItem,COLMAP> localManip;
  localManip.inversion(originGraph,graphInverse);
   
  //Pid setup
  for(UInt i=1; i <= graphInverse.rowSize(); ++i)
  { graphInverse.getRowMapL(i).setPid(proc); }
  
  graphInverse.updateRowFinder();
  
  //Normal indexing
  PINVGRAPH tempInverse(graphInverse);
  tempInverse.pushToGlobal();  
  
  pVectComm<ITEM,COLMAP> invComm(commDev);
  invComm.vectorNormal(tempInverse);
  
  //Global col indices merging
  PAIR_COL         pair;
  ITERATOR_COL     iter;
  BOOLPAIR_COL     boolPair;
  map<COLMAP,ITEM> tempContainer;
  
  for(UInt i=1; i <= tempInverse.size(); ++i)
  {
    tempInverse.getRowMapL(i).setPid(proc);
    
    pair.first  = tempInverse.getRowMapL(i);
    pair.second = tempInverse.getItemL(i);
    
    boolPair = tempContainer.insert(pair);
    
    //Item already existing
    if( ! boolPair.second )
    {
      iter = boolPair.first;
      iter->second.merge(tempInverse.getItemL(i));
    }
  }
  
  //Download from stl-map
  tempInverse.clear();
  tempInverse.reserve(tempContainer.size());
  
  for(iter = tempContainer.begin(); iter != tempContainer.end(); ++iter)
  { tempInverse.push_back(iter->first, iter->second); }
  
  tempInverse.updateRowFinder();
  
  //Change map
  pVectGlobalManip<ITEM,COLMAP> mapSwitcher(commDev);
  mapSwitcher.changeMap(tempInverse,graphInverse.getRowMap());
  
  //Serialization of the colMap  
  SETPAIR_ROW    colPair;
  set<pMapItem>  colSet;
  
  for(UInt i=1; i <= graphInverse.getColMap().size(); ++i)
  {
    graphInverse.getColMap().get(i).setPid(proc);
    colSet.insert(graphInverse.getColMap().get(i));
  }
  
  //Substitution global-cid to local-cid
  pMap<pMapItem> newColMap(graphInverse.getColMap());
  pMapItem colItem;
  ITEM     item;
  
  for(UInt i=1; i <= tempInverse.rowSize(); ++i)
  {
    item = tempInverse.getItemL(i);
    
    for(UInt j=1; j <= item.size(); ++j)
    {
      colItem.setGid(item(j));
      colPair = colSet.insert(colItem);     
            
      //ColItem not already existing
      if(colPair.second)
      {
	//Buildup of the new item
	colItem.setLid(newColMap.size() + 1);
	colItem.setPid(proc);
	
	//Insert in the new col map
	newColMap.push_back(colItem);
      }
    }
  }
  
  //Setting pid, final graph build and push to local lids
  targetGraph = tempInverse;
  targetGraph.setColMap(newColMap);
  
  for(UInt i=1; i <= targetGraph.rowSize(); ++i)
  { targetGraph.getRowMapL(i).setPid(proc); }

  targetGraph.updateRowFinder();
  targetGraph.updateColFinder();  
  targetGraph.pushToLocal();
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
overlap(Teuchos::RCP<PGRAPH> & Graph) const
{
  assert(commDevLoaded);
  assert(check(Graph));
  overlap(*Graph);
} 
 
template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
overlap(PGRAPH & Graph) const
{
  assert(commDevLoaded);
  assert(check(Graph));
  
  //Typedefs
  typedef typename set<COLMAP>::iterator  SETITER_COl;
  typedef pair<SETITER_COl,bool>          SETPAIR_COL;
  
  //Allocations
  UInt proc = commDev->rank();
  
  //Graph copy
  PGRAPH graphCopy(Graph);
  
  //Graph inversion
  PINVGRAPH invGraph;
  inversion(Graph,invGraph);
  
  //Global elements overlapping
  graphCopy.pushToGlobal();
  
  pVectGlobalManip<ITEM,pMapItem> rowComm(commDev);
  rowComm.changeMap(graphCopy,invGraph.getColMap());
  
  //ColMapFix
  CONTAINER_COLMAP colMap(Graph.getColMap());
  
  //Serialization of the colMap  
  SETPAIR_COL  colPair;
  set<COLMAP>  colSet;
  
  for(UInt i=1; i <= colMap.size(); ++i)
  {
    colMap.get(i).setPid(proc);
    colSet.insert(colMap.get(i));
  }
  
  //Substitution global-cid to local-cid
  COLMAP   colItem;
  ITEM     item;
  
  for(UInt i=1; i <= graphCopy.rowSize(); ++i)
  {
    item = graphCopy.getItemL(i);
    
    for(UInt j=1; j <= item.size(); ++j)
    {
      colItem.setGid(item(j));
      colPair = colSet.insert(colItem);     
            
      //ColItem not already existing
      if(colPair.second)
      {
	//Buildup of the new item
	colItem.setLid(colMap.size() + 1);
	colItem.setPid(proc);
	
	//Insert in the new col map
	colMap.push_back(colItem);
      }
    }
  }
  
  //Setting pid, final graph build and push to local lids
  Graph = graphCopy;
  Graph.setColMap(colMap);
  
  for(UInt i=1; i <= Graph.rowSize(); ++i)
  { Graph.getRowMapL(i).setPid(proc); }

  Graph.updateRowFinder();
  Graph.updateColFinder();  
  Graph.pushToLocal();
  
  //Col fixing
  colMapFixer_overlap(Graph.getColMap(),*commDev);
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
changeRowMap(Teuchos::RCP<PGRAPH> & Graph, const Teuchos::RCP<const CONTAINER_ROWMAP> & NewMap) const
{
  //Pre-checking
  assert(commDevLoaded);
  changeRowMap(*Graph,*NewMap);
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
changeRowMap(PGRAPH & Graph, const CONTAINER_ROWMAP & NewMap) const
{
  //Pre-checking
  assert(commDevLoaded);
  pMapGlobalManip<pMapItem> preGlobal(commDev);
  
  assert(check(Graph));
  assert(preGlobal.check(NewMap));
  
  //Graph copy
  PGRAPH graphCopy(Graph);
  graphCopy.pushToGlobal();
  
  //ChangeMap - pVect
  pVectGlobalManip<ITEM,pMapItem> vectDistributor(commDev);
  vectDistributor.changeMap(graphCopy,NewMap);
  
  //New colMap
  typedef typename set<COLMAP>::iterator ITERATOR;
  
  ITEM item;
  COLMAP mapItem;
  ITERATOR iter;
  set<COLMAP> setMap;
  
  for(UInt i=1; i <= graphCopy.rowSize(); ++i)
  {
    item = graphCopy.getItemL(i);
    
    for(UInt j=1; j <= item.size(); ++j)
    { 
      mapItem.setGid(item(j));
      setMap.insert(mapItem);
    }
  }
  
  //New colMap serialization
  UInt k=1;
  pMap<COLMAP> newColMap(setMap.size());
  
  for(iter = setMap.begin(); iter != setMap.end(); ++iter)
  {
    mapItem = *iter;
    mapItem.setLid(k);
    
    newColMap(k) = mapItem;
    ++k;
  }
  
  //ColMap fixing
  colMapFixer_changeMap(newColMap,*commDev);
 
  //Col substitution
  graphCopy.setColMap(newColMap);
  graphCopy.pushToLocal();
  
  Graph = graphCopy;
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const PGRAPH>       & OldPgraph,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<PGRAPH>             & NewPgraph)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldPgraph,
                    *NewCommDev,
                    *NewPgraph);
}
    
template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
reduceCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const PGRAPH       & OldPgraph,
                   const communicator & NewCommDev,
                         PGRAPH       & NewPgraph)
{
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() < OldCommDev.size());
    
    //Alloc
    NewPgraph.clear();
    NewPgraph = OldPgraph;
    
    UInt pid = NewCommDev.rank();
    
    //Main loop
    for(UInt i=1; i <= NewPgraph.getRowMap().size(); ++i)
    { NewPgraph.getRowMap().get(i).setPid(pid); }
    
    for(UInt i=1; i <= NewPgraph.getColMap().size(); ++i)
    { NewPgraph.getColMap().get(i).setPid(pid); }
    
    NewPgraph.updateRowFinder();
    NewPgraph.updateColFinder();
  }
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const PGRAPH>       & OldPgraph,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<PGRAPH>             & NewPgraph)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldPgraph,
                    *NewCommDev,
                    *NewPgraph);
}
    
template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
expandCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const PGRAPH       & OldPgraph,
                   const communicator & NewCommDev,
                         PGRAPH       & NewPgraph)
{
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() > OldCommDev.size());
    
    //Alloc
    NewPgraph.clear();
    NewPgraph = OldPgraph;
    
    UInt pid = NewCommDev.rank();
    
    //Main loop
    for(UInt i=1; i <= NewPgraph.getRowMap().size(); ++i)
    { NewPgraph.getRowMap().get(i).setPid(pid); }
    
    for(UInt i=1; i <= NewPgraph.getColMap().size(); ++i)
    { NewPgraph.getColMap().get(i).setPid(pid); }
    
    NewPgraph.updateRowFinder();
    NewPgraph.updateColFinder();
  }
}

template<typename ITEM, typename COLMAP>
size_t
pGraphGlobalManip<ITEM,pMapItem,COLMAP>::
memSize() const
{
  return(sizeof(bool));
}


//_________________________________________________________________________________________________
// PMAPITEMSHARE
//-------------------------------------------------------------------------------------------------

/*! Parallel graph manipulator \c pMapItemShare specialization */
template<typename ITEM, typename COLMAP> class pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>
{
  /*! @name Typedefs */ //@{
  public:
    typedef pGraph<ITEM,pMapItemShare,COLMAP>    PGRAPH;
    typedef pGraph<ITEM,COLMAP,pMapItemShare>    PINVGRAPH;
    typedef pVect<ITEM,pMapItemShare>            PVECT;
    typedef pMap<pMapItemSendRecv>               SENDRECV;
    typedef pMap<pMapItemShare>                  CONTAINER_ROWMAP;
    typedef pMap<COLMAP>                         CONTAINER_COLMAP;
    //@}
    
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<const communicator> commDev;
    //@}
    
    
    /*! @name Constructors and set functions */ //@{
  public:
    pGraphGlobalManip();
    pGraphGlobalManip(const Teuchos::RCP<const communicator> & CommDev);
    pGraphGlobalManip(const communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    //@}
    
    
    /*! @name Data retriving functions */ //@{
  public:
    /*! Returns the maximum row \c gid on all the processors */
    UInt sizeRowG(const Teuchos::RCP<const PGRAPH> & Graph) const;
    
    /*! Returns the maximum row \c gid on all the processors */
    UInt sizeRowG(const PGRAPH & Graph) const;
    
    /*! Returns the maximum col \c gid on all the processors */
    UInt sizeColG(const Teuchos::RCP<const PGRAPH> & Graph) const;
    
    /*! Returns the maximum col \c gid on all the processors */
    UInt sizeColG(const PGRAPH & Graph) const;
    
    /*! Perform some checks on both the coloumns and the rows */
    bool check(const Teuchos::RCP<const PGRAPH> & Graph) const;
    
    /*! Perform some checks on both the coloumns and the rows */
    bool check(const PGRAPH & Graph) const;
    //@}
    
    
    /*! @name Global manipulations */ //@{
  public:
    /*! Builds the global row numbering of the graph and updates the ownership structure */
    void buildGlobalNumbering(Teuchos::RCP<PGRAPH> & Graph) const;
    
    /*! Builds the global row numbering of the graph and updates the ownership structure */
    void buildGlobalNumbering(PGRAPH & Graph) const;
    
    /*! A higher performance method of \c buildGlobalNumbering. The local elements are marked with 
     the \c isLocal vector and they are not broadcasted*/
    void buildGlobalNumbering(Teuchos::RCP<PGRAPH> & Graph, const Teuchos::RCP<const sVect<bool> > & isRowLocal) const;
    
    /*! A higher performance method of \c buildGlobalNumbering. The local elements are marked with 
     the \c isLocal vector and they are not broadcasted*/
    void buildGlobalNumbering(PGRAPH & Graph, const sVect<bool> & isRowLocal) const;
    
    /*! Global inversion of the graph. The \c targetGraph has an updated ownership of the colMap */
    void inversion(const Teuchos::RCP<const PGRAPH> & originGraph, Teuchos::RCP<PINVGRAPH> & targetGraph) const;
    
    /*! Global inversion of the graph. The \c targetGraph has an updated ownership of the colMap */
    void inversion(const PGRAPH & originGraph, PINVGRAPH & targetGraph) const;
    
    /*! Enlarges the graph adding the subdomain neighbouring elements. The row and col numbering is retained while the neighbouring elements being added 
    are renumbered. Both the ownership structure of the rows and coloums are updated */
    void overlap(Teuchos::RCP<PGRAPH> & Graph) const;
    
    /*! Enlarges the graph adding the subdomain neighbouring elements. The row and col numbering is retained while the neighbouring elements being added 
    are renumbered. Both the ownership structure of the rows and coloums are updated */
    void overlap(PGRAPH & Graph) const;
    
    /*! Given a new map \c NewMap moves accordingly the elements across the communicator. The \c colMap s are updated accordingly */
    void changeRowMap(Teuchos::RCP<PGRAPH> & Graph, const Teuchos::RCP<const CONTAINER_ROWMAP> & NewMap) const;
    
    /*! Given a new map \c NewMap moves accordingly the elements across the communicator. The \c colMap s are updated accordingly */
    void changeRowMap(PGRAPH & Graph, const CONTAINER_ROWMAP & NewMap) const;
    //@}
    
    
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const PGRAPH>       & OldPgraph,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<PGRAPH>             & NewPgraph);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const PGRAPH       & OldPgraph,
                            const communicator & NewCommDev,
                                  PGRAPH       & NewPgraph);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const PGRAPH>       & OldPgraph,
                            const Teuchos::RCP<const communicator> & NewCommDev,
                                  Teuchos::RCP<PGRAPH>             & NewPgraph);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool         & isActive,
                            const communicator & OldCommDev,
                            const PGRAPH       & OldPgraph,
                            const communicator & NewCommDev,
                                  PGRAPH       & NewPgraph);
    //@}
    
    /*! @name Memory functions */ //@{
  public:
    size_t memSize() const;
    //@}
};


template<typename ITEM, typename COLMAP>
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
pGraphGlobalManip()
{
  commDevLoaded = false;
}

template<typename ITEM, typename COLMAP>
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
pGraphGlobalManip(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}
    
template<typename ITEM, typename COLMAP>
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
pGraphGlobalManip(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}
    
template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}
    
template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
setCommunicator(const communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename ITEM, typename COLMAP>
UInt
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
sizeRowG(const Teuchos::RCP<const PGRAPH> & Graph) const
{
  pVectGlobalManip<ITEM,pMapItemShare> tempManip(commDev);
  return(tempManip.sizeG(Graph));
}

template<typename ITEM, typename COLMAP>
UInt
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
sizeRowG(const PGRAPH & Graph) const
{
  pVectGlobalManip<ITEM,pMapItemShare> tempManip(commDev);
  return(tempManip.sizeG(Graph));
}

template<typename ITEM, typename COLMAP>
UInt
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
sizeColG(const Teuchos::RCP<const PGRAPH> & Graph) const
{
  pMapManip<COLMAP> tempManip;
  return(tempManip.sizeG(Graph));
}

template<typename ITEM, typename COLMAP>
UInt
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
sizeColG(const PGRAPH & Graph) const
{
  pMapManip<COLMAP> tempManip;
  return(tempManip.sizeG(Graph));
}

template<typename ITEM, typename COLMAP>
bool
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
check(const Teuchos::RCP<const PGRAPH> & Graph) const
{
  assert(commDevLoaded);
  return(check(*Graph));
}

template<typename ITEM, typename COLMAP>
bool
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
check(const PGRAPH & Graph) const
{
  assert(commDevLoaded);
  
  //Pid
  UInt pid = commDev->rank();
  
  //Graph copy
  PGRAPH  graphCopy(Graph);
  
  //Row checking
  pVectGlobalManip<ITEM,pMapItemShare> rowChecker(commDev);
  bool flagRL = rowChecker.check(Graph);
  
  if(!flagRL)
  { cout << "ERROR: Graph row local check failed. Pid: " << pid << endl; }
  
  //Col checking
  pMapGlobalManip<COLMAP> colChecker(commDev);
  bool flagCG = colChecker.check(Graph.getColMap());
  
  if(!flagCG)
  { cout << "ERROR: Graph col global check failed. Pid: " << pid << endl; }
  
  //Local colChecking
  pGraphManip<ITEM,pMapItemShare,COLMAP> localCheker;
  bool flagCL = localCheker.checkCol(graphCopy);
  
  if(!flagCL)
  { cout << "ERROR: Graph column local check failed. Pid: " << pid << endl; }
  
  //Check that elements with the same gid have the same pattern
  pVectComm<ITEM,pMapItemShare> vectComm(commDev);
  
  graphCopy.pushToGlobal();
  vectComm.vectorNormal(graphCopy);
  
  pVectManip<ITEM,pMapItemShare> graphOrder;
  graphOrder.orderData(graphCopy);
  
  UInt gid = 0;
  bool flagSM = true;
  
  if(graphCopy.rowSize() > 0)
  { gid = graphCopy.getRowMapL(1).getGid(); }
  
  for(UInt i=2; i <= graphCopy.rowSize(); ++i)
  {
    if(gid == graphCopy.getRowMapL(i).getGid())
    {
      flagSM = flagSM & (!(graphCopy.getItemL(i) !=  graphCopy.getItemL(i-1)));
    }
    
    gid = graphCopy.getRowMapL(i).getGid();
  }
  
  if(!flagSM)
  { cout << "ERROR: Graph shared items not consistent. Pid: " << pid << endl; }
  
  //Flag bathering
  bool flag = flagRL & flagCG & flagCL & flagSM;
  UInt size = commDev->size();
  
  sVect<UInt> flags;
  all_gather(*commDev,UInt(flag),flags);
  
  for(UInt i=1; i <= size; ++i)
  { flag = flag & bool(flags(i));}
  
  return(flag);
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
buildGlobalNumbering(Teuchos::RCP<PGRAPH> & Graph) const
{
  assert(commDevLoaded);
  assert(Graph->getColMap().size() != 0);
  buildGlobalNumbering(*Graph);
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
buildGlobalNumbering(PGRAPH & Graph) const
{
  assert(commDevLoaded);
  assert(Graph.getColMap().size() != 0);
  
  //Graph copy
  PGRAPH  graphCopy(Graph);
  
  //Going global
  graphCopy.pushToGlobal();
  
  //Numbering
  pVectGlobalManip<ITEM,pMapItemShare> vectManip(commDev);
  vectManip.buildGlobalNumbering(graphCopy);
  
  //New map setting
  Graph.setRowMap(graphCopy.getRowMap());
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
buildGlobalNumbering(Teuchos::RCP<PGRAPH> & Graph, const Teuchos::RCP<const sVect<bool> > & isRowLocal) const
{
  assert(commDevLoaded);
  buildGlobalNumbering(*Graph,*isRowLocal);
}
    
template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
buildGlobalNumbering(PGRAPH & Graph, const sVect<bool> & isRowLocal) const
{
  assert(commDevLoaded);
  
  //Graph copy
  PGRAPH  graphCopy(Graph);
  
  //Going global
  graphCopy.pushToGlobal();
  
  //Numbering
  pVectGlobalManip<ITEM,pMapItemShare> vectManip(commDev);
  vectManip.buildGlobalNumbering(graphCopy,isRowLocal);
  
  //New map setting
  Graph.setRowMap(graphCopy.getRowMap());
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
inversion(const Teuchos::RCP<const PGRAPH> & originGraph, Teuchos::RCP<PINVGRAPH> & targetGraph) const
{
  assert(commDevLoaded);
  assert(check(originGraph));
  inversion(*originGraph,*targetGraph);
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
inversion(const PGRAPH & originGraph, PINVGRAPH & targetGraph) const
{
  assert(commDevLoaded);
  assert(check(originGraph));
  
  //Typedef
  typedef pair<pMapItemShare,ITEM>                    PAIR_ROW;
  typedef pair<COLMAP,ITEM>                           PAIR_COL;
  typedef typename map<pMapItemShare,ITEM>::iterator  ITERATOR_ROW;
  typedef typename map<COLMAP,ITEM>::iterator         ITERATOR_COL;
  typedef pair<ITERATOR_ROW,bool>                     BOOLPAIR_ROW;
  typedef pair<ITERATOR_COL,bool>                     BOOLPAIR_COL;
  typedef typename set<pMapItemShare>::iterator       SETITER_ROW;
  typedef typename set<COLMAP>::iterator              SETITER_COL;
  typedef pair<SETITER_ROW,bool>                      SETPAIR_ROW;
  typedef pair<SETITER_COL,bool>                      SETPAIR_COL;
  
  //Clearing
  targetGraph.clear();
  
  //Allocations
  UInt proc = commDev->rank();
  
  //Local inversion
  PINVGRAPH graphInverse;
  pGraphManip<ITEM,pMapItemShare,COLMAP> localManip;
  localManip.inversion(originGraph,graphInverse);
  
  //Pid setup
  for(UInt i=1; i <= graphInverse.rowSize(); ++i)
  { graphInverse.getRowMapL(i).setPid(proc); }
  
  graphInverse.updateRowFinder();
  
  //Normal indexing
  PINVGRAPH tempInverse(graphInverse);
  tempInverse.pushToGlobal();  
  
  pVectComm<ITEM,COLMAP> invComm(commDev);
  invComm.vectorNormal(tempInverse);
  
  //Global col indices merging
  PAIR_COL         pair;
  ITERATOR_COL     iter;
  BOOLPAIR_COL     boolPair;
  map<COLMAP,ITEM> tempContainer;
  
  for(UInt i=1; i <= tempInverse.size(); ++i)
  {
    tempInverse.getRowMapL(i).setPid(proc);
    
    pair.first  = tempInverse.getRowMapL(i);
    pair.second = tempInverse.getItemL(i);
    
    boolPair = tempContainer.insert(pair);
    
    //Item already existing
    if( ! boolPair.second )
    {
      iter = boolPair.first;
      iter->second.merge(tempInverse.getItemL(i));
    }
  }
  
  //Download from stl-map
  tempInverse.clear();
  tempInverse.reserve(tempContainer.size());
  
  for(iter = tempContainer.begin(); iter != tempContainer.end(); ++iter)
  { tempInverse.push_back(iter->first, iter->second); }
  
  tempInverse.updateRowFinder();
  
  //Change map
  pVectGlobalManip<ITEM,COLMAP> mapSwitcher(commDev);
  mapSwitcher.changeMap(tempInverse,graphInverse.getRowMap());
  
  //Serialization of the colMap  
  SETPAIR_ROW         colPair;
  set<pMapItemShare>  colSet;
  
  for(UInt i=1; i <= graphInverse.getColMap().size(); ++i)
  {
    graphInverse.getColMap().get(i).setPid(proc);
    colSet.insert(graphInverse.getColMap().get(i));
  }
  
  //Substitution global-cid to local-cid
  pMap<pMapItemShare> newColMap(graphInverse.getColMap());
  pMapItemShare colItem;
  ITEM     item;
  
  for(UInt i=1; i <= tempInverse.rowSize(); ++i)
  {
    item = tempInverse.getItemL(i);
    
    for(UInt j=1; j <= item.size(); ++j)
    {
      colItem.setGid(item(j));
      colPair = colSet.insert(colItem);     
            
      //ColItem not already existing
      if(colPair.second)
      {
	//Buildup of the new item
	colItem.setLid(newColMap.size() + 1);
	colItem.setPid(proc);
	colItem.setShared(true); //! variation 
	
	//Insert in the new col map
	newColMap.push_back(colItem);
      }
    }
  }
  
  //! Update the ownership of the colMap - variation
  pMapGlobalManip<pMapItemShare> colMapFixer(commDev);
  colMapFixer.updateOwningSharing(newColMap);
  
  //Setting pid, final graph build and push to local lids
  targetGraph = tempInverse;
  targetGraph.setColMap(newColMap);
  
  for(UInt i=1; i <= targetGraph.rowSize(); ++i)
  { targetGraph.getRowMapL(i).setPid(proc); }

  targetGraph.updateRowFinder();
  targetGraph.updateColFinder();  
  targetGraph.pushToLocal();
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
overlap(Teuchos::RCP<PGRAPH> & Graph) const
{
  assert(commDevLoaded);
  assert(check(Graph));
  overlap(*Graph);
}
    
template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
overlap(PGRAPH & Graph) const
{
  assert(commDevLoaded);
  assert(check(Graph));
  
  //Typedefs
  typedef typename set<COLMAP>::iterator  SETITER_COl;
  typedef pair<SETITER_COl,bool>          SETPAIR_COL;
  
  //Allocations
  UInt proc = commDev->rank();
  
  //Graph copy
  PGRAPH graphCopy(Graph);
  
  //Graph inversion
  PINVGRAPH invGraph;
  inversion(Graph,invGraph);
  
  //Global elements overlapping
  graphCopy.pushToGlobal();
  
  pVectGlobalManip<ITEM,pMapItemShare> rowComm(commDev);
  rowComm.changeMap(graphCopy,invGraph.getColMap());
  
  //ColMapFix
  CONTAINER_COLMAP colMap(Graph.getColMap());
  
  //Serialization of the colMap  
  SETPAIR_COL  colPair;
  set<COLMAP>  colSet;
  
  for(UInt i=1; i <= colMap.size(); ++i)
  {
    colMap.get(i).setPid(proc);
    colSet.insert(colMap.get(i));
  }
  
  //Substitution global-cid to local-cid
  COLMAP   colItem;
  ITEM     item;
  
  for(UInt i=1; i <= graphCopy.rowSize(); ++i)
  {
    item = graphCopy.getItemL(i);
    
    for(UInt j=1; j <= item.size(); ++j)
    {
      colItem.setGid(item(j));
      colPair = colSet.insert(colItem);     
            
      //ColItem not already existing
      if(colPair.second)
      {
	//Buildup of the new item
	colItem.setLid(colMap.size() + 1);
	colItem.setPid(proc);
	
	//Insert in the new col map
	colMap.push_back(colItem);
      }
    }
  }
  
  //Setting pid, final graph build and push to local lids
  Graph = graphCopy;
  Graph.setColMap(colMap);
  
  for(UInt i=1; i <= Graph.rowSize(); ++i)
  { Graph.getRowMapL(i).setPid(proc); }

  Graph.updateRowFinder();
  Graph.updateColFinder();  
  Graph.pushToLocal();
  
  //Col fixing
  colMapFixer_overlap(Graph.getColMap(),*commDev);
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
changeRowMap(Teuchos::RCP<PGRAPH> & Graph, const Teuchos::RCP<const CONTAINER_ROWMAP> & NewMap) const
{
  //Pre-checking
  assert(commDevLoaded);
  changeRowMap(*Graph,*NewMap);
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
changeRowMap(PGRAPH & Graph, const CONTAINER_ROWMAP & NewMap) const
{
  //Pre-checking
  assert(commDevLoaded);
  pMapGlobalManip<pMapItemShare> preGlobal(commDev);
  
  assert(check(Graph));
  assert(preGlobal.check(NewMap));
  
  //Graph copy
  PGRAPH graphCopy(Graph);
  graphCopy.pushToGlobal();
  
  //ChangeMap - pVect
  pVectGlobalManip<ITEM,pMapItemShare> vectDistributor(commDev);
  vectDistributor.changeMap(graphCopy,NewMap);
  
  //New colMap
  typedef typename set<COLMAP>::iterator ITERATOR;
  
  ITEM item;
  COLMAP mapItem;
  ITERATOR iter;
  set<COLMAP> setMap;
  
  for(UInt i=1; i <= graphCopy.rowSize(); ++i)
  {
    item = graphCopy.getItemL(i);
    
    for(UInt j=1; j <= item.size(); ++j)
    { 
      mapItem.setGid(item(j));
      setMap.insert(mapItem);
    }
  }
  
  //New colMap serialization
  UInt k=1;
  pMap<COLMAP> newColMap(setMap.size());
  
  for(iter = setMap.begin(); iter != setMap.end(); ++iter)
  {
    mapItem = *iter;
    mapItem.setLid(k);
    
    newColMap(k) = mapItem;
    ++k;
  }
  
  //ColMap fixing
  colMapFixer_changeMap(newColMap,*commDev);
 
  //Col substitution
  graphCopy.setColMap(newColMap);
  graphCopy.pushToLocal();
  
  Graph = graphCopy;
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const PGRAPH>       & OldPgraph,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<PGRAPH>             & NewPgraph)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldPgraph,
                    *NewCommDev,
                    *NewPgraph);
}
    
template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
reduceCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const PGRAPH       & OldPgraph,
                   const communicator & NewCommDev,
                         PGRAPH       & NewPgraph)
{
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() < OldCommDev.size());
    
    //Alloc
    NewPgraph.clear();
    NewPgraph = OldPgraph;
    
    UInt pid = NewCommDev.rank();
    
    //Main loop
    for(UInt i=1; i <= NewPgraph.getRowMap().size(); ++i)
    { NewPgraph.getRowMap().get(i).setPid(pid); }
    
    for(UInt i=1; i <= NewPgraph.getColMap().size(); ++i)
    { NewPgraph.getColMap().get(i).setPid(pid); }
    
    NewPgraph.updateRowFinder();
    NewPgraph.updateColFinder();
  }
}

template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const PGRAPH>       & OldPgraph,
                   const Teuchos::RCP<const communicator> & NewCommDev,
                         Teuchos::RCP<PGRAPH>             & NewPgraph)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldPgraph,
                    *NewCommDev,
                    *NewPgraph);
}
    
template<typename ITEM, typename COLMAP>
void
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
expandCommunicator(const bool         & isActive,
                   const communicator & OldCommDev,
                   const PGRAPH       & OldPgraph,
                   const communicator & NewCommDev,
                         PGRAPH       & NewPgraph)
{
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() > OldCommDev.size());
    
    //Alloc
    NewPgraph.clear();
    NewPgraph = OldPgraph;
    
    UInt pid = NewCommDev.rank();
    
    //Main loop
    for(UInt i=1; i <= NewPgraph.getRowMap().size(); ++i)
    { NewPgraph.getRowMap().get(i).setPid(pid); }
    
    for(UInt i=1; i <= NewPgraph.getColMap().size(); ++i)
    { NewPgraph.getColMap().get(i).setPid(pid); }
    
    NewPgraph.updateRowFinder();
    NewPgraph.updateColFinder();
  }
}


template<typename ITEM, typename COLMAP>
size_t
pGraphGlobalManip<ITEM,pMapItemShare,COLMAP>::
memSize() const
{
  return(sizeof(bool));
}

#endif
