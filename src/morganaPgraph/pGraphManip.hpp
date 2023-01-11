/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PGRAPHMANIP_HPP
#define PGRAPHMANIP_HPP

#include <set>
#include <assert.h>
#include <iostream>

#include "Teuchos_RCP.hpp"

#include "pGraph.hpp"
#include "pMapManip.hpp"
#include "pVectManip.hpp"

using namespace std;


/*! Serial manipulation for \c pGraph */
template<typename ITEM, typename ROWMAP, typename COLMAP> class pGraphManip
{
    /*! @name Typedefs */ //@{
  public:
    typedef pGraph<ITEM,ROWMAP,COLMAP>   PGRAPH;
    typedef pMap<pMapItemSendRecv>     SENDRECV;
    typedef pMap<ROWMAP>       CONTAINER_ROWMAP;
    typedef pMap<COLMAP>       CONTAINER_COLMAP;
    typedef pGraph<ITEM,COLMAP,ROWMAP> INVGRAPH;
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<PGRAPH> graph;
    //@}
    
    /*! @name Internal data - finder */ //@{
  public:
    bool graphLoaded;
    bool finderOk;
    map<ITEM,ROWMAP> container;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    pGraphManip();
    pGraphManip(const Teuchos::RCP<PGRAPH> & Graph);
    pGraphManip(PGRAPH & Graph);
    void setMap(const Teuchos::RCP<PGRAPH> & Graph);
    void setMap(PGRAPH & Graph);
    //@}
    
    /*! @name Search functions */ //@{
  public:
    /*! Build the graph finder */
    void buildFinder();
    
    /*! Reset the finder */
    void resetFinder();
    
    /*! Given a \c pGraphItem determines whether this is present in the \c pGraph */
    bool isItem(const ITEM & Item) const;
    
    /*! Given a \c pGraphItem returns the \c lid of this graphItem in the graph */
    UInt getLidItem(const ITEM & Item) const;
    
    /*! Given a \c pGraphItem returns the \c gid of this graphItem in the graph */
    UInt getGidItem(const ITEM & Item) const;
    
    /*! Given a \c pGraphItem returns the \c mapItem of this graphItem in the graph */
    ROWMAP getMapItem(const ITEM & Item) const;
    //@}
    
    
    /*! @name Set-Get and indexing functions */ //@{
  public:
    /*! Set \c map to normal indexing */
    void setNormalIndexing();
    void setNormalIndexing(PGRAPH & Graph) const;
    void setNormalIndexing(Teuchos::RCP<PGRAPH> & Graph) const;
    
    /*! Re-order the position of the items in the vector using the \c lid data */
    void setIndexing();
    void setIndexing(PGRAPH & Graph) const;
    void setIndexing(Teuchos::RCP<PGRAPH> & Graph) const;
    
    /*! Returns the maximum row \c gid of \c graph */
    UInt getMaxRowGid() const;
    UInt getMaxRowGid(const PGRAPH & Graph) const;
    UInt getMaxRowGid(const Teuchos::RCP<const PGRAPH> & Graph) const;
    
    /*! Returns the maximum col \c gid of \c graph */
    UInt getMaxColGid() const;
    UInt getMaxColGid(const PGRAPH & Graph) const;
    UInt getMaxColGid(const Teuchos::RCP<const PGRAPH> & Graph) const;
    
    /*! Check that all the local cids in the graph are mapped by the col map */
    bool checkCol() const;
    bool checkCol(const PGRAPH & Graph) const;
    bool checkCol(const Teuchos::RCP<const PGRAPH> & Graph) const;
    //@}
    
    
    /*! @name Manipulation functions */ //@{
  public:
    /*! Creates the inverse of the local graph  */
    void inversion(INVGRAPH & Graph) const;
    void inversion(const PGRAPH & inGraph, INVGRAPH & outGraph) const;
    void inversion(const Teuchos::RCP<const PGRAPH> & inGraph, Teuchos::RCP<INVGRAPH> & outGraph) const;
    
    /*! Given a list \c lidIndices of row \c lid s generates a \c pGraph containing the union of all the items present in \c lidIndices
    and the items directly connected to them */
    void getConnections(PGRAPH & Graph, const sVect<UInt> & lidIndices) const;
    void getConnections(const PGRAPH & inGraph, PGRAPH & outGraph, const sVect<UInt> & lidIndices) const;
    void getConnections(const Teuchos::RCP<const PGRAPH> & inGraph, Teuchos::RCP<INVGRAPH> & outGraph, const sVect<UInt> & lidIndices) const;
    
    /*! Merges to graphs, the row and coloums numbering of the internal \c graph are retained whyle the ones of \c Graph are modified to obtain
    a merged graph. The final result is available into \c graph */
    void mergeGraph(const PGRAPH & addGraph);
    void mergeGraph(PGRAPH & inGraph, const PGRAPH & addGraph) const;
    void mergeGraph(Teuchos::RCP<PGRAPH> & inGraph, const Teuchos::RCP<const PGRAPH> & addGraph) const;
    //@}
    
    /*! @name Mem functions */ //@{
  public:
    size_t memSize() const;
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
pGraphManip<ITEM,ROWMAP,COLMAP>::
pGraphManip()
{
  graphLoaded = false;
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
pGraphManip<ITEM,ROWMAP,COLMAP>::
pGraphManip(const Teuchos::RCP<PGRAPH> & Graph)
{
  graphLoaded = true;
  graph = Graph;
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
pGraphManip<ITEM,ROWMAP,COLMAP>::
pGraphManip(PGRAPH & Graph)
{
  graphLoaded = true;
  graph = Teuchos::rcpFromRef(Graph);
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
setMap(const Teuchos::RCP<PGRAPH> & Graph)
{
  graphLoaded = true;
  graph = Graph;
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
setMap(PGRAPH & Graph)
{
  graphLoaded = true;
  graph = Teuchos::rcpFromRef(Graph);
}



//_________________________________________________________________________________________________
// SEARCH FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
buildFinder()
{
  typedef pair<ITEM,ROWMAP> VALUETYPE;
  
  assert(graphLoaded);
  
  finderOk = true;
  VALUETYPE value;
  
  container.clear();
  
  for(UInt i=1; i<=graph->size(); ++i)
  {
    value.first  = graph->getDataL(i);
    value.second = graph->getMapL(i);
    
    container.insert(value);
  }
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
resetFinder()
{
  finderOk = false;
  container.clear();
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
bool
pGraphManip<ITEM,ROWMAP,COLMAP>::
isItem(const ITEM & Item) const
{
  assert(graphLoaded);
  assert(finderOk);
  
  return((bool)container.count(Item));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt
pGraphManip<ITEM,ROWMAP,COLMAP>::
getLidItem(const ITEM & Item) const
{
  assert(graphLoaded);
  assert(finderOk);
  assert((bool)container.count(Item));
  
  typedef typename map<ITEM,ROWMAP>::const_iterator ITERATOR;
  typedef          pair<ITEM,ROWMAP>                VALUETYPE;
  
  ITERATOR  iteratore = container.find(Item);
  VALUETYPE value     = *iteratore;
  ROWMAP    mItem     = value.second;
  
  return(mItem.getLid());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt
pGraphManip<ITEM,ROWMAP,COLMAP>::
getGidItem(const ITEM & Item) const
{
  assert(graphLoaded);
  assert(finderOk);
  assert((bool)container.count(Item));
  
  typedef typename map<ITEM,ROWMAP>::const_iterator ITERATOR;
  typedef          pair<ITEM,ROWMAP>                VALUETYPE;
  
  ITERATOR  iteratore = container.find(Item);
  VALUETYPE value     = *iteratore;
  ROWMAP    mItem     = value.second;
  
  return(mItem.getGid());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
ROWMAP
pGraphManip<ITEM,ROWMAP,COLMAP>::
getMapItem(const ITEM & Item) const
{
  assert(graphLoaded);
  assert(finderOk);
  assert((bool)container.count(Item));
  
  typedef typename map<ITEM,ROWMAP>::const_iterator ITERATOR;
  typedef          pair<ITEM,ROWMAP>                VALUETYPE;
  
  ITERATOR  iteratore = container.find(Item);
  VALUETYPE value     = *iteratore;
  ROWMAP    mItem     = value.second;
  
  return(mItem);
}



//_________________________________________________________________________________________________
// SET-GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
setNormalIndexing()
{
  assert(graphLoaded);
  
  //Row normal numbering
  pVectManip<ITEM,ROWMAP> manipulatorRow;
  manipulatorRow.setNormalIndexing(*graph);
  
  //Col normal numbering
  pMapManip<COLMAP> manipulatorCol;
  manipulatorCol.setNormalIndexing(graph->getColMap());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
setNormalIndexing(PGRAPH & Graph) const
{
  //Row normal numbering
  pVectManip<ITEM,ROWMAP> manipulatorRow;
  manipulatorRow.setNormalIndexing(Graph);
  
  //Col normal numbering
  pMapManip<COLMAP> manipulatorCol;
  manipulatorCol.setNormalIndexing(Graph.getColMap());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
setNormalIndexing(Teuchos::RCP<PGRAPH> & Graph) const
{
  setNormalIndexing(*Graph);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
setIndexing()
{
  assert(graphLoaded);
  
  //Row indexing
  pVectManip<ITEM,ROWMAP> manipulator;
  manipulator.setIndexing(*graph);
  
  //Col indexing
  pMapManip<COLMAP> manipulatorCol;
  manipulatorCol.setIndexing(graph->getColMap());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
setIndexing(PGRAPH & Graph) const
{
  //Row indexing
  pVectManip<ITEM,ROWMAP> manipulator;
  manipulator.setIndexing(Graph);
  
  //Col indexing
  pMapManip<COLMAP> manipulatorCol;
  manipulatorCol.setIndexing(Graph.getColMap());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
setIndexing(Teuchos::RCP<PGRAPH> & Graph) const
{
  setIndexing(*Graph);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt
pGraphManip<ITEM,ROWMAP,COLMAP>::
getMaxRowGid() const
{
  assert(graphLoaded);
  
  pVectManip<ITEM,ROWMAP> manipulator;
  return(manipulator.getMaxGid(*graph));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt
pGraphManip<ITEM,ROWMAP,COLMAP>::
getMaxRowGid(const PGRAPH & Graph) const
{
  pVectManip<ITEM,ROWMAP> manipulator;
  return(manipulator.getMaxGid(Graph));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt
pGraphManip<ITEM,ROWMAP,COLMAP>::
getMaxRowGid(const Teuchos::RCP<const PGRAPH> & Graph) const
{
  return(getMaxRowGid(*Graph));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt
pGraphManip<ITEM,ROWMAP,COLMAP>::
getMaxColGid() const
{
  assert(graphLoaded);
  
  pMapManip<COLMAP> manipulator;
  return(manipulator.getMaxGid(graph->getColMap()));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt
pGraphManip<ITEM,ROWMAP,COLMAP>::
getMaxColGid(const PGRAPH & Graph) const
{
  pMapManip<COLMAP> manipulator;
  return(manipulator.getMaxGid(Graph.getColMap()));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt
pGraphManip<ITEM,ROWMAP,COLMAP>::
getMaxColGid(const Teuchos::RCP<const PGRAPH> & Graph) const
{
  return(getMaxColGid(*Graph));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
bool
pGraphManip<ITEM,ROWMAP,COLMAP>::
checkCol() const
{
  assert(graphLoaded);
  
  UInt lid;
  bool flag = true;
  
  for(UInt i=1; i <= graph->size(); ++i)
  {
    for(UInt j=1; j <= graph->rowSizeL(i); ++j)
    {
      lid = graph->getCid_LL(i,j);
      flag = flag & (lid >= 1) & (lid <= graph->colSize());
    }
  }
  
  return(flag);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
bool
pGraphManip<ITEM,ROWMAP,COLMAP>::
checkCol(const PGRAPH & Graph) const
{
  UInt lid;
  bool flag = true;
  
  for(UInt i=1; i <= Graph.size(); ++i)
  {
    for(UInt j=1; j <= Graph.rowSizeL(i); ++j)
    {
      lid = Graph.getCid_LL(i,j);
      flag = flag & (lid >= 1) & (lid <= Graph.colSize());
    }
  }
  
  return(flag);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
bool
pGraphManip<ITEM,ROWMAP,COLMAP>::
checkCol(const Teuchos::RCP<const PGRAPH> & Graph) const
{
  return(checkCol(*Graph));
}



//_________________________________________________________________________________________________
// MANIPULATION FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
inversion(INVGRAPH & Graph) const
{
  assert(graphLoaded);
  
  //Graph dimension
  UInt maxDim = 0;
  UInt itemSize, k;
  
  for(UInt i=1; i<=graph->size(); ++i)
  {
    itemSize = graph->get(i).size();
    
    for(UInt j=1; j <= itemSize; ++j)
    { maxDim = max(maxDim,graph->getCid_LL(i,j)); }
  }
  
  Graph.clear();
  Graph.resize(maxDim);
  
  //Switching maps
  Graph.setRowMap(graph->getColMap());
  Graph.setColMap(graph->getRowMap()); 
  
  //Graph construction
  for(UInt i=1; i<=graph->size(); ++i)
  {
    itemSize = graph->get(i).size();
    
    for(UInt j=1; j <= itemSize; ++j)
    {
       k = graph->getCid_LL(i,j);
       Graph.get(k).push_back(i);
    }
  }
  
  Graph.updateFinder();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
inversion(const PGRAPH & inGraph, INVGRAPH & outGraph) const
{
  //Graph dimension
  UInt maxDim = 0;
  UInt itemSize, k;
  
  for(UInt i=1; i <= inGraph.size(); ++i)
  {
    itemSize = inGraph.get(i).size();
    
    for(UInt j=1; j <= itemSize; ++j)
    { maxDim = max(maxDim, inGraph.getCid_LL(i,j)); }
  }
  
  outGraph.clear();
  outGraph.resize(maxDim);
  
  //Switching maps
  outGraph.setRowMap(inGraph.getColMap());
  outGraph.setColMap(inGraph.getRowMap()); 
  
  //Graph construction
  for(UInt i=1; i <= inGraph.size(); ++i)
  {
    itemSize = inGraph.get(i).size();
    
    for(UInt j=1; j <= itemSize; ++j)
    {
       k = inGraph.getCid_LL(i,j);
       outGraph.get(k).push_back(i);
    }
  }
  
  outGraph.updateFinder();
}
 
template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
inversion(const Teuchos::RCP<const PGRAPH> & inGraph, Teuchos::RCP<INVGRAPH> & outGraph) const
{
  inversion(inGraph,*outGraph);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
getConnections(PGRAPH & Graph, const sVect<UInt> & lidIndices) const
{
  //Graph clearing
  Graph.clear();
  
  //Build codomain indices
  UInt lid, itemSize;
  set<UInt> codomIndices;
  
  for(UInt i=1; i <= lidIndices.size(); ++i)
  {
    lid      = lidIndices(i);
    itemSize = graph->get(lid).size();
    
    for(UInt j=1; j <= itemSize; ++j)
    { codomIndices.insert(graph->getCid_LL(lid,j)); }
  }
  
  //Build sub-Graph
  typedef typename set<ITEM>::iterator ITERATOR;
  typedef pair<ITERATOR,bool>          PAIR;
  
  UInt      cid;
  set<ITEM> itemList;
  ITERATOR  iter;
  PAIR      pair;
  CONTAINER_ROWMAP newRowMap;
  
  for(UInt i=1; i <= graph->size(); ++i)
  {
    itemSize = graph->get(i).size();
    
    for(UInt j=1; j <= itemSize; ++j)
    {
      cid = graph->getCid_LL(i,j);
      
      if(codomIndices.count(cid) >= 1)
      { 
	pair = itemList.insert(graph->get(i));
	
	if(pair.second)
	{ newRowMap.push_back(graph->getMapL(i)); }
      }
    }
  }
  
  assert(newRowMap.size() == itemList.size());
  
  //Final building
  UInt k=1;
  
  for(iter = itemList.begin(); iter != itemList.end(); iter++)
  {
    Graph.push_back(*iter,newRowMap(k),false);
    ++k;
  }
  
  Graph.updateFinder();
  Graph.setColMap(graph->getColMap());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
getConnections(const PGRAPH & inGraph, PGRAPH & outGraph, const sVect<UInt> & lidIndices) const
{
  //Graph clearing
  outGraph.clear();
  
  //Build codomain indices
  UInt lid, itemSize;
  set<UInt> codomIndices;
  
  for(UInt i=1; i <= lidIndices.size(); ++i)
  {
    lid      = lidIndices(i);
    itemSize = inGraph.get(lid).size();
    
    for(UInt j=1; j <= itemSize; ++j)
    { codomIndices.insert(inGraph.getCid_LL(lid,j)); }
  }
  
  //Build sub-Graph
  typedef typename set<ITEM>::iterator ITERATOR;
  typedef pair<ITERATOR,bool>          PAIR;
  
  UInt      cid;
  set<ITEM> itemList;
  ITERATOR  iter;
  PAIR      pair;
  CONTAINER_ROWMAP newRowMap;
  
  for(UInt i=1; i <= inGraph.size(); ++i)
  {
    itemSize = inGraph.get(i).size();
    
    for(UInt j=1; j <= itemSize; ++j)
    {
      cid = inGraph.getCid_LL(i,j);
      
      if(codomIndices.count(cid) >= 1)
      { 
	pair = itemList.insert(inGraph.get(i));
	
	if(pair.second)
	{ newRowMap.push_back(inGraph.getMapL(i)); }
      }
    }
  }
  
  assert(newRowMap.size() == itemList.size());
  
  //Final building
  UInt k=1;
  
  for(iter = itemList.begin(); iter != itemList.end(); iter++)
  {
    outGraph.push_back(*iter,newRowMap(k),false);
    ++k;
  }
  
  outGraph.updateFinder();
  outGraph.setColMap(inGraph.getColMap());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
getConnections(const Teuchos::RCP<const PGRAPH> & inGraph, Teuchos::RCP<INVGRAPH> & outGraph, const sVect<UInt> & lidIndices) const
{
  getConnections(*inGraph,*outGraph,lidIndices);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
mergeGraph(const PGRAPH & Graph)
{
  //Data ordering
  typedef typename set<ROWMAP>::iterator  ITERATOR_ROWMAP;
  typedef typename set<COLMAP>::iterator  ITERATOR_COLMAP;
  typedef pair<ITERATOR_ROWMAP,bool>      PAIR_ROWMAP;
  typedef pair<ITERATOR_COLMAP,bool>      PAIR_COLMAP;
  
  UInt maxRowLid = graph->rowSize();
  UInt maxColLid = graph->colSize();
  
  set<ROWMAP> rowContainer;
  set<COLMAP> colContainer;
  
  for(UInt i=1; i <= graph->rowSize(); ++i)
  { rowContainer.insert(graph->getRowMap().get(i)); }
  
  for(UInt i=1; i <= graph->colSize(); ++i)
  { colContainer.insert(graph->getColMap().get(i)); }
  
  //Graph insertion
  PAIR_ROWMAP pairRowMap;
  PAIR_COLMAP pairColMap;
  
  ITEM   dataItem;
  ROWMAP rowMapItem;
  COLMAP colMapItem;
  
  for(UInt i=1; i <= Graph.size(); ++i)
  {
    rowMapItem = Graph.getMapL(i);
    pairRowMap = rowContainer.insert(rowMapItem);
    
    //The rowItem is not present
    if(pairRowMap.second)
    {
      //Change rowMap settings
      maxRowLid++;
      rowMapItem.setLid(maxRowLid);
      
      //ColItems cycle
      dataItem = Graph.getDataL(i);
      
      for(UInt j=1; j <= dataItem.size(); ++j)
      {
	assert(dataItem(j) >= 1);
	assert(dataItem(j) <= Graph.colSize());
	
	colMapItem = Graph.getColMap().get(dataItem(j));
	pairColMap = colContainer.insert(colMapItem);
	
	if(pairColMap.second) //ColItem not already present
	{
	  //Cancel the inserted item
	  colContainer.erase(colMapItem);
	  
	  //Update lid
	  maxColLid++;
	  colMapItem.setLid(maxColLid);
	  
	  //Update in the graphItem
	  dataItem(j) = maxColLid;
	  
	  //Insert in the colMap 
	  graph->getColMap().push_back(colMapItem);
	  
	  //Update the item in the set
	  colContainer.insert(colMapItem);
	}
	else //ColItem already present
	{
	  //Identify the item already present
	  colMapItem = *pairColMap.first;
	  
	  //Update in the graphItem
	  dataItem(j) = colMapItem.getLid();
	}
      }
      
      //Insert the new rowMapItem and dataItem
      graph->push_back(dataItem,rowMapItem);
    }
  }
  
  //Update the finder of the new graph 
  graph->updateFinder();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
mergeGraph(PGRAPH & inGraph, const PGRAPH & addGraph) const
{
  //Data ordering
  typedef typename set<ROWMAP>::iterator  ITERATOR_ROWMAP;
  typedef typename set<COLMAP>::iterator  ITERATOR_COLMAP;
  typedef pair<ITERATOR_ROWMAP,bool>      PAIR_ROWMAP;
  typedef pair<ITERATOR_COLMAP,bool>      PAIR_COLMAP;
  
  UInt maxRowLid = inGraph.rowSize();
  UInt maxColLid = inGraph.colSize();
  
  set<ROWMAP> rowContainer;
  set<COLMAP> colContainer;
  
  for(UInt i=1; i <= inGraph.rowSize(); ++i)
  { rowContainer.insert(inGraph.getRowMap().get(i)); }
  
  for(UInt i=1; i <= inGraph.colSize(); ++i)
  { colContainer.insert(inGraph.getColMap().get(i)); }
  
  //Graph insertion
  PAIR_ROWMAP pairRowMap;
  PAIR_COLMAP pairColMap;
  
  ITEM   dataItem;
  ROWMAP rowMapItem;
  COLMAP colMapItem;
  
  for(UInt i=1; i <= addGraph.size(); ++i)
  {
    rowMapItem = addGraph.getMapL(i);
    pairRowMap = rowContainer.insert(rowMapItem);
    
    //The rowItem is not present
    if(pairRowMap.second)
    {
      //Change rowMap settings
      maxRowLid++;
      rowMapItem.setLid(maxRowLid);
      
      //ColItems cycle
      dataItem = addGraph.getDataL(i);
      
      for(UInt j=1; j <= dataItem.size(); ++j)
      {
        assert(dataItem(j) >= 1);
        assert(dataItem(j) <= addGraph.colSize());

        colMapItem = addGraph.getColMap().get(dataItem(j));
        pairColMap = colContainer.insert(colMapItem);

        if(pairColMap.second) //ColItem not already present
        {
          //Cancel the inserted item
          colContainer.erase(colMapItem);

          //Update lid
          maxColLid++;
          colMapItem.setLid(maxColLid);

          //Update in the graphItem
          dataItem(j) = maxColLid;

          //Insert in the colMap 
          inGraph.getColMap().push_back(colMapItem);

          //Update the item in the set
          colContainer.insert(colMapItem);
        }
        else //ColItem already present
        {
          //Identify the item already present
          colMapItem = *pairColMap.first;

          //Update in the graphItem
          dataItem(j) = colMapItem.getLid();
        }
      }
      
      //Insert the new rowMapItem and dataItem
      inGraph.push_back(dataItem,rowMapItem);
    }
  }
  
  //Update the finder of the new graph 
  inGraph.updateFinder();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraphManip<ITEM,ROWMAP,COLMAP>::
mergeGraph(Teuchos::RCP<PGRAPH> & inGraph, const Teuchos::RCP<const PGRAPH> & addGraph) const
{
  mergeGraph(*inGraph,*addGraph);
}


//_________________________________________________________________________________________________
// MEM FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
size_t
pGraphManip<ITEM,ROWMAP,COLMAP>::
memSize() const
{
  if(container.size() == 0)
  { return(2 * sizeof(bool)); }
  else
  { 
    return( 2 * sizeof(bool)
            + (container.begin()->memSize()
             + 2 * sizeof(Real*) ) * container.size());
  }
}

#endif
