/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef _PGRAPH_HPP
#define _PGRAPH_HPP

#include <vector>
#include <assert.h>
#include <iostream>

#include "pMap.hpp"
#include "pVect.hpp"

#include "morganaTypes.hpp"
#include "simpleFormats.hpp"

using namespace std;


/*! Container for parallel Graph implementation */
template<typename ITEM, typename ROWMAP, typename COLMAP> class pGraph : public pVect<ITEM,ROWMAP>
{
  public:
    typedef pVect<ITEM,ROWMAP> PVECT;
    typedef pMap<ROWMAP>       CONTAINER_ROWMAP;
    typedef pMap<COLMAP>       CONTAINER_COLMAP;
    typedef sVect<ITEM>        CONTAINER_DATA;
    typedef pMapManip<COLMAP>  MANIPULATOR;
  
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    /*! @name Internal data */ //@{
  public:
    bool             isLocal;       //! Determines whether the graph-items have a local or global numbering
    bool             colStartupOk;  //! Determines whether the column finder has been bootstrapped
    CONTAINER_COLMAP colMap;        //! The coloumn map
    MANIPULATOR      colManip;      //! The manipulator for the column map
    //@}
    
    /*! @name Constructors and operators */ //@{
  public:
    /*! Constructor */
    pGraph();
    
    /*! Constructor
    \param n length of the graph */
    pGraph(const UInt & n);
    
    /*! Copy constructor */
    pGraph(const pGraph & G);
    
    /*! Constructor - \c startupOk = true, apply normal indexing, finder build */
    pGraph(const CONTAINER_ROWMAP & RowMap,
           const CONTAINER_COLMAP & ColMap,
           const CONTAINER_DATA   & Data,
           const bool             & IsLocal = true);
    
    /*! Constructor - \c startupOk = true, apply normal indexing, finder build */
    pGraph(const Teuchos::RCP<CONTAINER_ROWMAP> & RowMap,
           const Teuchos::RCP<CONTAINER_COLMAP> & ColMap,
           const Teuchos::RCP<CONTAINER_DATA>   & Data,
           const bool                           & IsLocal = true);
    
    /*! Equality operator */
    pGraph & operator=(const pGraph & G);
    //@}
    
    /*! @name Global-Local col-numbering */ //@{
  public:
    /*! For each row substitutes the local \c cids with the global ones. The local \c cid s are searched in the \c colMap and the corresponding global \c cid is written back 
    in the \c graphItem . If any local \c cid is not found in the \c colMap the algorithm fails */
    void pushToGlobal();
    
    /*! Is the opposite algorithm of \c pushToGlobal. The \c cid s in the \c graphItem s are interpreted as global and the corrensponding local \c cid s are searched
    in the \c colMap. If any global \c cid is not found in the \c colMap the algorithm fails */
    void pushToLocal();
    
    /*! If \c true the \c cid s in the \c graphItem s are intended as local. */
    const bool & colIsLocal() const;
    
    /*! If \c true the \c cid s in the \c graphItem s are intended as local. */
    bool & colIsLocal();
    
    /*! Resets the row finder */
    void resetRowFinder();
    
    /*! Updates the row finder */
    void updateRowFinder();
    
    /*! Reset the col finder */
    void resetColFinder();
    
    /*! Update the col finder */
    void updateColFinder();
    
    /*! Update the sorting of the items */
    void updateSorting();
    
    /*! Clear all the data */
    void clear();
    //@}
    
    /*! @name Set and Add functions */ //@{
  public:
    /*! Set data - \c startupOk = true, apply normal indexing, finder build */
    void setData(const CONTAINER_ROWMAP & RowMap,
                 const CONTAINER_COLMAP & ColMap,
                 const CONTAINER_DATA   & Data,
                 const bool             & IsLocal = true);
    
    /*! Add data block - \c startupOk = true, apply normal indexing, finder build */
    void setData(const Teuchos::RCP<CONTAINER_ROWMAP> & RowMap,
                 const Teuchos::RCP<CONTAINER_COLMAP> & ColMap,
                 const Teuchos::RCP<CONTAINER_DATA>   & Data,
                 const bool                           & IsLocal = true);
    
    /*! Add data block - \c startupOk = true, apply normal indexing, finder re-building */
    void addData(const CONTAINER_ROWMAP & RowMap,
                 const CONTAINER_COLMAP & ColMap,
                 const CONTAINER_DATA   & Data);
    
    /*! Add data block - \c startupOk = true, apply normal indexing, finder re-building */
    void addData(const Teuchos::RCP<CONTAINER_ROWMAP> & RowMap,
                 const Teuchos::RCP<CONTAINER_COLMAP> & ColMap,
                 const Teuchos::RCP<CONTAINER_DATA>   & Data);
    //@}
    
    /*! @name Interface functions */ //@{
  public:
    UInt & getLid(const UInt & i);
    UInt & getGid(const UInt & i);
    UInt & getCid_LL(const UInt & i, const UInt & j);
    UInt & getCid_LG(const UInt & i, const UInt & j);
    UInt & getCid_GL(const UInt & i, const UInt & j);
    UInt & getCid_GG(const UInt & i, const UInt & j);
    const UInt & getLid(const UInt & i) const;
    const UInt & getGid(const UInt & i) const;
    const UInt & getCid_LL(const UInt & i, const UInt & j) const;
    const UInt & getCid_LG(const UInt & i, const UInt & j) const;
    const UInt & getCid_GL(const UInt & i, const UInt & j) const;
    const UInt & getCid_GG(const UInt & i, const UInt & j) const;
    UInt rowSizeL(const UInt & i) const;
    UInt rowSizeG(const UInt & i) const;
    UInt rowSize() const;
    UInt colSize() const;
    //@}
    
    /*! @name Map interface functions */ //@{
  public:
    ITEM & getItemL(const UInt & lid);
    const ITEM & getItemL(const UInt & lid) const;
    ITEM & getItemG(const UInt & gid);
    const ITEM & getItemG(const UInt & gid) const;
    bool isRowG(const UInt & gid) const;
    ROWMAP & getRowMapL(const UInt & lid);
    const ROWMAP & getRowMapL(const UInt & lid) const;
    ROWMAP & getRowMapG(const UInt & gid);
    const ROWMAP & getRowMapG(const UInt & gid) const;
    bool isColG(const UInt & gid) const;
    COLMAP & getColMapL(const UInt & lid);
    const COLMAP & getColMapL(const UInt & lid) const;
    COLMAP & getColMapG(const UInt & gid);
    const COLMAP & getColMapG(const UInt & gid) const;
    //@}
    
    /*! @name Block map functions */ //@{
  public:
    pMap<ROWMAP> & getRowMap();
    pMap<COLMAP> & getColMap();
    const pMap<ROWMAP> & getRowMap() const;
    const pMap<COLMAP> & getColMap() const;
    void setRowMap(const pMap<ROWMAP> & RowMap);
    void setColMap(const pMap<COLMAP> & ColMap);
    //@}
    
    /*! @name Printout */ //@{
  public:
    template<typename G, typename RM, typename CM>
    friend ostream & operator<<(ostream & f, const pGraph<G,RM,CM> & M);
    //@}
};


//_______________________________________________________________________________________________________
// CONSTRUCTORS AND OPERATORS
//-------------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
pGraph<ITEM,ROWMAP,COLMAP>::
pGraph() : pVect<ITEM,ROWMAP>()
{
  isLocal      = true;
  colStartupOk = false;
  
  colManip.setMap(colMap);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
pGraph<ITEM,ROWMAP,COLMAP>::
pGraph(const UInt & n) : pVect<ITEM,ROWMAP>(n)
{
  isLocal      = true;
  colStartupOk = false;
  
  colManip.setMap(colMap);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
pGraph<ITEM,ROWMAP,COLMAP>::
pGraph(const pGraph & G) : pVect<ITEM,ROWMAP>(G)
{
  isLocal = G.isLocal;
  colMap  = G.colMap;
  
  colStartupOk = true;
  colMap       = G.colMap;
  
  colManip.setMap(colMap);
  colManip.setNormalIndexing();
  colManip.buildFinder();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
pGraph<ITEM,ROWMAP,COLMAP>::
pGraph(const CONTAINER_ROWMAP & RowMap,
       const CONTAINER_COLMAP & ColMap,
       const CONTAINER_DATA   & Data,
       const bool             & IsLocal) : pVect<ITEM,ROWMAP>(RowMap,Data)
{
  isLocal = IsLocal;
  
  colMap       = ColMap;
  colStartupOk = true;
  
  colManip.setMap(colMap);
  colManip.setNormalIndexing();
  colManip.buildFinder();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
pGraph<ITEM,ROWMAP,COLMAP>::
pGraph(const Teuchos::RCP<CONTAINER_ROWMAP> & RowMap,
       const Teuchos::RCP<CONTAINER_COLMAP> & ColMap,
       const Teuchos::RCP<CONTAINER_DATA>   & Data,
       const bool                           & IsLocal) : pVect<ITEM,ROWMAP>(RowMap,Data)
{
  isLocal = IsLocal;
  
  colMap       = *ColMap;
  colStartupOk = true;
  
  colManip.setMap(colMap);
  colManip.setNormalIndexing();
  colManip.buildFinder();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
pGraph<ITEM,ROWMAP,COLMAP> &
pGraph<ITEM,ROWMAP,COLMAP>::
operator=(const pGraph & G)
{ 
  PVECT::operator=(G);
  isLocal = G.isLocal;
  
  colMap       = G.colMap;
  colStartupOk = true;
  
  colManip.resetFinder();
  colManip.setNormalIndexing();
  colManip.buildFinder();
  
  return(*this);
}



//_______________________________________________________________________________________________________
// GLOBAL-LOCAL NUMBERING
//-------------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraph<ITEM,ROWMAP,COLMAP>::
pushToGlobal()
{
  assert(isLocal);
  isLocal = false;
  
  UInt lid, gid;
  
  for(UInt i=1; i <= PVECT::data->size(); ++i)
  {
    for(UInt j=1; j <= PVECT::data->get(i).size(); ++j)
    {
      lid = PVECT::data->get(i).getCid(j);
      
      assert(lid >= 1);
      assert(lid <= colMap.size());
      
      gid = colMap(lid).getGid();
      PVECT::data->get(i).getCid(j) = gid;
    }
    
    PVECT::data->get(i).updateSorting();
  }
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraph<ITEM,ROWMAP,COLMAP>::
pushToLocal()
{
  assert(!isLocal);
  isLocal = true;
  
  //Allocations
  COLMAP item;
  UInt lid, gid;
  
  //ColManipulator
  MANIPULATOR manipulator(colMap);
  manipulator.buildFinder();
  
  //Updating of the gids
  for(UInt i=1; i <= PVECT::data->size(); ++i)
  {
    for(UInt j=1; j <= PVECT::data->get(i).size(); ++j)
    {
      gid = PVECT::data->get(i).getCid(j);
      item.setGid(gid);
      
      assert(manipulator.isItem(item));
      
      lid = manipulator.getLidItem(item);
      PVECT::data->get(i).getCid(j) = lid;
    }
    
    PVECT::data->get(i).updateSorting();
  }
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
const bool &
pGraph<ITEM,ROWMAP,COLMAP>::
colIsLocal() const
{
  return(isLocal);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
bool &
pGraph<ITEM,ROWMAP,COLMAP>::
colIsLocal()
{
  return(isLocal);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraph<ITEM,ROWMAP,COLMAP>::
resetRowFinder()
{
  PVECT::resetFinder();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraph<ITEM,ROWMAP,COLMAP>::
updateRowFinder()
{
  PVECT::updateFinder();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraph<ITEM,ROWMAP,COLMAP>::
resetColFinder()
{
  colStartupOk = false;
  colManip.resetFinder();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraph<ITEM,ROWMAP,COLMAP>::
updateColFinder()
{
  colStartupOk = true;
  
  colManip.resetFinder();
  colManip.setNormalIndexing();
  colManip.buildFinder();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraph<ITEM,ROWMAP,COLMAP>::
updateSorting()
{
  for(UInt i=1; i <= PVECT::data->size(); ++i)
  { PVECT::data->get(i).updateSorting(); }
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraph<ITEM,ROWMAP,COLMAP>::
clear()
{
  colStartupOk = false;
  
  PVECT::clear();
  colManip.resetFinder();
  colMap.clear();
}


//_______________________________________________________________________________________________________
// SET AND ADD FUNCTIONS
//-------------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraph<ITEM,ROWMAP,COLMAP>::
setData(const CONTAINER_ROWMAP & RowMap,
        const CONTAINER_COLMAP & ColMap,
        const CONTAINER_DATA   & Data,
        const bool             & IsLocal)
{
  PVECT::setData(RowMap,Data);
  isLocal = IsLocal;
  
  colStartupOk = true;
  colMap       = ColMap;
  
  colManip.resetFinder();
  colManip.setNormalIndexing();
  colManip.buildFinder();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraph<ITEM,ROWMAP,COLMAP>::
setData(const Teuchos::RCP<CONTAINER_ROWMAP> & RowMap,
        const Teuchos::RCP<CONTAINER_COLMAP> & ColMap,
        const Teuchos::RCP<CONTAINER_DATA>   & Data,
        const bool                           & IsLocal)
{
  PVECT::setData(RowMap,Data);
  isLocal = IsLocal;
  
  colStartupOk = true;
  colMap       = *ColMap;
  
  colManip.resetFinder();
  colManip.setNormalIndexing();
  colManip.buildFinder();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraph<ITEM,ROWMAP,COLMAP>::
addData(const CONTAINER_ROWMAP & RowMap,
        const CONTAINER_COLMAP & ColMap,
        const CONTAINER_DATA   & Data)
{
  PVECT::addData(RowMap,Data);
  
  for(UInt i=1; i <= ColMap.size(); ++i)
  { colMap.push_back(ColMap(i)); }
  
  colStartupOk = true;
  
  colManip.resetFinder();
  colManip.setNormalIndexing();
  colManip.buildFinder();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraph<ITEM,ROWMAP,COLMAP>::
addData(const Teuchos::RCP<CONTAINER_ROWMAP> & RowMap,
        const Teuchos::RCP<CONTAINER_COLMAP> & ColMap,
        const Teuchos::RCP<CONTAINER_DATA>   & Data)
{
  PVECT::addData(RowMap,Data);
  
  for(UInt i=1; i <= ColMap.size(); ++i)
  { colMap.push_back(ColMap(i)); }
  
  colStartupOk = true;
  
  colManip.resetFinder();
  colManip.setNormalIndexing();
  colManip.buildFinder();
}



//_______________________________________________________________________________________________________
// INTERFACE FUNCTIONS
//-------------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt &
pGraph<ITEM,ROWMAP,COLMAP>::
getLid(const UInt & i)
{
  assert(i <= PVECT::size());
  return(PVECT::map->get(i).getLid());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt &
pGraph<ITEM,ROWMAP,COLMAP>::
getGid(const UInt & i)
{
  assert(i <= PVECT::size());
  return(PVECT::map->get(i).getGid());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt &
pGraph<ITEM,ROWMAP,COLMAP>::
getCid_LL(const UInt & i, const UInt & j)
{
  assert(isLocal);
  
  assert(i >= 1);
  assert(i <= PVECT::size());
  assert(j >= 1);
  assert(j <= PVECT::getL(i).size());
  
  return(PVECT::getL(i).getCid(j));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt &
pGraph<ITEM,ROWMAP,COLMAP>::
getCid_LG(const UInt & i, const UInt & j)
{
  assert(isLocal);
  
  assert(i >= 1);
  assert(i <= PVECT::size());
  assert(j >= 1);
  assert(j <= PVECT::getL(i).size());  
  assert(PVECT::getL(i).getCid(j) >= 1);
  assert(PVECT::getL(i).getCid(j) <= colMap.size());
  
  return( colMap(PVECT::getL(i).getCid(j)).getGid() );
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt &
pGraph<ITEM,ROWMAP,COLMAP>::
getCid_GL(const UInt & i, const UInt & j)
{
  assert(isLocal);
  
  assert(i >= 1);
  assert(PVECT::isG(i));
  assert(j >= 1);
  assert(j <= PVECT::getG(i).size());
  
  return( PVECT::getG(i).getCid(j) );
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt &
pGraph<ITEM,ROWMAP,COLMAP>::
getCid_GG(const UInt & i, const UInt & j)
{
  assert(isLocal);
  
  assert(i >= 1);
  assert(PVECT::isG(i));
  assert(j >= 1);
  assert(j <= PVECT::getG(i).size());  
  assert(PVECT::getG(i).getCid(j) >= 1);
  assert(PVECT::getG(i).getCid(j) <= colMap.size());
  
  return( colMap(PVECT::getG(i).getCid(j)).getGid() );
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
const UInt &
pGraph<ITEM,ROWMAP,COLMAP>::
getLid(const UInt & i) const
{ 
  assert(i <= PVECT::size());
  
  return( PVECT::map->get(i).getLid() );
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
const UInt &
pGraph<ITEM,ROWMAP,COLMAP>::
getGid(const UInt & i) const
{
  assert(i <= PVECT::size());
  
  return(PVECT::map->get(i).getGid());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
const UInt &
pGraph<ITEM,ROWMAP,COLMAP>::
getCid_LL(const UInt & i, const UInt & j) const
{
  assert(isLocal);
  
  assert(i >= 1);
  assert(i <= PVECT::size());
  assert(j >= 1);
  assert(j <= PVECT::getL(i).size());
  
  return(PVECT::getL(i).getCid(j));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
const UInt &
pGraph<ITEM,ROWMAP,COLMAP>::
getCid_LG(const UInt & i, const UInt & j) const
{
  assert(isLocal);
  
  assert(i >= 1);
  assert(i <= PVECT::size());
  assert(j >= 1);
  assert(j <= PVECT::getL(i).size());  
  assert(PVECT::getL(i).getCid(j) >= 1);
  assert(PVECT::getL(i).getCid(j) <= colMap.size());
  
  return(colMap(PVECT::getL(i).getCid(j)).getGid());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
const UInt &
pGraph<ITEM,ROWMAP,COLMAP>::
getCid_GL(const UInt & i, const UInt & j) const
{
  assert(isLocal);
  
  assert(i >= 1);
  assert(PVECT::isG(i));
  assert(j >= 1);
  assert(j <= PVECT::getG(i).size());
  
  return(PVECT::getG(i).getCid(j));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
const UInt &
pGraph<ITEM,ROWMAP,COLMAP>::
getCid_GG(const UInt & i, const UInt & j) const
{
  assert(isLocal);
  
  assert(i >= 1);
  assert(i <= PVECT::size());
  assert(j >= 1);
  assert(j <= PVECT::getG(i).size());  
  assert(PVECT::getG(i).getCid(j) >= 1);
  assert(PVECT::getG(i).getCid(j) <= colMap.size());
  
  return(colMap(PVECT::getG(i).getCid(j)).getGid());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt
pGraph<ITEM,ROWMAP,COLMAP>::
rowSizeL(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= PVECT::size());
  
  return(PVECT::getL(i).size());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt
pGraph<ITEM,ROWMAP,COLMAP>::
rowSizeG(const UInt & i) const
{
  assert(PVECT::isG(i));
  
  return(PVECT::getG(i).size());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt
pGraph<ITEM,ROWMAP,COLMAP>::
rowSize() const
{
  return(PVECT::size());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
UInt
pGraph<ITEM,ROWMAP,COLMAP>::
colSize() const
{
  return(colMap.size());
}



//_______________________________________________________________________________________________________
// BLOCK MAP FUNCTIONS
//-------------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
ITEM &
pGraph<ITEM,ROWMAP,COLMAP>::
getItemL(const UInt & lid)
{
  assert(lid >= 1);
  assert(lid <= PVECT::size());
  
  return(PVECT::getL(lid));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
const ITEM &
pGraph<ITEM,ROWMAP,COLMAP>::
getItemL(const UInt & lid) const
{
  assert(lid >= 1);
  assert(lid <= PVECT::size());
  
  return(PVECT::getL(lid));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
ITEM &
pGraph<ITEM,ROWMAP,COLMAP>::
getItemG(const UInt & gid)
{
  assert(PVECT::isG(gid));
  return(PVECT::getG(gid));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
const ITEM &
pGraph<ITEM,ROWMAP,COLMAP>::
getItemG(const UInt & gid) const
{
  assert(PVECT::isG(gid));
  return(PVECT::getG(gid));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
bool
pGraph<ITEM,ROWMAP,COLMAP>::
isRowG(const UInt & gid) const
{
  assert(PVECT::isG(gid));
  return(PVECT::isG(gid));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
ROWMAP &
pGraph<ITEM,ROWMAP,COLMAP>::
getRowMapL(const UInt & lid)
{
  assert(lid >= 1);
  assert(lid <= PVECT::size());
  return(PVECT::getMapL(lid));
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
const ROWMAP &
pGraph<ITEM,ROWMAP,COLMAP>::
getRowMapL(const UInt & lid) const
{
  assert(lid >= 1);
  assert(lid <= PVECT::size());
  return(PVECT::getMapL(lid));
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
ROWMAP &
pGraph<ITEM,ROWMAP,COLMAP>::
getRowMapG(const UInt & gid)
{
  assert(PVECT::isG(gid));
  return(PVECT::getMapG(gid));
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
const ROWMAP &
pGraph<ITEM,ROWMAP,COLMAP>::
getRowMapG(const UInt & gid) const
{
  assert(PVECT::isG(gid));
  return(PVECT::getMapG(gid));
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
bool
pGraph<ITEM,ROWMAP,COLMAP>::
isColG(const UInt & gid) const
{
  assert(colStartupOk);
  return( colManip->isItem(COLMAP(0,gid)) );
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
COLMAP &
pGraph<ITEM,ROWMAP,COLMAP>::
getColMapL(const UInt & lid)
{
  assert(colStartupOk);
  assert(lid >= 1);
  assert(lid <= colMap.size());
  
  return(colMap.get(lid));
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
const COLMAP &
pGraph<ITEM,ROWMAP,COLMAP>::
getColMapL(const UInt & lid) const
{
  assert(colStartupOk);
  assert(lid >= 1);
  assert(lid <= colMap.size());
  
  return(colMap.get(lid));
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
COLMAP &
pGraph<ITEM,ROWMAP,COLMAP>::
getColMapG(const UInt & gid)
{
  assert(colStartupOk);
  assert(isColG(gid));
  
  UInt lid = colMap.getLidItem(COLMAP(0,gid));
  return(colMap.get(lid));
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
const COLMAP &
pGraph<ITEM,ROWMAP,COLMAP>::
getColMapG(const UInt & gid) const
{
  assert(colStartupOk);
  assert(isColG(gid));
  
  UInt lid = colMap.getLidItem(COLMAP(0,gid));
  return(colMap.get(lid));
}



//_______________________________________________________________________________________________________
// BLOCK MAP FUNCTIONS
//-------------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
pMap<ROWMAP> &
pGraph<ITEM,ROWMAP,COLMAP>::
getRowMap()
{
  return(PVECT::getMapRef());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
pMap<COLMAP> &
pGraph<ITEM,ROWMAP,COLMAP>::
getColMap()
{
  return(colMap);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
const pMap<ROWMAP> &
pGraph<ITEM,ROWMAP,COLMAP>::
getRowMap() const
{
  return(PVECT::getMapRef());
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
const pMap<COLMAP> &
pGraph<ITEM,ROWMAP,COLMAP>::
getColMap() const
{
  return(colMap);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraph<ITEM,ROWMAP,COLMAP>::
setRowMap(const pMap<ROWMAP> & RowMap)
{
  *(PVECT::map) = RowMap;
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
pGraph<ITEM,ROWMAP,COLMAP>::
setColMap(const pMap<COLMAP> & ColMap)
{
  colMap = ColMap;
}



//_______________________________________________________________________________________________________
// PRINTOUT AND SERIALIZATION
//-------------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
template<class ARK>
void
pGraph<ITEM,ROWMAP,COLMAP>::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  PVECT::serialize(ar,version);
  ar & colMap;
}

template<typename G, typename RM, typename CM>
ostream & operator<<(ostream & f, const pGraph<G,RM,CM> & M)
{
  cout << "ROW" << endl;
  for(UInt i=1; i <= M.data->size(); ++i)
  {
    f << i << endl;
    f << " map:  " << M.map->get(i);
    f << " data: " << endl << M.data->get(i) << endl << endl;
  }
  
  cout << "COL" << endl;
  f << M.colMap << endl;
  
  return(f);
}

#endif
