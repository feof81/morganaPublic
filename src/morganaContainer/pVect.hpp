/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PVECT_HPP
#define PVECT_HPP

#include "simpleFormats.hpp"
#include "pMap.hpp"
#include "pMapManip.hpp"

using namespace std;


/*! Parallel vector data type, global and local id access support. Contains a finder device to find the global ids. */
template<typename DATA, typename MAP> class pVect
{  
    /*! @name Parallel support */ //@{
  public:
    typedef sVect<DATA>                  CONTAINER_DATA;
    typedef pMap<MAP>                    CONTAINER_MAP;
    typedef Teuchos::RCP<CONTAINER_DATA> RCP_CONTAINER_DATA;
    typedef Teuchos::RCP<CONTAINER_MAP>  RCP_CONTAINER_MAP;
    typedef pMapManip<MAP>               MANIPULATOR;
    //@}
    
    
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    
    /*! @name Internal data */ //@{
  public:
    Teuchos::RCP<CONTAINER_DATA>  data;
    Teuchos::RCP<CONTAINER_MAP>   map;
    Teuchos::RCP<MANIPULATOR>     mapManip;
    bool startupOk;
    //@}
    
    
    /*! @name Constructors and operators */ //@{
  public:
    /*! Constructor - \c startupOk = false */
    pVect();
    
    /*! Constructor - \c startupOk = false */
    pVect(const UInt & n);
    
    /*! Constructor - \c startupOk = true, apply normal indexing, finder build */
    pVect(const pVect & Vect);
    
    /*! Constructor - \c startupOk = true, apply normal indexing, finder build */
    pVect(const CONTAINER_MAP & Map, const CONTAINER_DATA & Data);
    
    /*! Constructor - \c startupOk = true, apply normal indexing, finder build */
    pVect(const Teuchos::RCP<CONTAINER_MAP> & Map, const Teuchos::RCP<CONTAINER_DATA> & Data);
    
    /*! Equality operator - \c startupOk = true, apply normal indexing, finder build */
    pVect operator=(const pVect & Vect);
    //@}
    
    
    /*! @name Set and Add functions */ //@{
  public:
    /*! Set data - \c startupOk = true, apply normal indexing, finder build */
    void setData(const CONTAINER_MAP & Map, const CONTAINER_DATA & Data);
    
    /*! Add data block - \c startupOk = true, apply normal indexing, finder build */
    void setData(const Teuchos::RCP<CONTAINER_MAP> & Map, const Teuchos::RCP<CONTAINER_DATA> & Data);
    
    /*! Add data block - \c startupOk = true, apply normal indexing, finder re-building */
    void addData(const CONTAINER_MAP & Map, const CONTAINER_DATA & Data);
    
    /*! Add data block - \c startupOk = true, apply normal indexing, finder re-building */
    void addData(const Teuchos::RCP<CONTAINER_MAP> & Map, const Teuchos::RCP<CONTAINER_DATA> & Data);
    
    /*! Cancel all the finder informations, \c startupOk = false */
    void resetFinder();
    
    /*! Re-build the finder, enforce normal indexing, \c startupOk = true */
    void updateFinder();
    
    /*! Buffering of all the lids */
    void bufferLids();
    
    /*! Restore the buffered lids */
    void restoreLids();
    
    /*!Clearing of all the informations - \c startupOk = false */
    void clear();
    //@}
    
    /*! @name Get - Set block functions */ //@{
  public:
    CONTAINER_DATA           & getDataRef();
    CONTAINER_MAP            & getMapRef();
    RCP_CONTAINER_DATA       & getDataRcp();
    RCP_CONTAINER_MAP        & getMapRcp();
    const RCP_CONTAINER_DATA & getDataRcp() const;
    const RCP_CONTAINER_MAP  & getMapRcp() const;
    const CONTAINER_DATA     & getDataRef() const;
    const CONTAINER_MAP      & getMapRef() const;
    void setMap(const CONTAINER_MAP & RowMap);
    void setData(const CONTAINER_DATA & Data);
    //@}
    
    /*! @name Get pedantic functions */ //@{
  public:
    UInt sizeL();
    DATA & getL(const UInt & lid);
    const DATA & getL(const UInt & lid) const;
    DATA & getG(const UInt & gid);
    const DATA & getG(const UInt & gid) const;
    bool isG(const UInt & gid) const;
    MAP & getMapL(const UInt & lid);
    const MAP & getMapL(const UInt & lid) const;
    MAP & getMapG(const UInt & gid);
    const MAP & getMapG(const UInt & gid) const;
    DATA & getDataL(const UInt & lid);
    const DATA & getDataL(const UInt & lid) const;
    DATA & getDataG(const UInt & gid);
    const DATA & getDataG(const UInt & gid) const;
    //@}
    
    /*! @name Get interface functions */ //@{
  public:
    UInt size() const;
    DATA & get(const UInt & lid);
    const DATA & get(const UInt & lid) const;
    DATA & operator() (const UInt & lid);
    const DATA & operator() (const UInt & lid) const;
    //@}
    
    /*! @name PushBack functions */ //@{
  public:
    /*! Add one item
    \param mapItem map for the item
    \param dataItem the datum to be added
    \param updateFinder if true the gid finder is updated,
    on the contrary is necessary to use the \c updateFinder method before accessing to the data (\c startupOk = false) */
    void push_back(const MAP & mapItem, const DATA & dataItem, bool updateFinder = false);
    
    /*! Add one item
    \param mapItem map for the item
    \param dataItem the datum to be added
    \param updateFinder if true the gid finder is updated,
    on the contrary is necessary to use the \c updateFinder method before accessing to the data (\c startupOk = false) */
    void push_back(const DATA & dataItem, const MAP & mapItem, bool updateFinder = false);
    
    /*! Resize */
    void resize(const UInt & n);
    
    /*! Memory allocation */
    void reserve(const UInt & n);
    //@}
        
    
    /*! @name Printout */ //@{
  public:
    template<typename D, typename M>
    friend ostream & operator<<(ostream & f, const pVect<D,M> & V);
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename DATA, typename MAP>
pVect<DATA,MAP>::
pVect() : startupOk(false)
{  
  data     = Teuchos::rcp(new CONTAINER_DATA);
  map      = Teuchos::rcp(new CONTAINER_MAP);
  mapManip = Teuchos::rcp(new MANIPULATOR);
  
  mapManip->setMap(map);
}

template<typename DATA, typename MAP>
pVect<DATA,MAP>::
pVect(const UInt & n) : startupOk(false)
{  
  data     = Teuchos::rcp(new CONTAINER_DATA(n));
  map      = Teuchos::rcp(new CONTAINER_MAP(n));
  mapManip = Teuchos::rcp(new MANIPULATOR);
  
  mapManip->setMap(map);
}

template<typename DATA, typename MAP>
pVect<DATA,MAP>::
pVect(const pVect & Vect) : startupOk(true)
{  
  data     = Teuchos::rcp(new CONTAINER_DATA);
  map      = Teuchos::rcp(new CONTAINER_MAP);
  mapManip = Teuchos::rcp(new MANIPULATOR);
  
  *data = *(Vect.data);
  *map  = *(Vect.map);

  mapManip->setMap(map);
  mapManip->setNormalIndexing();
  mapManip->buildFinder();
}
 
template<typename DATA, typename MAP>
pVect<DATA,MAP>::
pVect(const CONTAINER_MAP & Map, const CONTAINER_DATA & Data) : startupOk(true)
{
  assert(Map.size() == Data.size());
    
  data     = Teuchos::rcp(new CONTAINER_DATA);
  map      = Teuchos::rcp(new CONTAINER_MAP);
  mapManip = Teuchos::rcp(new MANIPULATOR);
  
  *map  = Map;
  *data = Data;
  
  mapManip->setMap(map);
  mapManip->setNormalIndexing();
  mapManip->buildFinder();
}

template<typename DATA, typename MAP>
pVect<DATA,MAP>::
pVect(const Teuchos::RCP<CONTAINER_MAP> & Map, const Teuchos::RCP<CONTAINER_DATA> & Data) : startupOk(true)
{
  assert(Map->size() == Data->size());
    
  data     = Teuchos::rcp(new CONTAINER_DATA);
  map      = Teuchos::rcp(new CONTAINER_MAP);
  mapManip = Teuchos::rcp(new MANIPULATOR);
  
  *map  = *Map;
  *data = *Data;
  
  mapManip->setMap(map);
  mapManip->setNormalIndexing();
  mapManip->buildFinder();
}

template<typename DATA, typename MAP>
pVect<DATA,MAP>
pVect<DATA,MAP>::
operator=(const pVect & Vect)
{
  startupOk = true;
  
  *data = *(Vect.data);
  *map  = *(Vect.map);
   
  mapManip->setMap(map);
  mapManip->setNormalIndexing();
  mapManip->buildFinder();
  
  return(*this);
}



//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
setData(const CONTAINER_MAP & Map, const CONTAINER_DATA & Data)
{
  assert(Map.size() == Data.size());
  
  startupOk = true;
  
  *map  = Map;
  *data = Data;
  
  mapManip->resetFinder();
  mapManip->setNormalIndexing();
  mapManip->buildFinder();
}

template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
setData(const Teuchos::RCP<CONTAINER_MAP> & Map, const Teuchos::RCP<CONTAINER_DATA> & Data)
{
  assert(Map.size() == Data.size());

  startupOk = true;
  
  *map  = Map;
  *data = Data;
  
  mapManip->resetFinder();
  mapManip->setNormalIndexing();
  mapManip->buildFinder();
}

template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
addData(const CONTAINER_MAP & Map, const CONTAINER_DATA & Data)
{
  assert(Map.size() == Data.size());
  
  startupOk = true;
  
  for(UInt i=1; i <= Map.size(); ++i)
  {
    map->push_back(Map(i));
    data->push_back(Data(i));
  }
  
  mapManip->resetFinder();
  mapManip->setNormalIndexing();
  mapManip->buildFinder();
}

template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
addData(const Teuchos::RCP<CONTAINER_MAP> & Map, const Teuchos::RCP<CONTAINER_DATA> & Data)
{
  assert(Map->size() == Data->size());
  
  startupOk = true;
  
  for(UInt i=1; i <= Map->size(); ++i)
  {
    map->push_back(Map(i));
    data->push_back(Data(i));
  }
  
  mapManip->resetFinder();
  mapManip->setNormalIndexing();
  mapManip->buildFinder();
}

template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
resetFinder()
{
  startupOk = false;
  mapManip->resetFinder();
}

template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
updateFinder()
{
  assert(data->size() == map->size());
  startupOk = true;
  
  mapManip->resetFinder();
  mapManip->setNormalIndexing();
  mapManip->buildFinder();
}

template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
clear()
{
  startupOk = false;
  
  mapManip->resetFinder();
  map->clear();
  data->clear();
}

template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
bufferLids()
{
  map->bufferLids();
}

template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
restoreLids()
{
  map->restoreLids();
}



//_________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DATA, typename MAP>
typename pVect<DATA,MAP>::CONTAINER_DATA &
pVect<DATA,MAP>::
getDataRef()
{
  return(*data);
}

template<typename DATA, typename MAP>
typename pVect<DATA,MAP>::CONTAINER_MAP &
pVect<DATA,MAP>::
getMapRef()
{
  return(*map);
}

template<typename DATA, typename MAP>
typename pVect<DATA,MAP>::RCP_CONTAINER_DATA &
pVect<DATA,MAP>::
getDataRcp()
{
  return(data);
}

template<typename DATA, typename MAP>
typename pVect<DATA,MAP>::RCP_CONTAINER_MAP &
pVect<DATA,MAP>::
getMapRcp()
{
  return(map);
}

template<typename DATA, typename MAP>
const typename pVect<DATA,MAP>::RCP_CONTAINER_DATA &
pVect<DATA,MAP>::
getDataRcp() const
{
  return(data);
}

template<typename DATA, typename MAP>
const typename pVect<DATA,MAP>::RCP_CONTAINER_MAP &
pVect<DATA,MAP>::
getMapRcp() const
{
  return(map);
}

template<typename DATA, typename MAP>
const typename pVect<DATA,MAP>::CONTAINER_DATA &
pVect<DATA,MAP>::
getDataRef() const
{
  return(*data);
}

template<typename DATA, typename MAP>
const typename pVect<DATA,MAP>::CONTAINER_MAP &
pVect<DATA,MAP>::
getMapRef() const
{
  return(*map);
}

template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
setMap(const CONTAINER_MAP & RowMap)
{
  *map = RowMap;
}

template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
setData(const CONTAINER_DATA & Data)
{
  *data = Data;
}

    
    
//_________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DATA, typename MAP>
UInt
pVect<DATA,MAP>::
sizeL()
{
  assert(data->size() == map->size());
  return(data->size());
}

template<typename DATA, typename MAP>
DATA &
pVect<DATA,MAP>::
getL(const UInt & lid)
{
  assert(lid>=1);
  assert(lid<=data->size());
  assert(data->size() == map->size());
  
  return(data->get(lid));
}

template<typename DATA, typename MAP>
const DATA &
pVect<DATA,MAP>::
getL(const UInt & lid) const
{
  assert(lid>=1);
  assert(lid<=data->size());
  assert(data->size() == map->size());
  
  return(data->get(lid));
}

template<typename DATA, typename MAP>
DATA &
pVect<DATA,MAP>::
getG(const UInt & gid)
{
  assert(startupOk);
  assert(isG(gid));
  assert(data->size() == map->size());
  
  UInt lid = mapManip->getLidItem(MAP(0,gid));
  return(data->get(lid));
}

template<typename DATA, typename MAP>
const DATA &
pVect<DATA,MAP>::
getG(const UInt & gid) const
{
  assert(startupOk);
  assert(isG(gid));
  assert(data->size() == map->size());
  
  UInt lid = mapManip->getLidItem(MAP(0,gid));
  return(data->get(lid));
}

template<typename DATA, typename MAP>
bool
pVect<DATA,MAP>::
isG(const UInt & gid) const
{
  assert(startupOk);
  assert(data->size() == map->size());
  return( mapManip->isItem(MAP(0,gid)) );
}

template<typename DATA, typename MAP>
MAP & 
pVect<DATA,MAP>::
getMapL(const UInt & lid)
{
  assert(lid>=1);
  assert(lid<=data->size());
  assert(data->size() == map->size());
  
  return(map->get(lid));
}

template<typename DATA, typename MAP>
const MAP &
pVect<DATA,MAP>::
getMapL(const UInt & lid) const
{
  assert(lid>=1);
  assert(lid<=data->size());
  assert(data->size() == map->size());
  
  return(map->get(lid));
}

template<typename DATA, typename MAP>
MAP &
pVect<DATA,MAP>::
getMapG(const UInt & gid)
{
  assert(startupOk);
  assert(isG(gid));
  assert(data->size() == map->size());
  
  UInt lid = mapManip->getLidItem(MAP(0,gid));
  return(map->get(lid));
}

template<typename DATA, typename MAP>
const MAP &
pVect<DATA,MAP>::
getMapG(const UInt & gid) const
{
  assert(startupOk);
  assert(isG(gid));
  assert(data->size() == map->size());
  
  UInt lid = mapManip->getLidItem(MAP(0,gid));
  return(map->get(lid));
}

template<typename DATA, typename MAP>
DATA &
pVect<DATA,MAP>::
getDataL(const UInt & lid)
{
  assert(lid>=1);
  assert(lid<=data->size());
  assert(data->size() == map->size());
  
  return(data->get(lid));
}

template<typename DATA, typename MAP>
const DATA &
pVect<DATA,MAP>::
getDataL(const UInt & lid) const
{
  assert(lid>=1);
  assert(lid<=data->size());
  assert(data->size() == map->size());
  
  return(data->get(lid));
}

template<typename DATA, typename MAP>
DATA &
pVect<DATA,MAP>::
getDataG(const UInt & gid)
{
  assert(startupOk);
  assert(isG(gid));
  assert(data->size() == map->size());
  
  UInt lid = mapManip->getLidItem(MAP(0,gid));
  return(data->get(lid));
}

template<typename DATA, typename MAP>
const DATA &
pVect<DATA,MAP>::
getDataG(const UInt & gid) const
{
  assert(startupOk);
  assert(isG(gid));
  assert(data->size() == map->size());
  
  UInt lid = mapManip->getLidItem(MAP(0,gid));
  return(data->get(lid));
}



//_________________________________________________________________________________________________
// EASY FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DATA, typename MAP>
UInt
pVect<DATA,MAP>::
size() const
{
  assert(data->size() == map->size());
  return(data->size());
}

template<typename DATA, typename MAP>
DATA &
pVect<DATA,MAP>::
get(const UInt & lid)
{
  assert(lid>=1);
  assert(lid<=data->size());
  assert(data->size() == map->size());
  
  return(data->get(lid));
}

template<typename DATA, typename MAP>
const DATA &
pVect<DATA,MAP>::
get(const UInt & lid) const
{
  assert(lid>=1);
  assert(lid<=data->size());
  assert(data->size() == map->size());
  
  return(data->get(lid));
}

template<typename DATA, typename MAP>
DATA &
pVect<DATA,MAP>::
operator() (const UInt & lid)
{
  assert(lid>=1);
  assert(lid<=data->size());
  assert(data->size() == map->size());
  
  return(data->get(lid));
}

template<typename DATA, typename MAP>
const DATA &
pVect<DATA,MAP>::
operator() (const UInt & lid) const
{
  assert(lid>=1);
  assert(lid<=data->size());
  assert(data->size() == map->size());
  
  return(data->get(lid));
}

template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
push_back(const MAP & mapItem, const DATA & dataItem, bool updateFinder)
{
  if(updateFinder)
  {
    startupOk = true;
  
    map->push_back(mapItem);
    data->push_back(dataItem);
  
    mapManip->resetFinder();
    mapManip->setNormalIndexing();
    mapManip->buildFinder();
  }
  else
  {
    startupOk = false;
    
    map->push_back(mapItem);
    data->push_back(dataItem);
  }
}

template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
push_back(const DATA & dataItem, const MAP & mapItem, bool updateFinder)
{
  if(updateFinder)
  {
    startupOk = true;
  
    map->push_back(mapItem);
    data->push_back(dataItem);
  
    mapManip->resetFinder();
    mapManip->setNormalIndexing();
    mapManip->buildFinder();
  }
  else
  {
    startupOk = false;
    
    map->push_back(mapItem);
    data->push_back(dataItem);
  }
}

template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
resize(const UInt & n)
{
  startupOk = false;
  
  data->resize(n);
  map->resize(n);
}

template<typename DATA, typename MAP>
void
pVect<DATA,MAP>::
reserve(const UInt & n)
{
  data->reserve(n);
  map->reserve(n);
}



//_________________________________________________________________________________________________
// PRINTOUT
//-------------------------------------------------------------------------------------------------
template<typename D, typename M>
ostream & operator<<(ostream & f, const pVect<D,M> & V)
{
  for(UInt i=1; i <= V.data->size(); ++i)
  {
    f << i << endl;
    f << " map:  " << V.map->get(i);
    f << " data: " << V.data->get(i) << endl;
  }
  
  return(f);
}

template<typename DATA, typename MAP>
template<class ARK>
void
pVect<DATA,MAP>::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & *map;
  ar & *data;
}

#endif
