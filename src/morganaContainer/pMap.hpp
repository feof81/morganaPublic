/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef PMAP_HPP
#define PMAP_HPP

#include <assert.h>
#include <iostream>

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "simpleFormats.hpp"
#include "morganaTypes.hpp"

using namespace std;


/*! Mapping the local-global correspondence */
template<typename ITEM> class pMap
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
  
  public:
    /*! Items of the graph */
    sVect<ITEM> items;
  
    
    /*! @name Constructors and equality operator */ //@{
  public:
    /*! Constructor */
    pMap();
    
    /*! Constructor -> initializes to normal indexing
    \param n length of the graph */
    pMap(const UInt & n);
    
    /*! Copy constructor */
    pMap(const pMap & G);
    
    /*! Equality operator */
    pMap & operator=(const pMap & G);
    
    /*! Data clean */
    void clear();
    //@}
  
    
    /*! @name Access functions */ //@{
  public:
    /*! Resize */
    void resize(const UInt & n);
    
    /*! Size */
    UInt size() const;
    
    /*! Get - memory cheched */
    ITEM & get(UInt i);
    
    /*! Get - memory cheched */
    const ITEM & get(UInt i) const;
    
    /*! Get - NOT memory cheched */
    inline ITEM & operator() (UInt const i);
    
    /*! Get - NOT memory cheched */
    inline const ITEM & operator() (UInt const i) const;
    
    /*! Set */
    void set(const UInt & i, const ITEM & I);
    
    /*! Push Back */
    void push_back(const ITEM & I);
    
    /*! Memory allocation */
    void reserve(const UInt & n);
    //@}
    
    
    /*! @name Other functions */ //@{
  public:
    /*! For each row sets \c lid -> \c bufLid */
    void bufferLids();
    
    /*! For each row sets \c bufLid -> \c lid */
    void restoreLids();
    //@}
    
    
    /*! @name Outstream operators */ //@{
  public:
    template<typename T>
    friend ostream & operator<<(ostream & f, const pMap<T> & M);
    //@}
};


//_______________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------------
template<typename ITEM>
pMap<ITEM>::
pMap()
{
}

template<typename ITEM>
pMap<ITEM>::
pMap(const UInt & n)
{
  items.resize(n);
  
  for(UInt i=1; i<=items.size(); ++i)
  { items(i).setLid(i); }
}

template<typename ITEM>
pMap<ITEM>::
pMap(const pMap & G) : items(G.items)
{
}

template<typename ITEM>
pMap<ITEM> & 
pMap<ITEM>::
operator=(const pMap & G)
{
  items = G.items;
  return *this;
}

template<typename ITEM>
void
pMap<ITEM>::
clear()
{
  items.clear();
}

template<typename ITEM>
void
pMap<ITEM>::
resize(const UInt & n)
{
  items.resize(n);
}

template<typename ITEM>
UInt
pMap<ITEM>::
size() const
{
  return(items.size());
}

template<typename ITEM>
ITEM &
pMap<ITEM>::
get(UInt i)
{
  assert(i >= 1);
  assert(i <= items.size());
  return(items(i));
}

template<typename ITEM>
const ITEM &
pMap<ITEM>::
get(UInt i) const
{
  assert(i >= 1);
  assert(i <= items.size());
  return(items(i));
}

template<typename ITEM>
ITEM &
pMap<ITEM>::
operator() (UInt const i)
{
  assert(i >= 1);
  assert(i <= items.size());
  return(items(i));
}

template<typename ITEM>
const ITEM &
pMap<ITEM>::
operator() (UInt const i) const
{
  assert(i >= 1);
  assert(i <= items.size());
  return(items(i));
}

template<typename ITEM>
void
pMap<ITEM>::
set(const UInt & i, const ITEM & I)
{
  assert(i >= 1);
  assert(i <= items.size());
  items(i) = I;
}

template<typename ITEM>
void
pMap<ITEM>::
push_back(const ITEM & I)
{
  items.push_back(I);
}

template<typename ITEM>
void
pMap<ITEM>::
reserve(const UInt & n)
{
  items.reserve(n);
}

template<typename T>
ostream & operator<<(ostream & f, const pMap<T> & M)
{
  for(UInt i=1; i <= M.size(); ++i)
  {
    f << M(i);
  }
  
  return(f);
}

template<typename ITEM>
void
pMap<ITEM>::
bufferLids()
{
  for(UInt i=1; i <= items.size(); ++i)
  {
    items(i).bufferLid();
  }
}

template<typename ITEM>
void
pMap<ITEM>::
restoreLids()
{
  for(UInt i=1; i <= items.size(); ++i)
  {
    items(i).restoreLid();
  }
}

template<typename ITEM>
template<class ARK>
void
pMap<ITEM>::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  UInt n = items.size();
  ar & n;
  
  if(n != items.size())
  { this->resize(n); }
  
  for(UInt i=1; i <= items.size(); ++i)
  {
    items(i).serialize(ar,version);
  }
}

#endif
