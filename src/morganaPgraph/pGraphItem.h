/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PGRAPHITEM_H
#define PGRAPHITEM_H

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "morganaTypes.hpp"
#include "simpleFormats.hpp"

#include <set>
#include <assert.h>
#include <iostream>

using namespace std;


/*! This class represents the connection between an element and a set of elements. A simple example is the connection between an element and its set of nodes.
This class contains an Id and a list of integers representing the connected elements. It can represents a generic graph.
This class can be ordered using the connected ids. This are ordered and each ordered id is compared. This class, therefore, contains an ordered list of ids
that is updated at each operation. The \c push_back can, optionally, disable the authomatic update. An internal flag keep track of the updates of to the 
ordered list. This flag prevents to access to an old ordered list throgh some assert checks.*/
class pGraphItem
{ 
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
  
    /*! @name Internal data */ //@{
  public:
    /*! List of connected ids */
    sVect<UInt> connected;
    
    /*! List of the ordered nodes */
    sVect<UInt> orderedNodes;

    /*! /c orderedNodes are update */
    bool nodesOrdered;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    /*! Constructor */
    pGraphItem();
    
    /*! Constructor
    \param N number of elements */
    pGraphItem(UInt N);
    
    /*! Constructor
    \param values connected elements
    \param id element id */
    pGraphItem(const sVect<UInt> & values);
    
    /*! Copy constructor */
    pGraphItem(const pGraphItem & G);
    //@}
    
    /*! @name Sizing functions */ //@{
  public:
    /*! Return the number of the points */
    UInt size() const;
    
    /*! Resize */
    void resize(const UInt & dim);
    
    /*! Memory allocation */
    void reserve(const UInt & n);
    
    /*! Add an element. The \c update flag forces the rebuild of the ordered list of cids. 
    If set to false the list is not updated.*/
    void push_back(const UInt & value, const bool & update = true);
    
    /*! Add an element only if it is not already present */
    void merge(const pGraphItem & G);
    
    /*! Clear function */
    void clear();
    //@}
    
    /*! @name Get functions */ //@{
  public:
    /*! Return the id i-th node: does not update the ordering */
    UInt & getCid(const UInt & i);
    
    /*! Return the id i-th node */
    const UInt & getCid(const UInt & i) const;
    
    /*! Returns all the connected ids */
    const sVect<UInt> & getCids() const;
    
    /*! Access operator */
    inline UInt & operator()(const UInt & i);
    
    /*! Access operator */
    inline const UInt & operator()(const UInt & i) const;
    //@}
    
    /*! @name Set functions */ //@{
  public:    
    /*! Sets the id of the i-th node: updates the ordered cids list */
    void setCid(const UInt & i, const UInt & value);
    
    /*! Sets all the connected ids: updates the ordered cids list */
    void setCids(const sVect<UInt> & Connected);
    //@}

    /*! @name Node sorting functions */ //@{
  public:
    /*! Updates \c orderedNodes */
    void updateSorting();
    
    /*! Gets the i-th node of \c orderedNodes */
    const UInt & getSorted(const UInt & i) const;
    
    /*! Returns a flag that states whether the cids ordered list is up-to-date */
    const bool & getNodesOrdered() const;
    //@}
  
    /*! @name Operator functions */ //@{
  public:
    /*! The "less" operator: the nodes are ascending-ordered. The first "less" index wins */
    bool operator<(const pGraphItem & E) const;
    
    /*! The inequality operator*/
    bool operator!=(const pGraphItem & E) const;
    
    /*! Equality operator */
    pGraphItem & operator=(const pGraphItem & E);
    //@}
    
    /*! @name Outstream operator */ //@{
  public:
    void   printSorted();
    size_t memSize() const;
    
    friend ostream & operator<<(ostream & f, const pGraphItem & G);
    //@}
};


//_________________________________________________________________________________________________
// INLINE FUNCTIONS
//-------------------------------------------------------------------------------------------------
UInt & 
pGraphItem::
operator()(const UInt & i)
{
  assert(i >= 1);
  assert(i <= connected.size());
  
  nodesOrdered = false;
  
  return(connected(i));
}

const UInt & 
pGraphItem::
operator()(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= connected.size());
  return(connected(i));
}

template<class ARK>
void
pGraphItem::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  int  N = connected.size();
  ar & N;
  
  if(UInt(N) != connected.size())
  { connected.resize(N); }
  
  for(int i=1; i <= N; ++i)
  { ar & connected(i); }
  
  updateSorting();
}


//_________________________________________________________________________________________________
// ASSOCIATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
size_t dynamicSizeOf(const pGraphItem & A);

#endif

