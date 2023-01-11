/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef SGRAPHITEMSUBLOC_H
#define SGRAPHITEMSUBLOC_H

#include "morganaTypes.hpp"
#include "simpleFormats.hpp"
#include "pGraphItem.h"

#include <assert.h>
#include <iostream>

using namespace std;


/*! Inherits \c pGraphItem and adds a list of indices. For instance every face is connected to two elements: this data structure allows to store the local indices that the face
has in the two elements. */
class pGraphItemSubLoc : public pGraphItem
{ 
  /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
  
    /*! @name Internal data */ //@{
  public: 
    /*! Lista indici locali  */
    sVect<UInt> subLocId;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    /*! Constructor */
    pGraphItemSubLoc();
    
    /*! Constructor
    \param N number of connected elements */
    pGraphItemSubLoc(const UInt & N);
    
    /*! Constructor
    \param values vector of the connected ids
    \param id element id */
    pGraphItemSubLoc(const sVect<UInt> & values);
    
    /*! Copy constructor */
    pGraphItemSubLoc(const pGraphItemSubLoc & G);
    //@}
    
    
    /*! @name Sizing functions */ //@{
  public:
    /*! Number of connected elements */
    UInt size() const;
    
    /*! Resizing */
    void resize(const UInt & dim);
    
    /*! Memory allocation */
    void reserve(const UInt & n);
    
    /*! Push back */
    void push_back(const UInt & value, const UInt & subIndex);
    //@}
    
    /*! @name Get-Set functions */ //@{
  public:
    /*! Set of the local index
    \param i index
    \param value sublocal index */
    void setSubLocId(const UInt & i, const UInt & value);
    
    /*! Set all the connected ids */
    void setCids(const sVect<UInt> & Connected);
    
    /*! Set cid-subId */
    void setData(const UInt & i, const UInt & Id, const UInt & subId);
    
    /*! Set cid-subId */
    void setData(const sVect<UInt> & Connected, const sVect<UInt> & subIds);
    
    /*! Get sublocal index
    \param i index
    \return subLocal index */
    UInt & getSubLocId(const UInt & i);
    
    /*! Get sublocal index
    \param i index
    \return subLocal index */
    const UInt & getSubLocId(const UInt & i) const;
    //@}
    
    
    /*! @name Operator functions */ //@{
  public:
    /*! Equality operator */
    pGraphItemSubLoc & operator=(const pGraphItemSubLoc & E);
    //@}
    
    
    /*! @name Outstream operator */ //@{
  public:
    size_t memSize() const;
    
    friend ostream & operator<<(ostream & f, const pGraphItemSubLoc & G);
    //@}
};


template<class ARK>
void
pGraphItemSubLoc::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  //Communication
  assert(subLocId.size() == connected.size());
  
  int  N = connected.size();
  ar & N;
  
  if(UInt(N) != connected.size())
  { 
    connected.resize(N);
    subLocId.resize(N);
  }
  
  for(int i=1; i <= N; ++i)
  { 
    ar & connected(i);
    ar & subLocId(i);
  }
}


//_________________________________________________________________________________________________
// ASSOCIATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
size_t dynamicSizeOf(const pGraphItem & A);

#endif
