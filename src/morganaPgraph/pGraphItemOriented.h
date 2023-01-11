/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PGRAPHITEMORIENTED_H
#define PGRAPHITEMORIENTED_H

#include "morganaTypes.hpp"
#include "simpleFormats.hpp"
#include "pGraphItem.h"

#include <assert.h>
#include <iostream>

using namespace std;


/*! Inherits \c pGraphItem and adds for each connected element a bool flag. This usually represents the connection orientation. For instance in the 
node-edge connection each edge could point inward or otward the node.*/
class pGraphItemOriented : public pGraphItem
{ 
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    
    /*! @name Internal data */ //@{
  public: 
    sVect<bool> orientation;  /*! Orientation vector */
    //@}
    
    
    /*! @name Constructors */ //@{
  public:
    /*! Constructor */
    pGraphItemOriented();
    
    /*! Constructor
    \param N number of the connected elements */
    pGraphItemOriented(const UInt & N);
    
    /*! Constructor
    \param values the ids of the connected elements
    \param id element id */
    pGraphItemOriented(const sVect<UInt> & values);
    
    /*! Copy constructor */
    pGraphItemOriented(const pGraphItemOriented & G);
    //@}
    
    
    /*! @name Sizing functions */ //@{
  public:
    /*! Number of connected elements */
    UInt size() const;
    
    /*! Resize */
    void resize(const UInt & dim);
    
    /*! Memory allocation */
    void reserve(const UInt & n);
    
    /*! Push back
    \param value connected element id
    \param orient orientation */
    void push_back(const UInt & value, bool orient = true);
    //@}
    
    
    /*! @name Get-Set functions */ //@{
  public:
    /*! Orientation set
    \param i sublocal id
    \param value orientation */
    void setOrientation(const UInt & i, const bool & value);
    
    /*! Set all the connected ids */
    void setCids(const sVect<UInt> & Connected);
    
    /*! Set cid-orientation */
    void setData(const UInt & i, const UInt & Id, const bool & logic);
    
    /*! Set cid-orientation */
    void setData(const sVect<UInt> & Connected, const sVect<bool> & logics);
    
    /*! Get orientation
    \param i sublocal id
    \return orientation */
    vector<bool>::reference getOrientation(const UInt & i);
    
    /*! Get orientation
    \param i sublocal id
    \return orientation */
    vector<bool>::const_reference getOrientation(const UInt & i) const;
    //@}
    
    
    /*! @name Operator functions */ //@{
  public:
    /*! Equality operator */
    pGraphItemOriented & operator=(const pGraphItemOriented & E);
    //@}
    
    
    /*! @name Outstream operator */ //@{
  public:
    size_t memSize() const;
    
    friend ostream & operator<<(ostream & f, const pGraphItemOriented & G);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
pGraphItemOriented::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  //Serialize to Int
  sVect<UInt> buffer(orientation.size());
  
  for(UInt i=1; i <= orientation.size(); ++i)
  { buffer(i) = UInt(orientation(i)); }
  
  //Communication
  assert(orientation.size() == connected.size());
  
  int  N = connected.size();
  ar & N;
  
  if(UInt(N) != connected.size())
  { 
    connected.resize(N);
    orientation.resize(N);
    buffer.resize(N);
  }
  
  for(int i=1; i <= N; ++i)
  { 
    ar & connected(i);
    ar & buffer(i);
  }
  
  //Pushing back to bool
  for(UInt i=1; i <= orientation.size(); ++i)
  { orientation(i) = bool(buffer(i)); }
}


//_________________________________________________________________________________________________
// ASSOCIATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
size_t dynamicSizeOf(const pGraphItem & A);


#endif
