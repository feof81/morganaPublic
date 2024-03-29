/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PMAPITEMSENDRECV_H
#define PMAPITEMSENDRECV_H

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "pMapItem.h"
#include "morganaTypes.hpp"

#include "morganaPmapItems.hpp"

#include <assert.h>
#include <iostream>


/*! Sending and receiving map.

<b> Stored internal data </b>
<ol>
<li> \c lid local id
<li> \c gid global id
<li> \c pid process id
<li> \c bufLid buffered lid
<li> \c sid sending id
<li> \c rid reading id
</ol>

*/
class pMapItemSendRecv : public pMapItem
{
  /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const pMapItems parallelType = pMapSendRecv;
    UInt sid, rid;
    //@}
    
    /*! @name Constructors and equality operator */ //@{
  public:
    pMapItemSendRecv();
    pMapItemSendRecv(const UInt & Lid, const UInt & Gid);
    pMapItemSendRecv(const UInt & Lid, const UInt & Gid, const UInt & Sid, const UInt & Rid);
    pMapItemSendRecv(const pMapItemSendRecv & M);
    pMapItemSendRecv & operator=(const pMapItemSendRecv & M);
    //@}
    
    /*! @name Set-Get functions */ //@{
  public:
    void setSid(const UInt & Sid);
    void setRid(const UInt & Rid);
    UInt & getSid();
    UInt & getRid();
    const UInt & getSid() const;
    const UInt & getRid() const;
    //@}
    
    /*! @name Logic operators */ //@{
  public:
    bool operator<(const pMapItemSendRecv & E) const;
    bool operator!=(const pMapItemSendRecv & E) const;
    //@}
    
    /*! @name Outstream operators */ //@{
  public:
    friend ostream & operator<<(ostream & f, const pMapItemSendRecv & M);
    
    size_t memSize() const;
    //@}
};



//_________________________________________________________________________________________________
// INLINED AND TEMPLATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
pMapItemSendRecv::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & lid;
  ar & gid;
  ar & pid;
  ar & sid;
  ar & rid;
  ar & bufLid;
}


//_________________________________________________________________________________________________
// ASSOCIATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
size_t dynamicSizeOf(const pMapItemSendRecv & A);


#endif
