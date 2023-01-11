/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PMAPITEMSHARE_H
#define PMAPITEMSHARE_H

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "pMapItem.h"
#include "morganaTypes.hpp"

#include "morganaPmapItems.hpp"

#include <assert.h>
#include <iostream>


/*! Local-Global map item with ownership control. Inherits \c pMapItem and adds the \c shared and \c owned data.

<b> Stored internal data </b>
<ol>
<li> \c lid local id
<li> \c gid global id
<li> \c pid process id
<li> \c bufLid buffered lid
<li> \c shared sharing bool
<li> \c owned owning bool
</ol>

*/
class pMapItemShare : public pMapItem
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const pMapItems parallelType = pMapShare;
    bool shared, owned;
    //@}
    
    /*! @name Constructors and equality operator */ //@{
  public:
    pMapItemShare();
    pMapItemShare(const UInt & Lid, const UInt & Gid);
    pMapItemShare(const UInt & Lid, const UInt & Gid, const bool & Shared, const bool & Owned);
    pMapItemShare(const UInt & Lid, const UInt & Gid, const UInt & Pid, const bool & Shared, const bool & Owned);
    pMapItemShare(const pMapItemShare & M);
    pMapItemShare & operator=(const pMapItemShare & M);
    //@}
    
    /*! @name Set-Get functions */ //@{
  public:
    void setShared(const bool & Shared);
    void setOwned(const bool & Owned);
    bool & getShared();
    bool & getOwned();
    const bool & getShared() const;
    const bool & getOwned() const;
    //@}
    
    /*! @name Outstream operators */ //@{
  public:
    friend ostream & operator<<(ostream & f, const pMapItemShare & M);
    
    size_t memSize() const;
    //@}
};



//_________________________________________________________________________________________________
// INLINED AND TEMPLATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
pMapItemShare::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & lid;
  ar & gid;
  ar & pid;
  ar & shared;
  ar & owned;
  ar & bufLid;
}


//_________________________________________________________________________________________________
// ASSOCIATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
size_t dynamicSizeOf(const pMapItemShare & A);


#endif
