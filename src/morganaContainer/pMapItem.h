/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PMAPITEM_H
#define PMAPITEM_H

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "morganaPmapItems.hpp"
#include "morganaTypes.hpp"

#include <assert.h>
#include <iostream>

using namespace std;

/*! Local-Global map item

<b> Stored internal data </b>
<ol>
<li> \c lid local id
<li> \c gid global id
<li> \c pid process id
<li> \c bufLid buffered lid
</ol>

Since for serach purposes the \c lid is often authomatically changed, when the \c bufferLid method is called the \c lid in stored in \c bufLid 
*/
class pMapItem
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    
    /*! @name Internal data */ //@{
  public:
    static const pMapItems parallelType = pMapPlain;
    UInt lid, gid, pid, bufLid;
    //@}
    
    
    /*! @name Constructors */ //@{
  public:
    /*! Constructor */
    pMapItem();
    
    /*! Constructor
    \param Lid local id
    \param Gid global id*/
    pMapItem(const UInt & Lid, const UInt & Gid);
    
    /*! Constructor
    \param Lid local id
    \param Gid global id
    \param Pid process id*/
    pMapItem(const UInt & Lid, const UInt & Gid, const UInt & Pid);
    
    /*! Copy constructor */
    pMapItem(const pMapItem & M);
    //@}
    
    
    /*! @name Logic operators */ //@{
  public:
    /*! Equality operator */
    pMapItem & operator=(const pMapItem & M);
    
    /*! Less operator */
    bool operator<(const pMapItem & E) const;
    
    /*! Not equal operator */
    bool operator!=(const pMapItem & E) const;
    //@}
    
    
    /*! @name Get-Set functions */ //@{
  public:
    void setLid(const UInt & Lid);
    void setGid(const UInt & Gid);
    void setPid(const UInt & Pid);
    void setBufLid(const UInt & BufLid);
    UInt & getLid();
    UInt & getGid();
    UInt & getPid();
    UInt & getBufLid();
    const UInt & getLid() const;
    const UInt & getGid() const;
    const UInt & getPid() const;
    const UInt & getBufLid() const;
    //@}
    
    
    /*! @name Other functions */ //@{
  public:
    /*! \c lid -> \c bufLid */
    void bufferLid();
    
    /*! \c bufLid -> \c lid */
    void restoreLid();
    //@}
    
    
    /*! @name Outstream operators */ //@{
  public:
    friend ostream & operator<<(ostream & f, const pMapItem & M);
    
    size_t memSize() const;
    //@}
};



//_________________________________________________________________________________________________
// INLINED AND TEMPLATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
pMapItem::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & lid;
  ar & gid;
  ar & pid;
  ar & bufLid; 
}


//_________________________________________________________________________________________________
// ASSOCIATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
size_t dynamicSizeOf(const pMapItem & A);

#endif
