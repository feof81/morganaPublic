/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef DOFMAPSTATIC1D_OPTIONS_H
#define DOFMAPSTATIC1D_OPTIONS_H

#include <set>
#include "typesInterface.hpp"
#include "simpleFormats.hpp"


/*! Options for the 1d static mapping class. Contains the list of the active \c geoId */
class dofMapStatic1d_options
{
    /*! @name Internal data */ //@{
  public:
    set<UInt> activeIds;
    //@}
    
    /*! @name Input functions */ //@{
  public:
    dofMapStatic1d_options();
    dofMapStatic1d_options(const dofMapStatic1d_options & O);
    dofMapStatic1d_options operator=(const dofMapStatic1d_options & O);
    void addGeoId(const UInt & Id);
    void clear();
    //@}
    
    /*! @name Output functions */ //@{
  public:
    UInt numGeoIds() const;
    bool isGeoId(const UInt & Id) const;
    sVect<UInt> geoIdsList() const;
    set<UInt>   getIdsSet() const;
    //@}
    
    /*! @name Outstream operator */ //@{
  public:
    friend ostream & operator<<(ostream & f, const dofMapStatic1d_options & P);
    //@}
};

#endif
