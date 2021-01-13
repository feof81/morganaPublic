/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef DOFMAPSTATIC2D_OPTIONS_H
#define DOFMAPSTATIC2D_OPTIONS_H

#include <set>
#include "typesInterface.hpp"
#include "simpleFormats.hpp"


/*! Options for the 2d static mapping class. Contains the list of the active geoIds */
class dofMapStatic2d_options
{
    /*! @name Internal data */ //@{
  public:
    set<UInt> activeIds;
    //@}
    
    /*! @name Input functions */ //@{
  public:
    dofMapStatic2d_options();
    dofMapStatic2d_options(const dofMapStatic2d_options & O);
    dofMapStatic2d_options operator=(const dofMapStatic2d_options & O);
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
    friend ostream & operator<<(ostream & f, const dofMapStatic2d_options & P);
    //@}
};

#endif
