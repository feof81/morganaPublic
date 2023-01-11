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
    typedef set<UInt> SET;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    SET        activeIds;
    sVect<SET> blockIds;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    dofMapStatic2d_options();
    dofMapStatic2d_options(const dofMapStatic2d_options & O);
    dofMapStatic2d_options operator=(const dofMapStatic2d_options & O);
    void clear();
    //@}
    
    /*! @name GeoId Functions */ //@{
  public:
    void addGeoId(const UInt & Id);
    UInt numGeoIds() const;
    bool isGeoId(const UInt & Id) const;
    sVect<UInt> geoIdsList() const;
    set<UInt>   getIdsSet() const;
    //@}
    
    /*! @name BlockIds Functions */ //@{
  public:
    void setBlockNum(const UInt & numBlocks);
    void addBlockGeoId(const SET & Ids);
    void addBlockGeoId(const UInt & k, const UInt & Id);
    void setBlockGeoId(const UInt & k, const SET         & Ids);
    void setBlockGeoId(const UInt & k, const sVect<UInt> & Ids);
    bool isBlockGeoId(const UInt & k, const UInt & Id) const;
    sVect<UInt> getBlockIdsList(const UInt & k) const;
    set<UInt>   getBlockIdsSet(const UInt & k) const;
    UInt        getNumBlocks() const;
    //@}
    
    /*! @name Outstream operator */ //@{
  public:
    friend ostream & operator<<(ostream & f, const dofMapStatic2d_options & P);
    //@}
};

#endif
