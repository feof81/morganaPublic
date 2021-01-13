/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FEDYNAMICDOFCARD1D_H
#define FEDYNAMICDOFCARD1D_H

#include <ostream>

#include "morganaTypes.hpp"
#include "morganaGeometry.hpp"

using namespace std;


/*! Contains all the information neded to define the base of a dynamic 1d finite element */
class feDynamicDofCard1d
{
    /*! @name Internal data */ //@{
  public:
    ReferenceGeometry geoType;
    UInt lev, locId, elId;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    feDynamicDofCard1d();
    feDynamicDofCard1d(const ReferenceGeometry & GeoType, const UInt & Lev, const UInt & LocId);
    feDynamicDofCard1d(const feDynamicDofCard1d & F);
    feDynamicDofCard1d operator=(const feDynamicDofCard1d & F);
    //@}
    
    /*! @name Set functions */ //@{
  public:
    void setGeoType(const ReferenceGeometry & GeoType);
    void setLevel(const UInt & Lev);
    void setLocalId(const UInt & LocId);
    void setLocalElId(const UInt & LocElId);
    //@}
    
    /*! @name Get functions */ //@{
  public:
    ReferenceGeometry & getGeoType();
    UInt & getLevel();
    UInt & getLocalId();
    UInt & getLocalElId();
    //@}
    
    /*! @name Const reference get functions */ //@{
  public:
    const ReferenceGeometry & getGeoType() const;
    const UInt & getLevel() const;
    const UInt & getLocalId() const;
    const UInt & getLocalElId() const;
    //@}
  
    /*! @name Printout */ //@{
  public:
    friend ostream & operator<<(ostream & f, const feDynamicDofCard1d & V);
    //@}
};

#endif
