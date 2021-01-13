/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FESTATICDOFCARD2D_H
#define FESTATICDOFCARD2D_H

#include <ostream>

#include "morganaTypes.hpp"
#include "morganaGeometry.hpp"

using namespace std;


/*! Information dof-card for the static 2d fields */
class feStaticDofCard2d
{
    /*! @name Internal data */ //@{
  public:
    ReferenceGeometry geoType;
    UInt lev, locId, elId;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    feStaticDofCard2d();
    feStaticDofCard2d(const ReferenceGeometry & GeoType, const UInt & Lev, const UInt & LocId);
    feStaticDofCard2d(const feStaticDofCard2d & F);
    feStaticDofCard2d operator=(const feStaticDofCard2d & F);
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
    friend ostream & operator<<(ostream & f, const feStaticDofCard2d & V);
    //@}
};

#endif
