/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef SEARCHBOUNDINGBOX_H
#define SEARCHBOUNDINGBOX_H

#include <set>
#include "typesInterface.hpp"
#include "simpleFormats.hpp"


/*! The single leaf of a tree */
class searchBoundingBox
{
    /*! @name Internal datas */ //@{
  public:
    sVect<UInt> elements;
    point3d Pmin, Pmax;
    //@}
  
    /*! @name Constructors */ //@{
  public:
    searchBoundingBox();
    searchBoundingBox(const searchBoundingBox & B);
    searchBoundingBox & operator=(const searchBoundingBox & B);
    //@}
    
    /*! @name Set functions */ //@{
  public:
    void setElements(const sVect<UInt> & Elements);
    void setBoundingBox(const point3d & PMIN, const point3d & PMAX);
    //@}
    
    /*! @name Get functions */ //@{
  public:
    const sVect<UInt> & getElements() const;
    const point3d     & getPmin() const;
    const point3d     & getPmax() const;
    //@}
    
    /*! @name Partitioning functions */ //@{
  public:
    bool    isInternal(const point3d & P);
    point3d getPmin(const bool & ix, const bool & iy, const bool & iz) const;
    point3d getPmax(const bool & ix, const bool & iy, const bool & iz) const;
    UInt    octPart(const point3d & P);
    void    octPart(bool & ix, bool & iy, bool & iz, const point3d & P);
    std::set<UInt> octPart(const point3d & Pmin, const point3d & Pmax);
    //@}
    
    /*! @name Distance functions */ //@{
  public:
    Real getMinDist(const point3d & P);
    Real getMaxDist(const point3d & P);
    //@}
    
    /*! @name Outstream operator */ //@{
  public:
    friend ostream & operator<<(ostream & f, const searchBoundingBox & B);
    //@}
};

#endif
