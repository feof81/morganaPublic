/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef SOCTTREEITEM_H
#define SOCTTREEITEM_H

#include "typesInterface.hpp"
#include "simpleFormats.hpp"


/*! The single leaf of a tree */
class sOctTreeItem
{
    /*! @name Internal data */ //@{
  public:
    UInt father;
    UInt sons[8];
    bool isLeaf;
    UInt id;
    //@}
    
    /*! @name Construcotrs */ //@{
  public:
    sOctTreeItem();
    sOctTreeItem(const UInt & Id);
    sOctTreeItem(const sOctTreeItem & O);
    sOctTreeItem & operator=(const sOctTreeItem & O);
    static UInt octMap(const bool & ix, const bool & iy, const bool & iz);
    //@}
    
    /*! @name Comparison functions */ //@{
  public:
    bool operator<(const sOctTreeItem & V) const;
    bool operator!=(const sOctTreeItem & V) const;
    //@}
    
    /*! @name Set functions */ //@{
  public:
    void setFather(const UInt & Father);
    void setSons(const UInt VVV, const UInt FVV, const UInt VFV, const UInt FFV,
                 const UInt VVF, const UInt FVF, const UInt VFF, const UInt FFF);
    void setSons(const sVect<UInt> Sons);
    void setIsLeaf(const bool & IsLeaf);
    void setId(const UInt & Id);
    //@}
    
    /*! @name Const get Functions */ //@{
  public:
    const UInt & getFather() const;
    const UInt & getSons(const bool & ix, const bool & iy, const bool & iz) const;
    const UInt & getSons(const UInt & i) const;
    const bool & getIsLeaf() const;
    const UInt & getId() const;
    //@}
    
    /*! @name Get functions */ //@{
  public:
    UInt & getFather();
    UInt & getSons(const bool & ix, const bool & iy, const bool & iz);
    UInt & getSons(const UInt & i);
    bool & getIsLeaf();
    UInt & getId();
    //@}
    
    /*! @name Outstream operator */ //@{
  public:
    friend ostream & operator<<(ostream & f, const sOctTreeItem & P);
    
    size_t memSize() const;
    //@}
};

#endif
