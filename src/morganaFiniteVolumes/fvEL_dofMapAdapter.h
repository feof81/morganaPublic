/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FVEL_DOFMAPADAPTER_H
#define FVEL_DOFMAPADAPTER_H

#include "morganaTypes.hpp"
#include "simpleFormats.hpp"


/*! Adapter for the \c dofMapStatic classes */
class fvEL_dofMapAdapter
{
    /*! @name Internal data */ //@{
  public:
    bool dataLoaded;
    sVect<bool> itemIsActive;
    sVect<UInt> newItemLid;
    sVect<UInt> newToOld;
    UInt numNewItems;
    //@}
    
    /*! @name Functions */ //@{
  public:
    fvEL_dofMapAdapter();
    fvEL_dofMapAdapter(const sVect<bool> & ItemIsActive, const sVect<UInt> & NewItemLid);
    void setData(const sVect<bool> & ItemIsActive, const sVect<UInt> & NewItemLid);
    void startup();
    bool getOldIsActive(const UInt & i) const;
    UInt getOldToNew(const UInt & i) const;
    UInt getNewToOld(const UInt & i) const;
    UInt getNumNewItems() const;
    //@}
};

#endif
