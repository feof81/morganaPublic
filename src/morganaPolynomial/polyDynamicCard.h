/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef POLYDYNAMICCARD_H
#define POLYDYNAMICCARD_H

#include "typesInterface.hpp"
#include "simpleFormats.hpp"
#include "polyDynamicSubCard.h"


/*! The container class that storage all the information to define a polynomial. */
class polyDynamicCard
{
    /*! @name Internal data */ //@{
  public:
    sVect<polyDynamicSubCard> subCards;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    polyDynamicCard();
    polyDynamicCard(const polyDynamicCard & card);
    polyDynamicCard operator=(const polyDynamicCard & card);
    void operator*=(const polyDynamicCard & card);
    void operator*=(const Real & val);
    void reorder();
    void simplify();
    //@}
    
    /*! @name Set and get functions */ //@{
  public:
    void addSlot(const UInt & cx, const UInt & cy, const UInt & cz, const Real & c);
    void setSlot(const UInt & i, const UInt & cx, const UInt & cy, const UInt & cz, const Real & c);
    UInt getCx(const UInt & i) const;
    UInt getCy(const UInt & i) const;
    UInt getCz(const UInt & i) const;
    Real getCoeff(const UInt & i) const;
    UInt getNumCoeff() const;
    //@}
    
    /*! @name Outstream operator */ //@{
  public:
    friend ostream & operator<<(ostream & f, const polyDynamicCard & G);
    //@}
};

#endif
