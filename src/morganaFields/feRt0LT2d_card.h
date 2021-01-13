/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FERT0LT2D_CARD_H
#define FERT0LT2D_CARD_H

#include "morganaTypes.hpp"

#include "geoMapInterface.hpp"

#include "morganaFields.hpp"
#include "elCard2d.hpp"
#include "feStaticDofCard2d.h"
#include "feStaticEvalIterators.hpp"

#include "feMapPiola2d.hpp"


//_________________________________________________________________________________________________
// FE CARD
//-------------------------------------------------------------------------------------------------

/*! Raviart Thomas 2d finite element card */
class feRt0LT2d_card
{
  /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}  
  
    /*! @name Internal data */ //@{
  public:
    static const FECardLabel cardLabel = CR_RT0LT_2d;
    bool loaded;
    sVect<bool> edgesOrientation;
    //@}
  
    /*! @name Constructors */ //@{
  public:
    feRt0LT2d_card();
    feRt0LT2d_card(const sVect<bool> & EdgesOrientation);
    feRt0LT2d_card(const feRt0LT2d_card & card);
    feRt0LT2d_card operator=(const feRt0LT2d_card & card);
    //@}
    
    /*! @name Set and get */ //@{
  public:
    void setEdgesOrientation(const sVect<bool> & EdgesOrientation);
    void setEdgeOrientation(const UInt & i, const bool & EdgeOrientation);
    const sVect<bool> & getEdgesOrientation() const;
    sVect<bool> & getEdgesOrientation();
    bool getEdgeOrientation(const UInt & i) const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    friend ostream & operator<<(ostream & f, const feRt0LT2d_card & V);
    //@}
};


#endif
