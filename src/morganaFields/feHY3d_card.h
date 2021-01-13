/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FEHY3D_CARD_H
#define FEHY3D_CARD_H

#include "morganaTypes.hpp"
#include "simpleFormats.hpp"
#include "morganaFields.hpp"


/*! Card for the hybrid FE */
class feHY3d_card
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const FECardLabel cardLabel = CR_HY_3d;
    sVect<bool> activeFaces;
    //@}
    
    /*! @name Constructors and functions */ //@{
  public:
    feHY3d_card();
    feHY3d_card(const feHY3d_card & C);
    feHY3d_card operator=(const feHY3d_card & C);
    bool operator!=(const feHY3d_card & C) const;
    //@}
    
    /*! @name Set and Get functions */ //@{
  public:
    void setActiveFaces(const sVect<bool> & ActiveFaces);
    void setActiveFace(const UInt & i, const bool & active);
    const sVect<bool> & getActiveFaces() const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    friend ostream & operator<<(ostream & f, const feHY3d_card & V);
    //@}
};



//_________________________________________________________________________________________________
// PARALLEL SUPPORT
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
feHY3d_card::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  //Serialize to Int
  sVect<UInt> buffer(4);
  
  buffer(1) = UInt(activeFaces(1));
  buffer(2) = UInt(activeFaces(2));
  buffer(3) = UInt(activeFaces(3));
  buffer(4) = UInt(activeFaces(4));
  
  //Send and receive
  ar & buffer;
  
  //Transform back
  activeFaces(1) = bool(buffer(1));
  activeFaces(2) = bool(buffer(2));
  activeFaces(3) = bool(buffer(3));
  activeFaces(4) = bool(buffer(4));
}


#endif
