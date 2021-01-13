/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FEHYB3D_CARD_H
#define FEHYB3D_CARD_H

#include "morganaTypes.hpp"
#include "simpleFormats.hpp"
#include "morganaFields.hpp"


/*! Card for the hybrid FE */
class feHYB3d_card
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const FECardLabel cardLabel = CR_HYB_3d;
    sVect<bool> activeFaces;
    sVect<bool> facesOrientation;
    //@}
    
    /*! @name Constructors and functions */ //@{
  public:
    feHYB3d_card();
    feHYB3d_card(const feHYB3d_card & C);
    feHYB3d_card operator=(const feHYB3d_card & C);
    bool operator!=(const feHYB3d_card & C) const;
    //@}
    
    /*! @name Set and Get functions */ //@{
  public:
    void setActiveFaces(const sVect<bool> & ActiveFaces);
    void setFacesOrientation(const sVect<bool> & FacesOrientation);
    
    void setActiveFace(const UInt & i, const bool & active);
    void setFaceOrientation(const UInt & i, const bool & FaceOrientation);
    
    const sVect<bool> & getActiveFaces() const;
    const sVect<bool> & getFacesOrientation() const;
    
    sVect<bool> & getActiveFaces();
    sVect<bool> & getFacesOrientation();
    
    bool getActiveFaces(const UInt & i) const;
    bool getFaceOrientation(const UInt & i) const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    friend ostream & operator<<(ostream & f, const feHYB3d_card & V);
    //@}
};



//_________________________________________________________________________________________________
// PARALLEL SUPPORT
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
feHYB3d_card::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & activeFaces;
  ar & facesOrientation;  
}

#endif
