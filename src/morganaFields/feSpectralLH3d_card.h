/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FESPECTRALLH3D_CARD_H
#define FESPECTRALLH3D_CARD_H

#include "glBase.h"
#include "polyDynamic.h"
#include "elCard3d.hpp"
#include "feDynamicDofCard3d.h"
#include "morganaFields.hpp"

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>


/*! Finite element card for the spectral3d-FE */
class feSpectralLH3d_card
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
  
  /*! @name Internal data */ //@{
  public:
    static const FECardLabel cardLabel = CR_SK_3d; 
    bool isActive;
    UInt rx, ry, rz;
    //@}
    
    /*! @name Constructor and functions */ //@{
  public:
    feSpectralLH3d_card();
    feSpectralLH3d_card(const feSpectralLH3d_card & C);
    feSpectralLH3d_card operator=(const feSpectralLH3d_card & C);
    bool operator!=(const feSpectralLH3d_card & C) const;
    //@}
    
    /*! @name Set functions */ //@{
  public:
    void setIsActive(const bool & IsActive);
    void setR(const UInt & RX, const UInt & RY, const UInt & RZ);
    //@}
    
    /*! @name Get const reference */ //@{
  public:
    const bool & getIsActive() const;
    const UInt & getRx() const;
    const UInt & getRy() const;
    const UInt & getRz() const;
    //@}
    
    /*! @name Get reference */ //@{
  public:
    bool & getIsActive();
    UInt & getRx();
    UInt & getRy();
    UInt & getRz();
    //@}
    
    /*! @name Printout */ //@{
  public:
    friend ostream & operator<<(ostream & f, const feSpectralLH3d_card & V);
    //@}
};


//_________________________________________________________________________________________________
// PARALLEL SUPPORT
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
feSpectralLH3d_card::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & isActive;
  ar & rx;
  ar & ry;
  ar & rz;
}

#endif
