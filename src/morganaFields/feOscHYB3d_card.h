/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FEOSCHYB3D_CARD_H
#define FEOSCHYB3D_CARD_H

#include "feOsc3d_card.h"
#include "feHYB3d_card.h"
#include "morganaFields.hpp"


/*! Card for the primal hybrid oscillating FE */
class feOscHYB3d_card : public feHYB3d_card, public feOsc3d_card
{
  /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const FECardLabel cardLabel = CR_OSHYB_3d;
    //@}
    
    /*! @name Constructors and functions */ //@{
  public:
    feOscHYB3d_card();
    feOscHYB3d_card(const UInt & N);
    feOscHYB3d_card(const feOscHYB3d_card & C);
    feOscHYB3d_card operator=(const feOscHYB3d_card & C);
    bool operator!=(const feOscHYB3d_card & C) const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    friend ostream & operator<<(ostream & f, const feOscHYB3d_card & V);
    //@}
};


//_________________________________________________________________________________________________
// PARALLEL SUPPORT
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
feOscHYB3d_card::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & feHYB3d_card::activeFaces;
  ar & feHYB3d_card::facesOrientation;
  
  int N = this->size();
  
  ar & N;
  
  if(N != this->size())
  { this->resize(N); }
  
  for(UInt i=1; i <= N; ++i)
  {
    ar & this->getY(i);
    ar & this->getH(i);
  }
}

#endif
