/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FEHY2D_CARD_H
#define FEHY2D_CARD_H

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "morganaTypes.hpp"
#include "simpleFormats.hpp"
#include "morganaFields.hpp"


/*! Card for the hybrid FE */
class feHY2d_card
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const FECardLabel cardLabel = CR_HY_2d;
    sVect<bool> activeEdges;
    //@}
    
    /*! @name Constructors and functions */ //@{
  public:
    feHY2d_card();
    feHY2d_card(const feHY2d_card & C);
    feHY2d_card operator=(const feHY2d_card & C);
    bool operator!=(const feHY2d_card & C) const;
    //@}
    
    /*! @name Set and Get functions */ //@{
  public:
    void setActiveEdges(const sVect<bool> & ActiveEdges);
    void setActiveEdge(const UInt & i, const bool & active);
    const sVect<bool> & getActiveEdges() const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    friend ostream & operator<<(ostream & f, const feHY2d_card & V);
    //@}
};



//_________________________________________________________________________________________________
// PARALLEL SUPPORT
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
feHY2d_card::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  //Serialize to Int
  sVect<UInt> buffer(3);
  
  buffer(1) = UInt(activeEdges(1));
  buffer(2) = UInt(activeEdges(2));
  buffer(3) = UInt(activeEdges(3));
  
  //Send and receive
  ar & buffer;
  
  //Transform back
  activeEdges(1) = bool(buffer(1));
  activeEdges(2) = bool(buffer(2));
  activeEdges(3) = bool(buffer(3));
}


#endif
