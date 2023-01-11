/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef BOUNDCARD_H
#define BOUNDCARD_H

#include <assert.h>
#include <iostream>
#include <fstream>
#include <set>

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "simpleFormats.hpp"


using namespace std;


/*! Card for bounds imposition */
class boundCard
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
  
    /*! @name Internal data */ //@{
  public:
    set<UInt> geoIds;
    Real value;
    UInt flag; 
    //@}
    
    /*! @name Constructors */ //@{
  public:
    boundCard();
    boundCard(const set<UInt> & GeoIds, const UInt & Flag, const Real & Value);
    boundCard(const boundCard & B);
    boundCard operator=(const boundCard & B);
    bool operator!=(const boundCard & B) const;
    bool operator==(const boundCard & B) const;
    //@}
    
    /*! @name Functions */ //@{
  public:
    void setGeoIds(const set<UInt> & GeoIds);
    void addGeoId(const UInt & GeoId);
    void resetGeoIds();
    void setValue(const Real & Value);
    void setFlag(const UInt & Flag);
    bool isGeoId(const UInt & GeoId);
    const set<UInt> & getGeoIds() const;
    sVect<UInt> getVectGeoIds() const;
    Real & getValue();
    const Real & getValue() const;
    UInt & getFlag();
    const UInt & getFlag() const;
    //@}
    
    /*! @name Outstream operator */ //@{
  public:
    friend ostream & operator<<(ostream & f, const boundCard & B);
    //@}
};



//_________________________________________________________________________________________________
// INLINE FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
boundCard::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & geoIds;
  ar & value;
  ar & flag;
}



//_________________________________________________________________________________________________
// FUNCTIONS
//-------------------------------------------------------------------------------------------------
namespace boundCard_functions
{
  sVect<UInt> boundList(const sVect<boundCard> & cards);
  sVect<UInt> unBoundList(const sVect<UInt> mainList, const sVect<boundCard> & cards);
  bool        isCompatible(const sVect<UInt> mainList, const sVect<boundCard> & cards);
}





#endif
