/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef BCCARDFE_H
#define BCCARDFE_H

#include <assert.h>
#include <iostream>
#include <fstream>

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "simpleFormats.hpp"


/*! Card for boundary conditions imposition */
class bcCardFE
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
 
    /*! @name Internal data */ //@{
  public:
    UInt geoId;
    UInt bcFlag;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    bcCardFE();
    bcCardFE(const UInt & GeoId, const UInt & BcFlag);
    bcCardFE(const bcCardFE & B);
    bcCardFE operator=(const bcCardFE & B);
    //@}
    
    /*! @name Functions */ //@{
  public:
    void setGeoId(const UInt & GeoId);
    void setBcFlag(const UInt & BcFlag);
    const UInt & getGeoId() const;
    const UInt & getBcFlag() const;
    UInt & getGeoId();
    UInt & getBcFlag();
    //@}
    
    /*! @name Ordinal functions */ //@{
  public:
    bool operator<(const bcCardFE & B) const;
    bool operator!=(const bcCardFE & B) const;
    //@}
    
    /*! @name Outstream operator */ //@{
  public:
    friend ostream & operator<<(ostream & f, const bcCardFE & B);
    //@}
};


template<class ARK>
void
bcCardFE::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & geoId;
  ar & bcFlag;
}


#endif
