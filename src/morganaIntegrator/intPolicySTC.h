/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTPOLICYSTC_H
#define INTPOLICYSTC_H

#include "elCard3d.hpp"
#include "geoMapInterface.hpp"
#include "morganaIntegrator.hpp"

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>


/*! Integration card STC */
class intPolicySTC
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
  
    /*! @name Internal data */ //@{
  public:
    bool isActive;
    //@}
    
    /*! @name Functions */ //@{
  public:
    intPolicySTC();
    intPolicySTC(const bool & IsActive);
    intPolicySTC(const intPolicySTC & Policy);
    intPolicySTC operator=(const intPolicySTC & Policy);
    void setIsActive(const bool & IsActive);
    const bool & getIsActive() const;
    bool & getIsActive();
    friend ostream & operator<<(ostream & f, const intPolicySTC & P);
    //@}
};


template<class ARK>
void
intPolicySTC::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & isActive;
}


#endif
