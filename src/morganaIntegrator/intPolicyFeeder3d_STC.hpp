/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTPOLICYFEEDER3D_STC_HPP
#define INTPOLICYFEEDER3D_STC_HPP

#include <set>
#include "intPolicySTC.h"
#include "mesh3d.hpp"
#include "simpleFormats.hpp"


/*! Feeder of the integration cards */
template<typename GEOSHAPE, typename PMAPTYPE>
class intPolicyFeeder3d_STC
{
    /*! @name Typedefs */ //@{
  public:
    typedef intPolicySTC                       INTCARD;
    typedef mesh3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> MESH3D;
    typedef pVect<INTCARD,PMAPTYPE>            INTCARDS;
    //@}
   
   /*! @name Data */ //@{
 public:
   bool gridLoaded;
   Teuchos::RCP<const MESH3D> grid3d;
   //@}
   
   /*! @name Functions */ //@{
 public:
   intPolicyFeeder3d_STC();
   intPolicyFeeder3d_STC(const Teuchos::RCP<const MESH3D> & Grid3d);
   void setMesh(const Teuchos::RCP<const MESH3D> & Grid3d);
   INTCARDS getIntCards(const set<UInt> geoIds);
   INTCARDS getIntCards(const sVect<UInt> geoIds);
   //@}
};


template<typename GEOSHAPE, typename PMAPTYPE>
intPolicyFeeder3d_STC<GEOSHAPE,PMAPTYPE>::
intPolicyFeeder3d_STC()
{
  gridLoaded = false;
}

template<typename GEOSHAPE, typename PMAPTYPE>
intPolicyFeeder3d_STC<GEOSHAPE,PMAPTYPE>::
intPolicyFeeder3d_STC(const Teuchos::RCP<const MESH3D> & Grid3d)
{
  gridLoaded = true;
  grid3d     = Grid3d;
}

template<typename GEOSHAPE, typename PMAPTYPE>
void
intPolicyFeeder3d_STC<GEOSHAPE,PMAPTYPE>::
setMesh(const Teuchos::RCP<const MESH3D> & Grid3d)
{
  gridLoaded = true;
  grid3d     = Grid3d;
}

template<typename GEOSHAPE, typename PMAPTYPE>
typename intPolicyFeeder3d_STC<GEOSHAPE,PMAPTYPE>::INTCARDS
intPolicyFeeder3d_STC<GEOSHAPE,PMAPTYPE>::
getIntCards(const set<UInt> geoIds)
{
  assert(gridLoaded);
  
  UInt geoId;
  bool active;
  INTCARDS outVect(grid3d->getNumElements());
  outVect.setMap(grid3d->getElements().getRowMap());
  
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    geoId  = grid3d->getElementL(i).getGeoId();
    active = (geoIds.count(geoId) >= 1);
    outVect(i).setIsActive(active);
  }
  
  return(outVect);
}


template<typename GEOSHAPE, typename PMAPTYPE>
typename intPolicyFeeder3d_STC<GEOSHAPE,PMAPTYPE>::INTCARDS
intPolicyFeeder3d_STC<GEOSHAPE,PMAPTYPE>::
getIntCards(const sVect<UInt> geoIds)
{
  assert(gridLoaded);
  
  set<UInt> geoSet;
  
  for(UInt i=1; i <= geoIds.size(); ++i)
  { geoSet.insert(geoIds(i)); }
  
  return(getIntCards(geoSet));
}


#endif
