/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef TRAITSRESUMEHDF5_H
#define TRAITSRESUMEHDF5_H

#include "pMapItem.h"
#include "pMapItemShare.h"
#include "pMapItemSendRecv.h"

#include "point2d.h"
#include "point3d.h"
#include "tensor2d.h"
#include "tensor3d.h"
#include "staticVector.hpp"

#include "geoElement.hpp"


/*! Trait class for resuming computations */
template<typename DATA>
class traitsResumeHDF5
{ };


/*! Trait class for resuming computations: \c pMapItem */
template<>
class traitsResumeHDF5<pMapItem>
{
  public:
    typedef pMapItem  DATA;
    typedef Real   OUTTYPE;
  
  public:
    static UInt    size();
    static OUTTYPE getValue(const DATA & data, const UInt & i);
    static void    setValue(const OUTTYPE & item, const UInt & i, DATA & data);
};


/*! Trait class for resuming computations: \c pMapItemShare */
template<>
class traitsResumeHDF5<pMapItemShare>
{
  public:
    typedef pMapItemShare  DATA;
    typedef Real        OUTTYPE;
  
  public:
    static UInt    size();
    static OUTTYPE getValue(const DATA & data, const UInt & i);
    static void    setValue(const OUTTYPE & item, const UInt & i, DATA & data);
};


/*! Trait class for resuming computations: \c pMapItemSendRecv */
template<>
class traitsResumeHDF5<pMapItemSendRecv>
{
  public:
    typedef pMapItemSendRecv  DATA;
    typedef Real           OUTTYPE;
  
  public:
    static UInt    size();
    static OUTTYPE getValue(const DATA & data, const UInt & i);
    static void    setValue(const OUTTYPE & item, const UInt & i, DATA & data);
};


/*! Trait class for resuming computations: \c UInt */
template<>
class traitsResumeHDF5<UInt>
{
  public:
    typedef UInt    DATA;
    typedef Real OUTTYPE;
  
  public:
    static UInt    size();
    static OUTTYPE getValue(const DATA & data, const UInt & i);
    static void    setValue(const OUTTYPE & item, const UInt & i, DATA & data);
};


/*! Trait class for resuming computations: \c Real */
template<>
class traitsResumeHDF5<Real>
{
  public:
    typedef Real    DATA;
    typedef Real OUTTYPE;
  
  public:
    static UInt    size();
    static OUTTYPE getValue(const DATA & data, const UInt & i);
    static void    setValue(const OUTTYPE & item, const UInt & i, DATA & data);
};


/*! Trait class for resuming computations: \c point2d */
template<>
class traitsResumeHDF5<point2d>
{
  public:
    typedef point2d  DATA;
    typedef Real  OUTTYPE;
  
  public:
    static UInt    size();
    static OUTTYPE getValue(const DATA & data, const UInt & i);
    static void    setValue(const OUTTYPE & item, const UInt & i, DATA & data);
};


/*! Trait class for resuming computations: \c point3d */
template<>
class traitsResumeHDF5<point3d>
{
  public:
    typedef point3d  DATA;
    typedef Real  OUTTYPE;
  
  public:
    static UInt    size();
    static OUTTYPE getValue(const DATA & data, const UInt & i);
    static void    setValue(const OUTTYPE & item, const UInt & i, DATA & data);
};


/*! Trait class for resuming computations: \c tensor2d */
template<>
class traitsResumeHDF5<tensor2d>
{
  public:
    typedef tensor2d  DATA;
    typedef Real   OUTTYPE;
  
  public:
    static UInt    size();
    static OUTTYPE getValue(const DATA & data, const UInt & i);
    static void    setValue(const OUTTYPE & item, const UInt & i, DATA & data);
};


/*! Trait class for resuming computations: \c tensor3d */
template<>
class traitsResumeHDF5<tensor3d>
{
  public:
    typedef tensor3d  DATA;
    typedef Real   OUTTYPE;
  
  public:
    static UInt    size();
    static OUTTYPE getValue(const DATA & data, const UInt & i);
    static void    setValue(const OUTTYPE & item, const UInt & i, DATA & data);
};


/*! Trait class for resuming computations: \c staticVector */
template<size_t N>
class traitsResumeHDF5<staticVector<N> >
{
  public:
    typedef staticVector<N>  DATA;
    typedef Real          OUTTYPE;
  
  public:
    static UInt    size();
    static OUTTYPE getValue(const DATA & data, const UInt & i);
    static void    setValue(const OUTTYPE & item, const UInt & i, DATA & data);
};


/*! Trait class for resuming computations: \c geoElement */
template<typename GEOSHAPE>
class traitsResumeHDF5<geoElement<GEOSHAPE> >
{
  public:
    typedef geoElement<GEOSHAPE>  DATA;
    typedef Real               OUTTYPE;
  
  public:
    static UInt    size();
    static OUTTYPE getValue(const DATA & data, const UInt & i);
    static void    setValue(const OUTTYPE & item, const UInt & i, DATA & data);
};



//_________________________________________________________________________________________________
// STATIC VECTOR
//-------------------------------------------------------------------------------------------------
template<size_t N>
UInt
traitsResumeHDF5<staticVector<N> >::
size()
{
  return(N);
}

template<size_t N>
typename traitsResumeHDF5<staticVector<N> >::OUTTYPE
traitsResumeHDF5<staticVector<N> >::
getValue(const DATA & data, const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  return(data(i));
}

template<size_t N>
void
traitsResumeHDF5<staticVector<N> >::
setValue(const OUTTYPE & item, const UInt & i, DATA & data)
{
  assert(i >= 1);
  assert(i <= N);
  
  data(i) = item;
}



//_________________________________________________________________________________________________
// GEOELEMENT
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE>
UInt
traitsResumeHDF5<geoElement<GEOSHAPE> >::
size()
{
  return(geoElement<GEOSHAPE>::numPoints + 1);
}

template<typename GEOSHAPE>
typename traitsResumeHDF5<geoElement<GEOSHAPE> >::OUTTYPE
traitsResumeHDF5<geoElement<GEOSHAPE> >::
getValue(const DATA & data, const UInt & i)
{
  assert(i >= 1);
  assert(i <= (geoElement<GEOSHAPE>::numPoints + 1));
  
  if(i <= geoElement<GEOSHAPE>::numPoints)
  { return(data.getCid(i)); }
  else
  { return(data.getGeoId()); }  
}

template<typename GEOSHAPE>
void
traitsResumeHDF5<geoElement<GEOSHAPE> >::
setValue(const OUTTYPE & item, const UInt & i, DATA & data)
{
  assert(i >= 1);
  assert(i <= (geoElement<GEOSHAPE>::numPoints + 1));
  
  if(i <= geoElement<GEOSHAPE>::numPoints)
  { data.setCid(i,item); }
  else
  { data.setGeoId(item); }  
}


#endif
