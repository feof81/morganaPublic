/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef GEOELEMENT_HPP
#define GEOELEMENT_HPP

#include "pGraphItem.h"
#include "geoShapes.h"

using namespace std;


/*! Implementation of \c pGraphItem that represent the topological connections between an element and the related nodes. */
template<typename GEOSHAPE>
class geoElement : public pGraphItem, public GEOSHAPE
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
  
    /*! @name Internal data */ //@{
  public:
    /*! Static data */
    static const UInt numVertices = GEOSHAPE::numVertices;
    static const UInt numPoints   = GEOSHAPE::numPoints;
    
    /*! Internal logic */
    bool fixedLength;
      
    /*! The geometrical id */
    UInt geoId;
    //@}
  
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    geoElement(bool FixedLength = false);
    
    /*! Constructor */
    geoElement(const UInt & N);
      
    /*! Costruttore di copia */
    geoElement<GEOSHAPE> (const geoElement<GEOSHAPE> & E);
      
    /*! Operatore uguaglianza */
    geoElement<GEOSHAPE> operator=(const geoElement<GEOSHAPE> & E);
      
    /*! Return the geometrical Id */
    const UInt & getGeoId() const;
    
    /*! Return the geometrical Id */
    UInt & getGeoId();
      
    /*! Set id geometrico */
    void setGeoId(const UInt & Id);
    //@}
    
    /*! @name Sizing functions */ //@{
  public:    
    /*! Resize */
    void resize(const UInt & dim);
    
    /*! Memory allocation */
    void reserve(const UInt & n);
    
    /*! Add and element */
    void push_back(const UInt & value, const bool & update = true);
    
    /*! Add and element only if it is not already present */
    void merge(const pGraphItem & G);
    //@}
    
    /*! @name Outstream operator */ //@{
  public:
    template<typename T>
    friend ostream & operator<<(ostream & f, const geoElement<T> & G);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE>
geoElement<GEOSHAPE>::
geoElement(bool FixedLength) : pGraphItem()
{
  fixedLength = FixedLength;
  
  if(FixedLength)
  {
    connected.resize(numPoints);
    updateSorting();
  }
}

template<typename GEOSHAPE>
geoElement<GEOSHAPE>::
geoElement(const UInt & N) : pGraphItem(N)
{
  assert(numPoints == N);
  fixedLength = true;
}

template<typename GEOSHAPE>
geoElement<GEOSHAPE>::
geoElement(const geoElement<GEOSHAPE> & E) : pGraphItem(E)
{
  geoId       = E.geoId;
  fixedLength = E.fixedLength;
}

template<typename GEOSHAPE>
geoElement<GEOSHAPE>
geoElement<GEOSHAPE>::
operator=(const geoElement<GEOSHAPE> & E)
{
  geoId       = E.geoId;
  connected   = E.connected;
  fixedLength = E.fixedLength;
  
  updateSorting();
  
  return *this;
}

template<typename GEOSHAPE>
const UInt &
geoElement<GEOSHAPE>::
getGeoId() const
{
  return(geoId);
}

template<typename GEOSHAPE>
UInt &
geoElement<GEOSHAPE>::
getGeoId()
{
  return(geoId);
}

template<typename GEOSHAPE>
void
geoElement<GEOSHAPE>::
setGeoId(const UInt & Id)
{
  geoId = Id;
}

template<typename GEOSHAPE>
void
geoElement<GEOSHAPE>::
resize(const UInt & dim)
{
  assert(!fixedLength);
  pGraphItem::resize(dim);
}

template<typename GEOSHAPE>
void
geoElement<GEOSHAPE>::
reserve(const UInt & n)
{
  assert(!fixedLength);
  pGraphItem::reserve(n);
}

template<typename GEOSHAPE>
void
geoElement<GEOSHAPE>::
push_back(const UInt & value, const bool & update)
{
  assert(!fixedLength);
  pGraphItem::push_back(value,update);
}

template<typename GEOSHAPE>
void
geoElement<GEOSHAPE>::
merge(const pGraphItem & G)
{
  assert(!fixedLength);
  pGraphItem::merge(G);
}

template<typename GEOSHAPE>
template<class ARK>
void
geoElement<GEOSHAPE>::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  //Number of items
  int  N = connected.size();
  ar & N;
  
  //Items
  if(UInt(N) != connected.size())
  { connected.resize(N); }
  
  for(int i=1; i <= N; ++i)
  { ar & connected(i); }
  
  //Id
  ar & geoId;
  
  //Updating
  updateSorting();
}

template<typename GEOSHAPE>
ostream & operator<<(ostream & f, const geoElement<GEOSHAPE> & G)
{
  f << "GeoId          : " << G.geoId << endl;
  f << "Num Connected  : " << G.size() << endl;
  f << "Connected Id's : ";
  
  for(UInt i=1; i <= G.size(); ++i)
  {
    f << G(i) << " ";
  }
  f << endl;
  
  return(f);
}


#endif
