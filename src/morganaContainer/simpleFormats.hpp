/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef SIMPLEFORMATS_HPP
#define SIMPLEFORMATS_HPP

#include <vector>
#include <assert.h>

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "morganaTypes.hpp"

using namespace std;


//_________________________________________________________________________________________________
//                                       serialVector CLASS
//_________________________________________________________________________________________________

/*! Simple vect Class. In order to include the parallel support the \c sVect is serializable.
Some traits specilizations are included for the \c Real and \c UInt types. In general all the 
not-user-defined types should have a proper serialization specialization. */
template <typename T, int OFFSETVEC = 1>
class sVect : public vector<T>
{ 
    /*! @name Typedefs */ //@{
  public:
    typedef std::vector<T> raw_container;
    typedef typename raw_container::size_type size_type;
    typedef typename raw_container::reference reference;
    typedef typename raw_container::const_reference const_reference;
    //@}
    

    /*! @name Parallel support */ //@{
  public:    
    friend class boost::serialization::access;
    
    UInt srlzSize; //Serialization size
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    
    /*! @name Constructors */ //@{
  public:
    /*! Constructor */
    explicit sVect( size_type i ) : raw_container( i ) {};
    
    /*! Constructor */
    explicit sVect( const raw_container & );
    
    /*! Constructor */
    sVect() : raw_container() {};
    
    /*! Copy constructor */
    sVect( const sVect<T, OFFSETVEC> & );
    
    /*! Destructor */
   ~sVect() {}
   //@}


    /*! @name Operators */ //@{
  public:
    /*! Equality operator */
    sVect<T, OFFSETVEC> & operator=( const sVect<T, OFFSETVEC> & );
    
    /*! Get - with memory check */
    T & get( size_type i );
    
    /*! Get - with memory check */
    const T & get( size_type i ) const;
    
    /*! Get - NOT memory check */
    inline reference operator() ( size_type const i );
    
    /*! Get - NOT memory check */
    inline const_reference operator() ( size_type const i ) const;
    
    /*! Data clean */
    void clean();

    /*! Memory check */
    inline bool bCheck( size_type const i ) const
    { return i >= OFFSETVEC && i < this->size() + OFFSETVEC ; }
    //@}
    
    /*! @name Printout */ //@{
  public:
    template<typename TT, int OO>
    friend ostream & operator<<(ostream & f, const sVect<TT,OO> & V);
    //@}
};


template <typename T, int OFFSETVEC>
sVect<T, OFFSETVEC>::
sVect( const sVect<T, OFFSETVEC> & v ) : raw_container( v )
{}

template <typename T, int OFFSETVEC>
sVect<T, OFFSETVEC> &
sVect<T, OFFSETVEC>::
operator=( const sVect<T, OFFSETVEC> & v )
{
    raw_container::operator=( v );
    return *this;
}

template <typename T, int OFFSETVEC>
typename sVect<T, OFFSETVEC>::reference
sVect<T, OFFSETVEC>::
operator() ( size_type const i )
{ 
  assert(bCheck( i ));
  return ( this->operator[] ( i - OFFSETVEC ) );
}

template <typename T, int OFFSETVEC>
typename sVect<T, OFFSETVEC>::const_reference 
sVect<T, OFFSETVEC>::
operator() ( size_type const i ) const
{ 
  assert(bCheck( i ));
  return ( this->operator[] ( i - OFFSETVEC ) ); 
}

template <typename T, int OFFSETVEC>
T & sVect<T, OFFSETVEC>::
get( size_type i )
{
  assert(bCheck( i ) );
  return *( this->begin() + ( i - OFFSETVEC ) );
}

template <typename T, int OFFSETVEC>
const T& sVect<T, OFFSETVEC>::
get( size_type i ) const
{
  assert(bCheck( i ) );
  return *( this->begin() + ( i - OFFSETVEC ) );
}

template <typename T, int OFFSETVEC>
void 
sVect<T, OFFSETVEC>::
clean()
{
  raw_container tmp;
  this->clear();
  this->swap( tmp );
}

template <typename T, int OFFSETVEC>
template<class ARK>
void
sVect<T, OFFSETVEC>::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  srlzSize = this->size();
  ar & srlzSize;
  
  if(srlzSize != this->size())
  { this->resize(srlzSize); }
  
  for(UInt i=1; i <= this->size(); ++i)
  { this->get(i).serialize(ar,version); }
}



/*! Specialization for Real */
template <>
template<class ARK>
void
sVect<Real,1>::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  srlzSize = this->size();
  ar & srlzSize;
  
  if(srlzSize != this->size())
  { this->resize(srlzSize); }
  
  for(UInt i=1; i <= this->size(); ++i)
  { ar & this->get(i); }
}

/*! Specialization for UInt */
template <>
template<class ARK>
void
sVect<UInt,1>::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  srlzSize = this->size();
  ar & srlzSize;
  
  if(srlzSize != this->size())
  { this->resize(srlzSize); }
  
  for(UInt i=1; i <= this->size(); ++i)
  { ar & this->get(i); }
}




template<typename TT, int OO>
ostream & operator<<(ostream & f, const sVect<TT,OO> & V)
{
  for(UInt i=0; i < V.size(); ++i)
  {
    f << V[i] << endl;
  }
  
  return(f);
}




//_________________________________________________________________________________________________
//                                       serialArray CLASS
//_________________________________________________________________________________________________
/*! Classe simpleArray 
<b> Esempio </b>
sArray<int> b(3,5) // an arrray with 3 rows and 5 columns
b(3,2)=5;
b.reshape(2,3) // now b is 2x3
*/
template <typename T, int OFFSETVEC = 1>
class sArray : public vector<T>
{
    /*! @name Typedefs */ //@{
  public:
    typedef vector<T> raw_container;
    typedef typename raw_container::size_type size_type;
    typedef typename raw_container::reference reference;
    typedef typename raw_container::const_reference const_reference;
    typedef typename raw_container::iterator iterator;
    typedef typename raw_container::const_iterator const_iterator;
    //@}
    
    /*! @name Internal data */ //@{
  private:
    size_type _M_nrows;
    size_type _M_ncols;
    //@}
    
    /*! @name Constructors */ //@{ 
  public:
    explicit sArray( size_type ntot );
    explicit sArray();
    explicit sArray( size_type nrows, size_type ncols );
    ~sArray() {}
    //@}

    /*! @name Operators */ //@{ 
  public:
    reference operator() ( size_type const i );
    const_reference operator() ( size_type const i ) const;
    reference operator() ( size_type const i, size_type const j );
    const_reference operator() ( size_type const i, size_type const j ) const;
    iterator columnIterator( size_type const col );
    //@}
    
    /*! @name Size functions */ //@{ 
  public:
    void reshape( sArray<T, OFFSETVEC>::size_type const n, size_type const m );
    void clean();
    bool bCheck( size_type const i, size_type const j ) const;
    inline size_type nrows() const;
    inline size_type ncols() const;
    //@}
};


template <typename T, int OFFSETVEC>
sArray<T, OFFSETVEC>::
sArray() : vector<T>(), _M_nrows( 0 ), _M_ncols( 1 )
{}

template <typename T, int OFFSETVEC>
sArray<T, OFFSETVEC>::
sArray( size_type ntot ) : vector<T>( ntot ), _M_nrows( ntot ), _M_ncols( 1 )
{}

template <typename T, int OFFSETVEC>
sArray<T, OFFSETVEC>::
sArray( size_type nrows, size_type ncols ) : vector<T>( nrows*ncols ), _M_nrows( nrows ), _M_ncols( ncols )
{}

template <typename T, int OFFSETVEC>
void
sArray<T, OFFSETVEC>::
reshape( size_type nrows, size_type ncols )
{
  raw_container::resize( nrows * ncols );
  _M_nrows = nrows;
  _M_ncols = ncols;
}

template <typename T, int OFFSETVEC>
typename sArray<T, OFFSETVEC>::reference
sArray<T, OFFSETVEC>::
operator() ( size_type const i )
{
  return * ( this->begin() + ( i - OFFSETVEC ) );
}

template <typename T, int OFFSETVEC>
typename sArray<T, OFFSETVEC>::const_reference
sArray<T, OFFSETVEC>::
operator() ( size_type const i ) const
{
  return * ( this->begin() + ( i - OFFSETVEC ) );
}

template <typename T, int OFFSETVEC>
typename sArray<T, OFFSETVEC>::reference
sArray<T, OFFSETVEC>::
operator() ( size_type const i, size_type const j )
{
  return * ( this->begin() + ( j - OFFSETVEC ) * _M_nrows + ( i - OFFSETVEC ) );
}

template <typename T, int OFFSETVEC>
typename sArray<T, OFFSETVEC>::const_reference 
sArray<T, OFFSETVEC>::
operator() ( size_type const i, size_type const j ) const
{
  return * ( this->begin() + ( j - OFFSETVEC ) * _M_nrows + ( i - OFFSETVEC ) );
}

template <typename T, int OFFSETVEC>
typename sArray<T, OFFSETVEC>::iterator 
sArray<T, OFFSETVEC>::
columnIterator( size_type const col )
{
  if ( col > _M_ncols ) {return typename sArray<T, OFFSETVEC>::iterator();}
  else                  {return this->begin() + ( col - OFFSETVEC ) * _M_nrows;}
}

template <typename T, int OFFSETVEC>
void 
sArray<T, OFFSETVEC>::
clean()
{
  raw_container tmp;
  this->clear();
  this->swap( tmp );
  _M_nrows = 0;
  _M_ncols = 0;
}

template <typename T, int OFFSETVEC>
bool 
sArray<T, OFFSETVEC>::
bCheck( size_type const i, size_type const j ) const
{
  return i >= OFFSETVEC && i - OFFSETVEC + ( j - OFFSETVEC ) * _M_nrows < this->size() ;
}

template <typename T, int OFFSETVEC>
typename sArray<T, OFFSETVEC>::size_type
sArray<T, OFFSETVEC>::
nrows() const
{
  return _M_nrows;
}

template <typename T, int OFFSETVEC>
typename sArray<T, OFFSETVEC>::size_type 
sArray<T, OFFSETVEC>::
ncols() const
{
  return _M_ncols;
}

#endif
