/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef STATICVECTOR_HPP
#define STATICVECTOR_HPP

#include <assert.h>
#include <iostream>
#include <fstream>

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>
#include <boost/array.hpp>

#include "morganaTypes.hpp"



using namespace std;


//_________________________________________________________________________________________________
// EXPRESSIONS
//-------------------------------------------------------------------------------------------------
template<typename LHS, typename RHS>
struct plusStaticVector                                     //The sum operation
{
  static Real apply(const UInt & i, const LHS & left, const RHS & right) 
  {
    assert(left.size() == right.size());
    return(left(i) + right(i));
  }
};

template<typename LHS, typename RHS>
struct minusStaticVector                                    //The minus operation
{
  static Real apply(const UInt & i, const LHS & left, const RHS & right)
  {
    assert(left.size() == right.size());
    return(left(i) - right(i));
  }
};

template<typename LHS>
struct multStaticVector                                     //The sum operation
{
  static Real apply(const UInt & i, const LHS & left, const Real & right)
  { return(left(i) * right); }
};

template<typename LHS>
struct divStaticVector                                      //The div operation
{
  static Real apply(const UInt & i, const LHS & left, const Real & right)
  { return(left(i) / right); }
};

template<typename LHS, typename RHS, typename OP>          //The expression template
struct expressionStaticVector
{
    static const morganaTypes myType = typeStaticVector;
  
    const LHS & left;
    const RHS & right;
  
    expressionStaticVector(const LHS & Left, const RHS & Right) : left(Left), right(Right) { }
    Real operator()(const UInt & i) const { return( OP::apply(i,left,right) ); }
    UInt  size() const {return(left.size());}
};



//_________________________________________________________________________________________________
// STATIC VECTOR CLASS
//-------------------------------------------------------------------------------------------------

/*! Implements an algebraic vector with variable leghth. The vector legth is defined statically */
template<size_t N>
class staticVector : public boost::array<Real,N>
{
    /*! Typedefs */  
  public:
    typedef boost::array<Real,N> BASE;
  
    
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    UInt srlzSize;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    
     /*! @name Internal data */ //@{
  public:
    static const morganaTypes myType = typeStaticVector;
    //@}
    
    
    /*! @name Internal constructors */ //@{
  public:
    /*! Standard constructor */
    staticVector();
    
    /*! Copy constructor */
    staticVector(const staticVector & V);
    
    /*! Expression template constructor */
    template<typename LHS, typename RHS, typename OP>
    staticVector(const expressionStaticVector<LHS,RHS,OP> & Expression);
    //@}
    
    
    /*! @name Operators */ //@{
  public:
    /*! Access operator */
    inline Real & operator()(const UInt & i);
    
    /*! Access operator */
    inline const Real & operator()(const UInt & i) const;
    
    /*! The equality operator */
    staticVector & operator=(const staticVector & C);
    
    /*! The equality operator for expressions */
    template<typename LHS, typename RHS, typename OP>
    void operator=(const expressionStaticVector<LHS,RHS,OP> & Expression);
    
    /*! Sum operator */
    void operator+=(const staticVector & C);
    
    /*! Subtract operator */
    void operator-=(const staticVector & C);

    /*! Division operator */
    void operator/=(const Real & a);
    
    /*! Product operator */
    void operator*=(const Real & a);
    
    /*! Scalar product */
    Real operator*(const staticVector & C) const;
    //@}
    
    
    /*! @name Internal functions */ //@{
  public:
    /*! Size of the vector */
    static UInt size();
    
    /*! Components sum */
    Real sum() const;
    
    /*! One-norm */
    Real norm1() const;
    
    /*! One-norm */
    static Real norm1(const staticVector<N> & V);
    
    /*! Inf-norm */
    static Real normInf(const staticVector<N> & V);
    
    /*! Set components to zero */
    void cancel();
    //@}
    
    /*! @name Ordinal functions */ //@{
  public:
    /*! Less operator: a vector is "less" than another if its first coordinate is less than the other vctor one.
    If the first component value is equal the second one is considered and so on */
    bool operator<(const staticVector & V) const;
    
    /*! Inequality operator: two vectors are equal if their components are equal to the geometric tolerance */
    bool operator!=(const staticVector & V) const;
    //@}
    
    
    /*! @name Other functions */ //@{
  public:
    /*! Outstream operator */
    template<size_t Nc>
    friend ostream & operator<<( ostream & f, const staticVector<Nc> & A);
    
    /*! Memory size */
    size_t memSize() const;
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<size_t N>
staticVector<N>::
staticVector()
{
  for(UInt i=0; i < N; ++i)
  { BASE::operator[](i) = 0.0; }
}

template<size_t N>
staticVector<N>::
staticVector(const staticVector<N> & V)
{
  for(UInt i=0; i < N; ++i)
  { BASE::operator[](i) = V[i]; }
}

template<size_t N>
staticVector<N> &
staticVector<N>::
operator=(const staticVector & C)
{
  for(UInt i=0; i < N; ++i)
  { BASE::operator[](i) = C[i]; }
  
  return(*this);
}

template<size_t N>
template<typename LHS, typename RHS, typename OP>
staticVector<N>::
staticVector(const expressionStaticVector<LHS,RHS,OP> & Expression)
{
  assert(Expression.size() == this->size());
  
  for(UInt i=1; i <= this->size(); ++i)
  { this->operator()(i) = Expression(i); }
}

template<size_t N>
template<class ARK>
void
staticVector<N>::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  srlzSize = N;
  ar & srlzSize;
  
  assert(srlzSize == N);
  
  for(UInt i=1; i <= srlzSize; ++i)
  { ar & this->operator()(i); }
}
    


//_________________________________________________________________________________________________
// OPERATORS
//-------------------------------------------------------------------------------------------------
template<size_t N>
Real & 
staticVector<N>::
operator()(const UInt & i)
{
  assert(i > 0);
  assert(i <= BASE::size());
  return(BASE::operator[](i-1));
}

template<size_t N>
const Real & 
staticVector<N>::
operator()(const UInt & i) const
{
  assert(i > 0);
  assert(i <= BASE::size());
  return(BASE::operator[](i-1));
}

template<size_t N>
template<typename LHS, typename RHS, typename OP>
void
staticVector<N>::
operator=(const expressionStaticVector<LHS,RHS,OP> & Expression)
{
  assert(Expression.size() == this->size());
  
  for(UInt i=1; i <= this->size(); ++i)
  { this->operator()(i) = Expression(i); }
}

template<size_t N>
void 
staticVector<N>::
operator+=(const staticVector<N> &C)
{  
  for(UInt i=1; i <= this->size(); ++i)
  { this->operator()(i) += C(i); }
}

template<size_t N>
void
staticVector<N>::
operator-=(const staticVector<N> & C)
{
  for(UInt i=1; i <= this->size(); ++i)
  { this->operator()(i) -= C(i); }
}

template<size_t N>
void
staticVector<N>::
operator*=(const Real &a)
{
  for(UInt i=1; i <= this->size(); ++i)
  { this->operator()(i) *= a; }
}

template<size_t N>
void 
staticVector<N>::
operator/=(const Real &a)
{
  for(UInt i=1; i <= this->size(); ++i)
  { this->operator()(i) /= a; }
}

template<size_t N>
Real
staticVector<N>::
operator*(const staticVector<N> &C) const
{
  Real tot=0.0;
  
  for(UInt i=1; i<=this->size(); ++i)
  { tot += this->operator()(i) * C(i); }
  
  return(tot);
}



//_________________________________________________________________________________________________
// INTERNAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<size_t N>
UInt 
staticVector<N>::
size()
{
  return(N);
}

template<size_t N>
Real 
staticVector<N>::
sum() const
{
  Real tot = 0.0;
  
  for(UInt i=1; i<=this->size(); ++i)
  { tot += this->operator()(i); }
  
  return(tot);
}

template<size_t N>
Real
staticVector<N>::
norm1() const
{
  Real tot = 0.0;
  
  for(UInt i=1; i<=this->size(); ++i)
  { tot += abs(this->operator()(i)); }
  
  return(tot);
}

template<size_t N>
Real
staticVector<N>::
norm1(const staticVector<N> & V)
{
  Real tot = 0.0;
  
  for(UInt i=1; i <= V.size(); ++i)
  { tot += abs(V(i)); }
  
  return(tot);
}

template<size_t N>
Real
staticVector<N>::
normInf(const staticVector<N> & V)
{
  Real tot = 0.0;
  
  for(UInt i=1; i <= V.size(); ++i)
  { tot = std::max(tot, abs(V(i)) ); }
  
  return(tot);
}

template<size_t N>
void
staticVector<N>::
cancel()
{
  for(UInt i=1; i <= this->size(); ++i)
  { this->operator()(i) = 0.0; }
}


//_________________________________________________________________________________________________
// ORDINAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<size_t N>
bool
staticVector<N>::
operator<(const staticVector & V) const
{
  for(int i=1; i <= N; ++i)
  {
    if( this->operator()(i) < (V.operator()(i) - geoToll) )
    {
      return(true);
    }
      
    if( this->operator()(i) > (V.operator()(i) - geoToll) )
    {
      return(false);
    }
  }
  
  return(false);
}

template<size_t N>
bool 
staticVector<N>::
operator!=(const staticVector & V) const
{  
  for(int i=1; i <= N; ++i)
  {
    if( ( this->operator()(i) > (V.operator()(i) + geoToll) ) || ( this->operator()(i) < (V.operator()(i) - geoToll) ) )
    {
      return(true);
    }
  }
  
  return(false);
}



//_________________________________________________________________________________________________
// OTHER FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<size_t Nc>
ostream & 
operator<<( ostream & f, const staticVector<Nc> & A)
{
  for(UInt i=1; i <= A.size(); ++i)
  { f << A(i) << " "; }
  
  f << endl;
  
  return(f);
}

template<size_t N>
size_t
staticVector<N>::
memSize() const
{
  return(N * sizeof(double) + sizeof(UInt));
}

#endif

