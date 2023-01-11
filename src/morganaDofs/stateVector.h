/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef STATEVECTOR_H
#define STATEVECTOR_H

#include <assert.h>
#include <iostream>
#include <fstream>

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "morganaTypes.hpp"

#include "Epetra_SerialDenseVector.h"

using namespace std;


//_________________________________________________________________________________________________
// EXPRESSIONS
//-------------------------------------------------------------------------------------------------
template<typename LHS, typename RHS>
struct plusStateVector                                     //The sum operation
{
  static Real apply(const UInt & i, const LHS & left, const RHS & right) 
  {
    assert(left.Length() == right.Length());
    return(left(i) + right(i));
  }
};

template<typename LHS, typename RHS>
struct minusStateVector                                    //The minus operation
{
  static Real apply(const UInt & i, const LHS & left, const RHS & right)
  {
    assert(left.Length() == right.Length());
    return(left(i) - right(i));
  }
};

template<typename LHS>
struct multStateVector                                     //The sum operation
{
  static Real apply(const UInt & i, const LHS & left, const Real & right)
  { return(left(i) * right); }
};

template<typename LHS>
struct divStateVector                                      //The div operation
{
  static Real apply(const UInt & i, const LHS & left, const Real & right)
  { return(left(i) / right); }
};

template<typename LHS, typename RHS, typename OP>          //The expression template
struct expressionStateVector
{
    static const morganaTypes myType = typeStateVector;
  
    const LHS & left;
    const RHS & right;
  
    expressionStateVector(const LHS & Left, const RHS & Right) : left(Left), right(Right) { }
    Real operator()(const UInt & i) const { return( OP::apply(i,left,right) ); }
    int  Length() const {return(left.Length());}
};



//_________________________________________________________________________________________________
// STATE VECTOR CLASS
//-------------------------------------------------------------------------------------------------
/*! Dynamically sized varible lengh vector. */
class stateVector : public Epetra_SerialDenseVector
{
  /*! @name Parallel support */ //@{
  private:
    friend class boost::serialization::access;
    
    UInt srlzSize; //Serialization size
    
    /*! Serialization function */
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    
    /*! @name Internal data */ //@{
  public:
    static const morganaTypes myType = typeStateVector;
    //@}
    
    
    /*! @name Constructors */ //@{
  public:
    /*! Constructor
    \param N number of components */
    stateVector(const UInt & N = 0);
    
    /*! Constructor
    \param Comp vector
    \param N number of components*/
    stateVector(Real * Comp, const UInt & N);
    
    /*! Copy constructor
    \param C vector */
    stateVector(const stateVector & C);
    
    /*! Costructor for expression template  */
    template<typename LHS, typename RHS, typename OP>
    stateVector(const expressionStateVector<LHS,RHS,OP> & Expression);
    //@}
    
    /*! @name Operators */ //@{
  public:
    /*! Access operator */
    inline Real & operator()(const UInt & i);
    
    /*! Access operator */
    inline const Real & operator()(const UInt & i) const;
    
    /*! The equality operator for expressions */
    template<typename LHS, typename RHS, typename OP>
    void operator=(const expressionStateVector<LHS,RHS,OP> & Expression);
    
    /*! Add operator */
    void operator+=(const stateVector & C);
    
    /*! Subtract operator */
    void operator-=(const stateVector & C);

    /*! Divide operator */
    void operator/=(const Real & a);
    
    /*! Multiplication operator */
    void operator*=(const Real & a);
    
    /*! Scalar product */
    Real operator*(const stateVector & C) const;
    //@}
    
    /*! @name Ordinal functions */ //@{
  public:
    /*! Less operator: a vector is "less" than another if its first coordinate is less than the other vctor one.
    If the first component value is equal the second one is considered and so on */
    bool operator<(const stateVector & V) const;
    
    /*! Inequality operator: two vectors are equal if their components are equal to the geometric tolerance */
    bool operator!=(const stateVector & V) const;
    //@}
    
    /*! @name Internal functions */ //@{
  public:
    /*! Sum of the components */
    Real sum() const;
    
    /*! One-Norm */
    Real norm1() const;
    
    /*! Set components to zero */
    void cancel();
    
    /*! Get the size */
    UInt size() const;
    //@}
    
    /*! @name Other functions */ //@{
  public:
    /*! Outstream operator */
    friend ostream & operator<<( ostream & f, const stateVector & A);
    
    /*! Memory size */
    size_t memSize() const;
    //@}
};



//_________________________________________________________________________________________________
// TEMPLATED AND INLINED FUNCTIONS
//-------------------------------------------------------------------------------------------------
Real & 
stateVector::
operator()(const UInt & i)
{
  return(Epetra_SerialDenseVector::operator()(i-1));
}

const Real & 
stateVector::
operator()(const UInt & i) const
{
  return(Epetra_SerialDenseVector::operator()(i-1));
}

template<typename LHS, typename RHS, typename OP>
stateVector::
stateVector(const expressionStateVector<LHS,RHS,OP> & Expression)
{
  if(Expression.Length() != this->Length())
  { this->Resize(Expression.Length()); }
  
  for(int i=1; i <= this->Length(); ++i)
  { this->operator()(i) = Expression(i); }
}

template<typename LHS, typename RHS, typename OP>
void
stateVector::
operator=(const expressionStateVector<LHS,RHS,OP> & Expression)
{
  if(Expression.Length() != this->Length())
  { this->Resize(Expression.Length()); }
  
  for(int i=1; i <= this->Length(); ++i)
  { this->operator()(i) = Expression(i); }
}

template<class ARK>
void
stateVector::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  srlzSize = this->Length();
  
  ar & srlzSize;
  
  if(srlzSize != this->Length())
  { this->Resize(srlzSize); }
  
  for(int i=1; i <= srlzSize; ++i)
  { ar & this->operator()(i); }
}

#endif
