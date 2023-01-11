/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef KOMPLEX_H
#define KOMPLEX_H

#include <assert.h>
#include <iostream>
#include <fstream>

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "morganaTypes.hpp"

using namespace std;


//_________________________________________________________________________________________________
// EXPRESSIONS
//-------------------------------------------------------------------------------------------------
template<typename LHS, typename RHS>                       //The sum operation
struct plusKomplex
{
  static Real applyReal(const LHS & left, const RHS & right) { return(left.getReal() + right.getReal()); }
  static Real applyImag(const LHS & left, const RHS & right) { return(left.getImag() + right.getImag()); }
};

template<typename LHS, typename RHS>                       //The sum operation
struct minusKomplex
{
  static Real applyReal(const LHS & left, const RHS & right) { return(left.getReal() - right.getReal()); }
  static Real applyImag(const LHS & left, const RHS & right) { return(left.getImag() - right.getImag()); }
};

template<typename LHS>                                     //The multiplication operation
struct multKomplex
{
  static Real applyReal(const LHS & left, const Real & right) { return(left.getReal() * right); }
  static Real applyImag(const LHS & left, const Real & right) { return(left.getImag() * right); }
};

template<typename LHS>                                     //The multiplication operation
struct divKomplex
{
  static Real applyReal(const LHS & left, const Real & right) { return(left.getReal() / right); }
  static Real applyImag(const LHS & left, const Real & right) { return(left.getImag() / right); }
};

template<typename LHS, typename RHS>                       //Complex-complex multiplication
struct multComplexComplex
{
  static Real applyReal(const LHS & left, const RHS & right) { return(left.getReal() * right.getReal() - left.getImag() * right.getImag()); }
  static Real applyImag(const LHS & left, const RHS & right) { return(left.getReal() * right.getImag() + left.getImag() * right.getReal()); }
};

template<typename LHS, typename RHS>                       //Complex-complex division
struct divComplexComplex
{
  static Real applyReal(const LHS & left, const RHS & right)
  { return( (left.getReal() * right.getReal() + left.getImag() * right.getImag()) / (right.getReal() * right.getReal() + right.getImag() * right.getImag())); }
  
  static Real applyImag(const LHS & left, const RHS & right)
  { return( (left.getImag() * right.getReal() - left.getReal() * right.getImag()) / (right.getReal() * right.getReal() + right.getImag() * right.getImag())); }
};

template<typename LHS, typename RHS, typename OP>          //The expression template
struct expressionKomplex
{
    static const morganaTypes myType = typeKomplex;
  
    const LHS & left;
    const RHS & right;
  
    expressionKomplex(const LHS & Left, const RHS & Right) : left(Left), right(Right) { }
    Real getReal() const { return(OP::applyReal(left,right)); }
    Real getImag() const { return(OP::applyImag(left,right)); }
};


//_________________________________________________________________________________________________
// POINT3D CLASS
//-------------------------------------------------------------------------------------------------
/*! The complex scalar number */
class komplex
{
  /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
  
  
    /*! @name Internal data */ //@{
  public:
    static const morganaTypes myType = typeKomplex;
    
    /// The components
    Real re, im;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    /*! Constructor */
    komplex(Real Re = 0.0, Real Im = 0.0);
    
    /*! Constructor */
    template<typename LHS, typename RHS, typename OP>
    komplex(const expressionKomplex<LHS,RHS,OP> & Expression);
    
    /*! Copy constructor */
    komplex(const komplex & V);

    /*! Destructor */
    virtual ~komplex();
    //@}
    
    /*! @name Set functions */ //@{
  public:
    /*! Set all components */
    void set(const Real & Re, const Real & Im);
    
    /*! Set the Real component */
    void setReal(const Real & Re);
    
    /*! Set the Imag component */
    void setImag(const Real & Im);

    /*! Set the i-th component */
    void setI(const UInt & i, const Real & val);
    //@}
    
    /*! @name Get functions */ //@{
  public:
    /*! Get the Real component */
    inline Real getReal();
    inline const Real & getReal() const;

    /*! Get the Imag component */
    inline Real getImag();
    inline const Real & getImag() const;

    /*! Get the i-th component */
    inline Real getI(UInt i);
    inline const Real & getI(const UInt & i) const;
    //@}
    
    /*! @name Operators */ //@{
  public:
    /*! Access operator */    
    Real operator()(const UInt & i);
    
    /*! Access operator */
    const Real & operator()(const UInt & i) const;
    
    /*! The equality operator */
    komplex & operator=(const komplex & C);
    
    /*! The equality operator for expressions */
    template<typename LHS, typename RHS, typename OP>
    void operator=(const expressionKomplex<LHS,RHS,OP> & Expression);
    
    /*! Sum operator */
    void operator+=(const komplex & C);
    
    /*! Subtract operator */
    void operator-=(const komplex & C);
    
    /*! Product operator */
    void operator*=(const Real & a);
    
    /*! Division operator */
    void operator/=(const Real & a);
    
    /*! Product operator */
    void operator*=(const komplex & C);
    
    /*! Division operator */
    void operator/=(const komplex & C);
    //@}
    
    /*! @name Ordinal functions */ //@{
  public:
    /*! Less operator: a vector is "less" than another if its first coordinate (Real then Imag) is less than the other vctor one.
    If the first component value is equal the second one is considered and so on */
    bool operator<(const komplex & C) const;
    
    /*! Inequality operator: two vectors are equal if their components (Real then Imag) are equal to the geometric tolerance */
    bool operator!=(const komplex & C) const;
    
    /*! Equality operator */
    bool operator==(const komplex & C) const;
    //@}
    
    /*! @name Other functions */ //@{
  public:
    /*! The norm of the R3 vector */
    Real norm() const;
    
    /*! The phase */
    Real phase() const;
    
    /*! Static norm function */
    static Real norm(const komplex & C);
    
    /*! The phase */
    static Real phase(const komplex & C);
    
    /*! Static Real function */
    static Real cReal(const komplex & C);
    static Real cReal(const Real & C);
    
    /*! Static Imag function */
    static Real cImag(const komplex & C);
    static Real cImag(const Real & C);
    
    /*! Static Access functions */
    static Real cComp(const komplex & C, const UInt & i);
    static Real cComp(const Real & C, const UInt & i);
    
    /*! Complex exponential */
    static komplex iexp(const Real & theta);
    
    /*! Complex conjugate */
    static komplex conj(const komplex & C);
    
    /*! Outstream operator */
    friend ostream & operator<<(ostream & f, const komplex & C);
    
    /*! Memory size */
    size_t memSize() const;
    //@}
};


//_________________________________________________________________________________________________
// INLINED AND TEMPLATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
Real
komplex::
getReal()
{
  return(re);
}

const Real &
komplex::
getReal() const
{
  return(re);
}

Real
komplex::
getImag()
{
  return(im);
}

const Real &
komplex::
getImag() const
{
  return(im);
}

Real
komplex::
getI(UInt i)
{
  assert(i <= 2);
  assert(i >= 1);
  
  if(i == 1) { return(re); }
  else       { return(im); }
}

const Real &
komplex::
getI(const UInt & i) const
{
  assert(i <= 2);
  assert(i >= 1);
  
  if(i == 1) { return(re); }
  else       { return(im); }
}

template<typename LHS, typename RHS, typename OP>
komplex::
komplex(const expressionKomplex<LHS,RHS,OP> & Expression)
{
  re = Expression.getReal();
  im = Expression.getImag();
}

template<typename LHS, typename RHS, typename OP>
void
komplex::
operator=(const expressionKomplex<LHS,RHS,OP> & Expression)
{
  re = Expression.getReal();
  im = Expression.getImag();
}

template<class ARK>
void
komplex::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & re;
  ar & im;
}


//_________________________________________________________________________________________________
// OPERATORS
//-------------------------------------------------------------------------------------------------
inline expressionKomplex<komplex,komplex,multComplexComplex<komplex,komplex> >
operator*(const komplex & left, const komplex & right)
{
  return(expressionKomplex<komplex,komplex,multComplexComplex<komplex,komplex> >(left,right));
}

inline expressionKomplex<komplex,komplex,divComplexComplex<komplex,komplex> >
operator/(const komplex & left, const komplex & right)
{
  return(expressionKomplex<komplex,komplex,divComplexComplex<komplex,komplex> >(left,right));
}

#endif
