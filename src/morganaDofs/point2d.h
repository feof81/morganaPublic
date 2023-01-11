/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef POINT2D_H
#define POINT2D_H

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
struct plusPoint2d
{
  static Real applyX(const LHS & left, const RHS & right) { return(left.getX() + right.getX()); }
  static Real applyY(const LHS & left, const RHS & right) { return(left.getY() + right.getY()); }
};

template<typename LHS, typename RHS>                       //The difference operation
struct minusPoint2d
{
  static Real applyX(const LHS & left, const RHS & right) { return(left.getX() - right.getX()); }
  static Real applyY(const LHS & left, const RHS & right) { return(left.getY() - right.getY()); }
};

template<typename LHS>                                     //The multiplication operation
struct multPoint2d
{
  static Real applyX(const LHS & left, const Real & right) { return(left.getX() * right); }
  static Real applyY(const LHS & left, const Real & right) { return(left.getY() * right); }
};

template<typename LHS>                                     //The multiplication operation
struct divPoint2d
{
  static Real applyX(const LHS & left, const Real & right) { return(left.getX() / right); }
  static Real applyY(const LHS & left, const Real & right) { return(left.getY() / right); }
};

template<typename LHS, typename RHS, typename OP>          //The expression template
struct expressionPoint2d
{
  static const morganaTypes myType = typePoint2d;
  
  const LHS & left;
  const RHS & right;
  
  expressionPoint2d(const LHS & Left, const RHS & Right) : left(Left), right(Right) { }
  Real getX() const { return(OP::applyX(left,right)); }
  Real getY() const { return(OP::applyY(left,right)); }
};



//_________________________________________________________________________________________________
// POINT3D CLASS
//-------------------------------------------------------------------------------------------------

/*! Two dimensional point */
class point2d
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
  
  
    /*! @name INternal data */ //@{
  public:
    static const morganaTypes myType = typePoint2d;
    
    /// The coordinates
    Real X[2];
    
    /// The node id
    UInt id;
    //@}
    
  
    /*! @name Costruttori */ //@{
  public:
    /*! Constructor */
    point2d(Real xx=0.0, Real yy=0.0);
    
    /*! Constructor */
    template<typename LHS, typename RHS, typename OP>
    point2d(const expressionPoint2d<LHS,RHS,OP> & Expression);
    
    /*! Copy constructor */
    point2d(const point2d & V);

    /*! Destructor */
    virtual ~point2d();
    //@}
    
    
    /*! @name Set functions */ //@{
  public:
    /*! Set all components of the point */
    void set(const Real & xx, const Real & yy);
    
    /*! Set the X component */
    void setX(const Real & xx);
    
    /*! Set the Y component */
    void setY(const Real & yy);

    /*! Set the i-th component */
    void setI(const UInt & i, const Real & val);
    
    /*! Set the Id*/
    void setId(const UInt & Id);
    //@}
    
    
    /*! @name Get functions */ //@{
  public:
    /*! Get the X component */
    inline Real getX();
    inline const Real & getX() const;

    /*! Get the Y component */
    inline Real getY();
    inline const Real & getY() const;

    /*! Get the i-th component */
    inline Real getI(UInt i);
    inline const Real & getI(const UInt & i) const;
      
    /*! Get the Id */
    UInt getId() const;
    //@}
    
    
    /*! @name Operators */ //@{
  public:
    /*! Access operator */    
    Real operator()(const UInt & i);
    
    /*! Access operator */
    const Real & operator()(const UInt & i) const;
    
    /*! The equality operator */
    point2d & operator=(const point2d & V);
    
    /*! The equality operator for expressions */
    template<typename LHS, typename RHS, typename OP>
    void operator=(const expressionPoint2d<LHS,RHS,OP> & Expression);
    
    /*! Sum operator */
    void operator+=(const point2d & V);
    
    /*! Subtract operator */
    void operator-=(const point2d & V);
    
    /*! Product operator */
    void operator*=(const Real & a);
    
    /*! Division operator */
    void operator/=(const Real & a);
    
    /*! Scalar product */
    Real operator*(const point2d &V) const;
    //@}
    
    
    /*! @name Ordinal functions */ //@{
  public:
    /*! Less operator: a vector is "less" than another if its first coordinate is less than the other vctor one.
    If the first component value is equal the second one is considered and so on */
    bool operator<(const point2d & V) const;
    
    /*! Inequality operator: two vectors are equal if their components are equal to the geometric tolerance */
    bool operator!=(const point2d & V) const;
    //@}
    
    
    /*! @name Other functions */ //@{
  public:
    /*! The norm of the R3 vector */
    Real norm2() const;
    
    /*! Static norm function */
    static Real norm2(const point2d & P);
    
    /*! Static dot function */
    static Real dot(const point2d & P1, const point2d & P2);
    
    /*! Clearing function */
    void clear();
    
    /*! Outstream operator */
    friend ostream & operator<<(ostream & f, const point2d & P);
    
    /*! Memory size */
    size_t memSize() const;
    //@}
};


//_________________________________________________________________________________________________
// INLINED AND TEMPLATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
Real 
point2d::
getX()
{ return(X[0]); }

const Real &
point2d::
getX() const
{ return(X[0]); }

Real
point2d::
getY()
{ return(X[1]); }

const Real &
point2d::
getY() const
{ return(X[1]); }

Real
point2d::
getI(UInt i)
{
  assert(i>=1);
  assert(i<=2);
  
  return(X[i-1]);
}

const Real &
point2d::
getI(const UInt & i) const
{
  assert(i>=1);
  assert(i<=2);
  
  return(X[i-1]);
}

template<typename LHS, typename RHS, typename OP>
point2d::
point2d(const expressionPoint2d<LHS,RHS,OP> & Expression)
{
  X[0] = Expression.getX();
  X[1] = Expression.getY();
}

template<typename LHS, typename RHS, typename OP>
void
point2d::
operator=(const expressionPoint2d<LHS,RHS,OP> & Expression)
{
  X[0] = Expression.getX();
  X[1] = Expression.getY();
}

template<class ARK>
void
point2d::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & X[0];
  ar & X[1];
  ar & id;
}

#endif
