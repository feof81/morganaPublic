/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef POINT3D_H
#define POINT3D_H

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
struct plusPoint3d
{
  static Real applyX(const LHS & left, const RHS & right) { return(left.getX() + right.getX()); }
  static Real applyY(const LHS & left, const RHS & right) { return(left.getY() + right.getY()); }
  static Real applyZ(const LHS & left, const RHS & right) { return(left.getZ() + right.getZ()); }
};

template<typename LHS, typename RHS>                       //The difference operation
struct minusPoint3d
{
  static Real applyX(const LHS & left, const RHS & right) { return(left.getX() - right.getX()); }
  static Real applyY(const LHS & left, const RHS & right) { return(left.getY() - right.getY()); }
  static Real applyZ(const LHS & left, const RHS & right) { return(left.getZ() - right.getZ()); }
};

template<typename LHS>                                     //The multiplication operation
struct multPoint3d
{
  static Real applyX(const LHS & left, const Real & right) { return(left.getX() * right); }
  static Real applyY(const LHS & left, const Real & right) { return(left.getY() * right); }
  static Real applyZ(const LHS & left, const Real & right) { return(left.getZ() * right); }
};

template<typename LHS>                                     //The multiplication operation
struct divPoint3d
{
  static Real applyX(const LHS & left, const Real & right) { return(left.getX() / right); }
  static Real applyY(const LHS & left, const Real & right) { return(left.getY() / right); }
  static Real applyZ(const LHS & left, const Real & right) { return(left.getZ() / right); }
};

template<typename LHS, typename RHS>                       //The rotor operation
struct rotPoint3d
{
  static Real applyX(const LHS & left, const RHS & right) { return(left.getY() * right.getZ() - right.getY() * left.getZ()); }
  static Real applyY(const LHS & left, const RHS & right) { return(left.getZ() * right.getX() - right.getZ() * left.getX()); }
  static Real applyZ(const LHS & left, const RHS & right) { return(left.getX() * right.getY() - right.getX() * left.getY()); }
};

template<typename LHS, typename RHS, typename OP>          //The expression template
struct expressionPoint3d
{
    static const morganaTypes myType = typePoint3d;
  
    const LHS & left;
    const RHS & right;
  
    expressionPoint3d(const LHS & Left, const RHS & Right) : left(Left), right(Right) { }
    Real getX() const { return(OP::applyX(left,right)); }
    Real getY() const { return(OP::applyY(left,right)); }
    Real getZ() const { return(OP::applyZ(left,right)); }
};



//_________________________________________________________________________________________________
// OPERATORS
//-------------------------------------------------------------------------------------------------
template<typename LHS, typename RHS>                       //The rotor operator
expressionPoint3d<LHS,RHS,rotPoint3d<LHS,RHS> >
operator^(const LHS & left, const RHS & right)
{
  return(expressionPoint3d<LHS,RHS,rotPoint3d<LHS,RHS> >(left,right));
}


//_________________________________________________________________________________________________
// POINT3D CLASS
//-------------------------------------------------------------------------------------------------
/*! The three dimensional point */
class point3d
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
  
  
    /*! @name Internal data */ //@{
  public:
    static const morganaTypes myType = typePoint3d;
    
    /// The coordinates
    Real X[3];
    
    /// The node id
    UInt id;
    //@}
    
  
    /*! @name Constructors */ //@{
  public:
    /*! Constructor */
    point3d(Real xx=0.0, Real yy=0.0, Real zz=0.0);
    
    /*! Constructor */
    template<typename LHS, typename RHS, typename OP>
    point3d(const expressionPoint3d<LHS,RHS,OP> & Expression);
    
    /*! Copy constructor */
    point3d(const point3d & V);

    /*! Destructor */
    virtual ~point3d();
    //@}
    
    
    /*! @name Set functions */ //@{
  public:
    /*! Set all components of the point */
    void set(const Real & xx, const Real & yy, const Real & zz);
    
    /*! Set the X component */
    void setX(const Real & xx);
    
    /*! Set the Y component */
    void setY(const Real & yy);

    /*! Set the Y component */
    void setZ(const Real & zz);

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

    /*! Get the Z component */
    inline Real getZ();
    inline const Real & getZ() const;

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
    point3d & operator=(const point3d & V);
    
    /*! The equality operator for expressions */
    template<typename LHS, typename RHS, typename OP>
    void operator=(const expressionPoint3d<LHS,RHS,OP> & Expression);
    
    /*! Sum operator */
    void operator+=(const point3d & V);
    
    /*! Subtract operator */
    void operator-=(const point3d & V);
    
    /*! Product operator */
    void operator*=(const Real & a);
    
    /*! Division operator */
    void operator/=(const Real & a);
    
    /*! Scalar product */
    Real operator*(const point3d &V) const;
    
    /*! Curl operator */
    void rotor(const point3d & Px, const point3d & Py, const point3d & Pz);
    //@}
    
    
    /*! @name Ordinal functions */ //@{
  public:
    /*! Less operator: a vector is "less" than another if its first coordinate is less than the other vctor one.
    If the first component value is equal the second one is considered and so on */
    bool operator<(const point3d & V) const;
    
    /*! Inequality operator: two vectors are equal if their components are equal to the geometric tolerance */
    bool operator!=(const point3d & V) const;
    
    /*! Equality operator */
    bool operator==(const point3d & V) const;
    //@}
    
    
    /*! @name Combinatorial functions */ //@{
  public:
    point3d combinationGS(const point3d & N) const;
    //@}
    
    
    /*! @name Other functions */ //@{
  public:
    /*! The norm of the R3 vector */
    Real norm2() const;
    
    /*! Static norm function */
    static Real norm2(const point3d & P);
    
    /*! Static dot function */
    static Real dot(const point3d & P1, const point3d & P2);
    
    /*! Clearing function */
    void clear();
    
    /*! Given a unit vector creates two other unit vectors that form an orthonormal basis */
    std::pair<point3d,point3d> orthoBasis();
    
    /*! Outstream operator */
    friend ostream & operator<<(ostream & f, const point3d & P);
    //@}
};


//_________________________________________________________________________________________________
// INLINED AND TEMPLATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
Real 
point3d::
getX()
{ return(X[0]); }

const Real &
point3d::
getX() const
{ return(X[0]); }

Real
point3d::
getY()
{ return(X[1]); }

const Real &
point3d::
getY() const
{ return(X[1]); }

Real
point3d::
getZ()
{ return(X[2]); }

const Real &
point3d::
getZ() const
{ return(X[2]); }

Real
point3d::
getI(UInt i)
{
  assert(i>=1);
  assert(i<=3);
  
  return(X[i-1]);
}

const Real &
point3d::
getI(const UInt & i) const
{
  assert(i>=1);
  assert(i<=3);
  
  return(X[i-1]);
}

template<typename LHS, typename RHS, typename OP>
point3d::
point3d(const expressionPoint3d<LHS,RHS,OP> & Expression)
{
  X[0] = Expression.getX();
  X[1] = Expression.getY();
  X[2] = Expression.getZ();
}

template<typename LHS, typename RHS, typename OP>
void
point3d::
operator=(const expressionPoint3d<LHS,RHS,OP> & Expression)
{
  X[0] = Expression.getX();
  X[1] = Expression.getY();
  X[2] = Expression.getZ();
} 

template<class ARK>
void
point3d::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & X[0];
  ar & X[1];
  ar & X[2];
  ar & id;
}

#endif
