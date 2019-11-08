/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef TENSOR2D_H
#define TENSOR2D_H

#include <assert.h>
#include <iostream>
#include <fstream>

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "point2d.h"
#include "morganaTypes.hpp"

using namespace std;



//_________________________________________________________________________________________________
// EXPRESSIONS
//-------------------------------------------------------------------------------------------------
template<typename LHS, typename RHS>
struct plusTensor2d                                        //The sum operation
{
  static Real apply(const UInt & i, const UInt & j, const LHS & left, const RHS & right)
  { return(left.getIJ(i,j) + right.getIJ(i,j)); }
};

template<typename LHS, typename RHS>
struct minusTensor2d                                       //The minus operation
{
  static Real apply(const UInt & i, const UInt & j, const LHS & left, const RHS & right)
  { return(left.getIJ(i,j) - right.getIJ(i,j)); }
};

template<typename LHS>
struct multTensor2d                                        //The sum operation
{
  static Real apply(const UInt & i, const UInt & j, const LHS & left, const Real & right)
  { return(left.getIJ(i,j) * right); }
};

template<typename LHS>
struct divTensor2d                                         //The div operation
{
  static Real apply(const UInt & i, const UInt & j, const LHS & left, const Real & right)
  { return(left.getIJ(i,j) / right); }
};

template<typename LHS, typename RHS, typename OP>          //The expression template
struct expressionTensor2d
{
    static const morganaTypes myType = typeTensor2d;
  
    const LHS & left;
    const RHS & right;
  
    expressionTensor2d(const LHS & Left, const RHS & Right) : left(Left), right(Right) { }
    Real getIJ(const UInt & i, const UInt & j) const { return( OP::apply(i,j,left,right) ); }
};



//_________________________________________________________________________________________________
// TENSOR3D CLASS
//-------------------------------------------------------------------------------------------------

/*! The two dimensional tensor */
class tensor2d
{
  /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    
    /*! @name Internal data */ //@{
  public:
    static const morganaTypes myType = typeTensor2d;
    
    /// Components
    Real T[2][2];
      
    /// The id
    UInt id;
    //@}
    
    
    /*! @name Constructors */ //@{
  public:
    /*! Constructor */
    tensor2d();
      
    /*! Constructor */
    template<typename LHS, typename RHS, typename OP>
    tensor2d(const expressionTensor2d<LHS,RHS,OP> & Expression);
      
    /*! Copy constructor */
    tensor2d(const tensor2d & v);
    
    /*! Destructor */
    virtual ~tensor2d();
    //@}
    
    
    /*! @name Set functions */ //@{
  public:
    /*! Set the i-j component*/
    void setIJ(const UInt & i, const UInt & j, const Real & value);
      
    /*! Sum into the i-j component*/
    void sumIJ(const UInt & i, const UInt & j, const Real & value);
      
    /*! Set dell'i-esima riga */
    void setRow(const UInt & i, const point2d & value);
      
    /*! Set dell'i-esima colonna */
    void setCol(const UInt & j , const point2d & value);
    
    /*! Sum the i-th row */
    void sumRow(const UInt & i, const point2d & value);
    
    /*! Sum the i-th coloumn */
    void sumCol(const UInt & j, const point2d & value);
      
    /*! Set dell'Id */
    void setId(const UInt & Id);
    //@}
    
    
    /*! @name Get functions */ //@{
  public:
    /*! Set the i-j component*/
    inline Real getIJ(const UInt & i, const UInt & j);
    inline const Real & getIJ(const UInt & i, const UInt & j) const;
      
    /*! Set dell'i-esima riga */
    point2d getRow(const UInt & i);
      
    /*! Set dell'i-esima colonna */
    point2d getCol(const UInt & j);
      
    /*! Set dell'Id */
    UInt getId() const;
    //@}
    
    
    /*! @name Operators */ //@{
  public:
    /*! Access operator */    
    Real operator()(const UInt & i, const UInt & j);
    
    /*! Access operator */
    const Real & operator()(const UInt & i, const UInt & j) const;
    
    /*! The equality operator */
    tensor2d & operator=(const tensor2d & v);
    
    /*! The equality operator for expressions */
    template<typename LHS, typename RHS, typename OP>
    void operator=(const expressionTensor2d<LHS,RHS,OP> & Expression);
    
    /*! Sum operator */
    void operator+=(const tensor2d & V);
    
    /*! Subtract operator */
    void operator-=(const tensor2d & V);
    
    /*! Product operator */
    void operator*=(const Real & a);
    
    /*! Division operator */
    void operator/=(const Real & a);
    //@}
    
    
    /*! @name Internal functions */ //@{
  public:  
    /*! Compute the matrix inverse */
    void computeInverse();
      
    /*! Compute the matrix traspose */
    void transpose();
    //@}
    
    
    /*! @name External functions */ //@{
  public:
    /*! Frobenius norm*/
    Real getFrobenius() const;
      
    /*! Left product*/
    point2d firstIndexSaturation(const point2d & v) const;
      
    /*! Right product*/
    point2d secondIndexSaturation(const point2d & v) const;
      
    /*! Prodotto diretto a * b per produrre un tensore */
    void directProduct(const point2d & a, const point2d & b);
      
    /*! Matrix product */
    tensor2d matrixProduct(const tensor2d & A) const;
      
    /*! Prodotto scalare fra tensori this*A */
    Real scalarProduct(const tensor2d & A) const;
    //@}
    
    
    /*! @name Ordinal functions */ //@{
  public:
    /*! Less operator: a vector is "less" than another if its first coordinate is less than the other vctor one.
    If the first component value is equal the second one is considered and so on */
    bool operator<(const tensor2d & V) const;
    
    /*! Inequality operator: two vectors are equal if their components are equal to the geometric tolerance */
    bool operator!=(const tensor2d & V) const;
    //@}
    
    /*! @name Other functions */ //@{
  public:
    /*! Outstream operator */
    friend ostream & operator << ( ostream & f, const tensor2d & A);
    //@}
};



//_________________________________________________________________________________________________
// TEMPLATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename LHS, typename RHS, typename OP>
tensor2d::
tensor2d(const expressionTensor2d<LHS,RHS,OP> & Expression)
{
  T[0][0] = Expression.getIJ(1,1);
  T[0][1] = Expression.getIJ(1,2);
  
  T[1][0] = Expression.getIJ(2,1);
  T[1][1] = Expression.getIJ(2,2);
}

template<typename LHS, typename RHS, typename OP>
void
tensor2d::
operator=(const expressionTensor2d<LHS,RHS,OP> & Expression)
{
  T[0][0] = Expression.getIJ(1,1);
  T[0][1] = Expression.getIJ(1,2);
  
  T[1][0] = Expression.getIJ(2,1);
  T[1][1] = Expression.getIJ(2,2);
}

template<class ARK>
void
tensor2d::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & T[0][0];
  ar & T[0][1];
  
  ar & T[1][0];
  ar & T[1][1];
    
  ar & id;
}

Real
tensor2d::
getIJ(const UInt & i, const UInt & j)
{
  assert((i>=1) & (j>=1));
  assert((i<=2) & (j<=2));
  
  return(T[i-1][j-1]);
}

const Real &
tensor2d::
getIJ(const UInt & i, const UInt & j) const
{
  assert((i>=1) & (j>=1));
  assert((i<=2) & (j<=2));
  
  return(T[i-1][j-1]);
}

#endif
