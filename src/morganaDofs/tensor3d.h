/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef TENSOR3D_H
#define TENSOR3D_H

#include <assert.h>
#include <iostream>
#include <fstream>

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "point3d.h"
#include "morganaTypes.hpp"

using namespace std;



//_________________________________________________________________________________________________
// EXPRESSIONS
//-------------------------------------------------------------------------------------------------
template<typename LHS, typename RHS>
struct plusTensor3d                                        //The sum operation
{
  static Real apply(const UInt & i, const UInt & j, const LHS & left, const RHS & right) 
  { return(left.getIJ(i,j) + right.getIJ(i,j)); }
};

template<typename LHS, typename RHS>
struct minusTensor3d                                       //The minus operation
{
  static Real apply(const UInt & i, const UInt & j, const LHS & left, const RHS & right) 
  { return(left.getIJ(i,j) - right.getIJ(i,j)); }
};

template<typename LHS>
struct multTensor3d                                        //The sum operation
{
  static Real apply(const UInt & i, const UInt & j, const LHS & left, const Real & right)
  { return(left.getIJ(i,j) * right); }
};

template<typename LHS>
struct divTensor3d                                         //The div operation
{
  static Real apply(const UInt & i, const UInt & j, const LHS & left, const Real & right)
  { return(left.getIJ(i,j) / right); }
};

template<typename LHS, typename RHS, typename OP>          //The expression template
struct expressionTensor3d
{
    static const morganaTypes myType = typeTensor3d;
  
    const LHS & left;
    const RHS & right;
  
    expressionTensor3d(const LHS & Left, const RHS & Right) : left(Left), right(Right) { }
    Real getIJ(const UInt & i, const UInt & j) const { return( OP::apply(i,j,left,right) ); }
};




//_________________________________________________________________________________________________
// TENSOR3D CLASS
//-------------------------------------------------------------------------------------------------
/*! The three dimensional tensor */
class tensor3d
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    
    /*! @name Dati interni */ //@{
  public:
    static const morganaTypes myType = typeTensor3d;
    
    /// Components
    Real T[3][3];
      
    /// The id
    UInt id;
    //@}
    
    
    /*! @name Constructors */ //@{
  public:
    /*! Constructor */
    tensor3d();
    
    /*! Constructor */
    tensor3d(const Real & XX, const Real & XY, const Real & XZ,
	     const Real & YX, const Real & YY, const Real & YZ,
	     const Real & ZX, const Real & ZY, const Real & ZZ);

    tensor3d(const point3d & P1, const point3d & P2, const point3d & P3);
      
    /*! Constructor */
    template<typename LHS, typename RHS, typename OP>
    tensor3d(const expressionTensor3d<LHS,RHS,OP> & Expression);
      
    /*! Copy constructor */
    tensor3d(const tensor3d & v);
    
    /*! Destructor */
    virtual ~tensor3d();
    //@}
    
    
    /*! @name Set functions */ //@{
  public:
    /*! Set the i-j component*/
    void setIJ(const UInt & i, const UInt & j, const Real & value);
      
    /*! Sum into the i-j component*/
    void sumIJ(const UInt & i, const UInt & j, const Real & value);
      
    /*! Set dell'i-esima riga */
    void setRow(const UInt & i, const point3d & value);
      
    /*! Set dell'i-esima colonna */
    void setCol(const UInt & j , const point3d & value);
    
    /*! Sum the i-th row */
    void sumRow(const UInt & i, const point3d & value);
    
    /*! Sum the i-th coloumn */
    void sumCol(const UInt & j, const point3d & value);
      
    /*! Set dell'Id */
    void setId(const UInt & Id);
    //@}
    
    
    /*! @name Get functions */ //@{
  public:
    /*! Set the i-j component*/
    inline Real getIJ(const UInt & i, const UInt & j);
    inline const Real & getIJ(const UInt & i, const UInt & j) const;
      
    /*! Set i-th row */
    point3d getRow(const UInt & i) const;
      
    /*! Set i-th coloumn */
    point3d getCol(const UInt & j) const;
      
    /*! Set Id */
    UInt getId() const;
    //@}
    
    
    /*! @name Operators */ //@{
  public:
    /*! Access operator */    
    Real & operator()(const UInt & i, const UInt & j);
    
    /*! Access operator */
    const Real & operator()(const UInt & i, const UInt & j) const;
    
    /*! The equality operator */
    tensor3d & operator=(const tensor3d & v);
    
    /*! The equality operator for expressions */
    template<typename LHS, typename RHS, typename OP>
    void operator=(const expressionTensor3d<LHS,RHS,OP> & Expression);
    
    /*! Sum operator */
    void operator+=(const tensor3d & V);
    
    /*! Subtract operator */
    void operator-=(const tensor3d & V);
    
    /*! Product operator */
    void operator*=(const Real & a);
    
    /*! Division operator */
    void operator/=(const Real & a);
    //@}
    
    
    /*! @name Internal functions */ //@{
  public:  
    /*! Completes the second coloumn with a unit vector orthogonal to the second coloumn */
    void completeSecondColoumn();
     
    /*! Completes the second row with a unit vector orthogonal to the second row */
    void completeSecondRow();
      
    /*! Completes the third coloumn with the cross product of the first two coloumns */
    void completeThirdColoumn();

    /*! Completes the third row with the cross product of the first two rows */
    void completeThirdRow();
    
    /*! Computes the matrix inverse */
    void computeInverse();
      
    /*! Computes the matrix traspose */
    void transpose();
    //@}
    
    
    /*! @name External functions */ //@{
  public:
    /*! First invariant*/
    Real getFirstInvariant() const;
      
    /*! Second invariant*/
    Real getSecondInvariant() const;
      
    /*! Third invariant*/
    Real getThirdInvariant() const;

    /*! Frobenius norm*/
    Real getFrobenius() const;
      
    /*! Left product*/
    point3d firstIndexSaturation(const point3d & v) const;
      
    /*! Right product*/
    point3d secondIndexSaturation(const point3d & v) const;
      
    /*! Prodotto diretto a * b per produrre un tensore */
    void directProduct(const point3d & a, const point3d & b);
      
    /*! Matrix product */
    tensor3d matrixProduct(const tensor3d & A) const;
      
    /*! Scalar product (this * A) */
    Real scalarProduct(const tensor3d & A) const;
      
    /*! Saturation with the Ricci tensor eps_{i j k} v_k */
    void ricciThirdSaturation(const point3d & v);
      
    /*! Computes the vector (v23-v32, v31-v13, v12-v21) */
    point3d antiSymmetric(const tensor3d & A) const;
    //@}
    
    
    /*! @name Ordinal functions */ //@{
  public:
    /*! Less operator: a vector is "less" than another if its first coordinate is less than the other vctor one.
    If the first component value is equal the second one is considered and so on */
    bool operator<(const tensor3d & V) const;
    
    /*! Inequality operator: two vectors are equal if their components are equal to the geometric tolerance */
    bool operator!=(const tensor3d & V) const;
    //@}
    
    /*! @name Other functions */ //@{
  public:
    /*! Outstream operator */
    friend ostream & operator << ( ostream & f, const tensor3d & A);
    //@}
};



//_________________________________________________________________________________________________
// TEMPLATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename LHS, typename RHS, typename OP>
tensor3d::
tensor3d(const expressionTensor3d<LHS,RHS,OP> & Expression)
{
  T[0][0] = Expression.getIJ(1,1);
  T[0][1] = Expression.getIJ(1,2);
  T[0][2] = Expression.getIJ(1,3);
  
  T[1][0] = Expression.getIJ(2,1);
  T[1][1] = Expression.getIJ(2,2);
  T[1][2] = Expression.getIJ(2,3);
  
  T[2][0] = Expression.getIJ(3,1);
  T[2][1] = Expression.getIJ(3,2);
  T[2][2] = Expression.getIJ(3,3);
}

template<typename LHS, typename RHS, typename OP>
void
tensor3d::
operator=(const expressionTensor3d<LHS,RHS,OP> & Expression)
{
  T[0][0] = Expression.getIJ(1,1);
  T[0][1] = Expression.getIJ(1,2);
  T[0][2] = Expression.getIJ(1,3);
  
  T[1][0] = Expression.getIJ(2,1);
  T[1][1] = Expression.getIJ(2,2);
  T[1][2] = Expression.getIJ(2,3);
  
  T[2][0] = Expression.getIJ(3,1);
  T[2][1] = Expression.getIJ(3,2);
  T[2][2] = Expression.getIJ(3,3);
}

template<class ARK>
void
tensor3d::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & T[0][0];
  ar & T[0][1];
  ar & T[0][2];
  
  ar & T[1][0];
  ar & T[1][1];
  ar & T[1][2];
  
  ar & T[2][0];
  ar & T[2][1];
  ar & T[2][2];
  
  ar & id;
}

Real
tensor3d::
getIJ(const UInt & i, const UInt & j)
{
  assert((i>=1) & (j>=1));
  assert((i<=3) & (j<=3));
  
  return(T[i-1][j-1]);
}

const Real &
tensor3d::
getIJ(const UInt & i, const UInt & j) const
{
  assert((i>=1) & (j>=1));
  assert((i<=3) & (j<=3));
  
  return(T[i-1][j-1]);
}

#endif
