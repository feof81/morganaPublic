/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef STATEMATRIX_H
#define STATEMATRIX_H

#include <assert.h>
#include <iostream>
#include <fstream>

#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>

#include "morganaTypes.hpp"

#include "stateVector.h"
#include "Epetra_SerialDenseMatrix.h"

using namespace std;



//_________________________________________________________________________________________________
// EXPRESSIONS
//-------------------------------------------------------------------------------------------------
template<typename LHS, typename RHS>
struct plusStateMatrix                                     //The sum operation
{
  static Real apply(const UInt & i, const UInt & j, const LHS & left, const RHS & right) 
  {
    assert(left.RowDim() == right.RowDim());
    assert(left.ColDim() == right.ColDim());
    return(left(i,j) + right(i,j));
  }
};

template<typename LHS, typename RHS>
struct minusStateMatrix                                    //The minus operation
{
  static Real apply(const UInt & i, const UInt & j, const LHS & left, const RHS & right)
  {
    assert(left.RowDim() == right.RowDim());
    assert(left.ColDim() == right.ColDim());
    return(left(i,j) - right(i,j));
  }
};

template<typename LHS>
struct multStateMatrix                                     //The sum operation
{
  static Real apply(const UInt & i, const UInt & j, const LHS & left, const Real & right)
  { return(left(i,j) * right); }
};

template<typename LHS>
struct divStateMatrix                                      //The div operation
{
  static Real apply(const UInt & i, const UInt & j, const LHS & left, const Real & right)
  { return(left(i,j) / right); }
};

template<typename LHS, typename RHS, typename OP>          //The expression template
struct expressionStateMatrix
{
    static const morganaTypes myType = typeStateMatrix;
  
    const LHS & left;
    const RHS & right;
  
    expressionStateMatrix(const LHS & Left, const RHS & Right) : left(Left), right(Right) { }
    Real operator() (const UInt & i, const UInt & j) const { return( OP::apply(i,j,left,right) ); }
    int  RowDim() const {return(left.RowDim());}
    int  ColDim() const {return(left.ColDim());}
};


//_________________________________________________________________________________________________
// STATE MATRIX CLASS
//-------------------------------------------------------------------------------------------------
/*! Dynamically sized dense matrix */
class stateMatrix : public Epetra_SerialDenseMatrix
{
  /*! @name Parallel support */ //@{
  private:
    friend class boost::serialization::access;
    
    int srlzRows; //Serialization rows
    int srlzCols; //Serialization cols
    
    /*! Serialization function */
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    
    /*! @name Internal data */ //@{
  public:
    static const morganaTypes myType = typeStateMatrix;
    //@}
    
    
    /*! @name Constructors */ //@{
  public:
    /*! Constructor */
    stateMatrix();
    
    /*! Constructor */
    stateMatrix(const UInt & NRows, const UInt & NCols);
    
    /*! Copy constructor */
    stateMatrix(const stateMatrix &M);
    
    /*! Constructor for expression templates */
    template<typename LHS, typename RHS, typename OP>
    stateMatrix(const expressionStateMatrix<LHS,RHS,OP> & Expression);
    //@}
    
    
    /*! @name Operators */ //@{
  public:
    /*! Access operator */
    inline Real & operator()(const UInt & i, const UInt & j);
    
    /*! Access operator */
    inline const Real & operator()(const UInt & i, const UInt & j) const;
    
    /*! The equality operator for expression templates */
    template<typename LHS, typename RHS, typename OP>
    void operator=(const expressionStateMatrix<LHS,RHS,OP> & Expression);
    
    /*! Add operator */ 
    void operator+=(const stateMatrix & M);
    
    /*! Subtract operator */
    void operator-=(const stateMatrix & M);
    
    /*! Multiplication operator */
    void operator*=(const Real &a);
    
    /*! Division operator */
    void operator/=(const Real &a);
    //@}
    
    
    /*! @name Set-Get functions */ //@{
  public:
    /*! Set i-th row */
    void setRow(const UInt & i, const stateVector & V);
      
    /*! Set i-th coloumn */
    void setCol(const UInt & j, const stateVector & V);
    
    /*! Get i-th row  */
    stateVector getRow(const UInt & i);
      
    /*! Get i-th coloumn */
    stateVector getCol(const UInt & j);
    //@}
    
    
    /*! @name Ordinal functions */ //@{
  public:
    /*! Less operator: a vector is "less" than another if its first coordinate is less than the other vctor one.
    If the first component value is equal the second one is considered and so on */
    bool operator<(const stateMatrix & V) const;
    
    /*! Inequality operator: two vectors are equal if their components are equal to the geometric tolerance */
    bool operator!=(const stateMatrix & V) const;
    //@}
      
      
    /*! @name Other functions */ //@{
  public:
    /*! Outstream operator */
    friend ostream & operator<<( ostream & f, const stateMatrix & A);
    
    /*! Memory size */
    size_t memSize() const;
    //@}
};



//_________________________________________________________________________________________________
// TEMPLATED AND INLINED FUNCTIONS
//-------------------------------------------------------------------------------------------------
Real & 
stateMatrix::
operator()(const UInt & i, const UInt & j)
{
  return(Epetra_SerialDenseMatrix::operator()(i-1,j-1));
}

const Real & 
stateMatrix::
operator()(const UInt & i, const UInt & j) const
{
  return(Epetra_SerialDenseMatrix::operator()(i-1,j-1));
}

template<typename LHS, typename RHS, typename OP>
stateMatrix::
stateMatrix(const expressionStateMatrix<LHS,RHS,OP> & Expression)
{
  if( (Expression.RowDim() != this->RowDim()) || (Expression.ColDim() != this->ColDim()) )
  { 
    this->Reshape(Expression.RowDim(),Expression.ColDim());
  }
  
  for(int i=1; i <= this->RowDim(); ++i)
  {
    for(int j=1; j <= this->ColDim(); ++j)
    {
      this->operator()(i,j) = Expression(i,j);
    }
  }
}

template<typename LHS, typename RHS, typename OP>
void
stateMatrix::
operator=(const expressionStateMatrix<LHS,RHS,OP> & Expression)
{
  if( (Expression.RowDim() != this->RowDim()) || (Expression.ColDim() != this->ColDim()) )
  { 
    this->Reshape(Expression.RowDim(),Expression.ColDim());
  }
  
  for(int i=1; i <= this->RowDim(); ++i)
  {
    for(int j=1; j <= this->ColDim(); ++j)
    {
      this->operator()(i,j) = Expression(i,j);
    }
  }
}

template<class ARK>
void
stateMatrix::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  srlzRows = this->RowDim();
  srlzCols = this->ColDim();
  
  ar & srlzRows;
  ar & srlzCols;
  
  if( (srlzRows != this->RowDim()) || (srlzCols != this->ColDim()) )
  { 
    this->Reshape(srlzRows,srlzCols);
  }
  
  for(int i=1; i <= srlzRows; ++i)
  {
    for(int j=1; j <= srlzCols; ++j)
    {
      ar & this->operator()(i,j);
    }
  }
}

#endif

