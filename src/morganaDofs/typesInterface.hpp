/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef TYPESINTERFACE_HPP
#define TYPESINTERFACE_HPP

#include "staticAssert.hpp"

#include "point2d.h"
#include "point3d.h"
#include "tensor2d.h"
#include "tensor3d.h"
#include "stateVector.h"
#include "stateMatrix.h"
#include "staticVector.hpp"
#include "komplex.h"

#include "traitsMpiOptimization.hpp"



//_________________________________________________________________________________________________
// PLUS TRAITS AND OPERATOR : TYPE + TYPE
//-------------------------------------------------------------------------------------------------
/*! Empty plus trait */
template<typename LHS, typename RHS, morganaTypes>         //Empty trait
class plusTrait;

/*! Point 2d plus trait */
template<typename LHS, typename RHS>                       //Point2d trait
class plusTrait<LHS,RHS,typePoint2d>
{
  public:
    typedef expressionPoint2d<LHS,RHS,plusPoint2d<LHS,RHS> > EXPRESSION;
};

/*! Point 3d plus trait */
template<typename LHS, typename RHS>                       //Point3d trait
class plusTrait<LHS,RHS,typePoint3d>
{
  public:
    typedef expressionPoint3d<LHS,RHS,plusPoint3d<LHS,RHS> > EXPRESSION;
};

/*! Tensor 3d plus trait */
template<typename LHS, typename RHS>                       //Tensor3d trait
class plusTrait<LHS,RHS,typeTensor2d>
{
  public:
    typedef expressionTensor2d<LHS,RHS,plusTensor2d<LHS,RHS> > EXPRESSION;
};

/*! Tensor 2d plus trait */
template<typename LHS, typename RHS>                       //Tensor3d trait
class plusTrait<LHS,RHS,typeTensor3d>
{
  public:
    typedef expressionTensor3d<LHS,RHS,plusTensor3d<LHS,RHS> > EXPRESSION;
};

/*! StateVector plus trait */
template<typename LHS, typename RHS>                       //StateVector trait
class plusTrait<LHS,RHS,typeStateVector>
{
  public:
    typedef expressionStateVector<LHS,RHS,plusStateVector<LHS,RHS> > EXPRESSION;
};

/*! State Matrix plus trait */
template<typename LHS, typename RHS>                       //StateMatrix trait
class plusTrait<LHS,RHS,typeStateMatrix>
{
  public:
    typedef expressionStateMatrix<LHS,RHS,plusStateMatrix<LHS,RHS> > EXPRESSION;
};

/*! Static Vector plus trait */
template<typename LHS, typename RHS>                       //StaticVector trait
class plusTrait<LHS,RHS,typeStaticVector>
{
  public:
    typedef expressionStaticVector<LHS,RHS,plusStaticVector<LHS,RHS> > EXPRESSION;
};

/*! Komplex plus trait */
template<typename LHS, typename RHS>                       //Komplex trait
class plusTrait<LHS,RHS,typeKomplex>
{
  public:
    typedef expressionKomplex<LHS,RHS,plusKomplex<LHS,RHS> > EXPRESSION;
};

/*! The sum operator */
template<typename LHS, typename RHS>                       //The sum operator
typename plusTrait<LHS,RHS,LHS::myType>::EXPRESSION
operator+(const LHS & left, const RHS & right)
{
  typedef typename plusTrait<LHS,RHS,LHS::myType>::EXPRESSION  EXPRESSION;
  assert(staticAssert<LHS::myType == RHS::myType>::returnValue);
  return(EXPRESSION(left,right));
}



//_________________________________________________________________________________________________
// MINUS TRAITS AND OPERATOR : TYPE - TYPE
//-------------------------------------------------------------------------------------------------
/*! Empty minus trait */
template<typename LHS, typename RHS, morganaTypes>         //Empty trait
class minusTrait;

/*! Point 2d minus trait */
template<typename LHS, typename RHS>                       //Point2d trait
class minusTrait<LHS,RHS,typePoint2d>
{
  public:
    typedef expressionPoint2d<LHS,RHS,minusPoint2d<LHS,RHS> > EXPRESSION;
};

/*! Point 3d minus trait */
template<typename LHS, typename RHS>                       //Point3d trait
class minusTrait<LHS,RHS,typePoint3d>
{
  public:
    typedef expressionPoint3d<LHS,RHS,minusPoint3d<LHS,RHS> > EXPRESSION;
};

/*! Tensor 2d minus trait */
template<typename LHS, typename RHS>                       //Tensor2d trait
class minusTrait<LHS,RHS,typeTensor2d>
{
  public:
    typedef expressionTensor2d<LHS,RHS,minusTensor2d<LHS,RHS> > EXPRESSION;
};

/*! Tensor 3d minus trait */
template<typename LHS, typename RHS>                       //Tensor3d trait
class minusTrait<LHS,RHS,typeTensor3d>
{
  public:
    typedef expressionTensor3d<LHS,RHS,minusTensor3d<LHS,RHS> > EXPRESSION;
};

/*! State Vector minus trait */
template<typename LHS, typename RHS>                       //StateVector trait
class minusTrait<LHS,RHS,typeStateVector>
{
  public:
    typedef expressionStateVector<LHS,RHS,minusStateVector<LHS,RHS> > EXPRESSION;
};

/*! State Matrix minus trait */
template<typename LHS, typename RHS>                       //StateMatrix trait
class minusTrait<LHS,RHS,typeStateMatrix>
{
  public:
    typedef expressionStateMatrix<LHS,RHS,minusStateMatrix<LHS,RHS> > EXPRESSION;
};

/*! Static Vector minus trait */
template<typename LHS, typename RHS>                       //StaticVector trait
class minusTrait<LHS,RHS,typeStaticVector>
{
  public:
    typedef expressionStaticVector<LHS,RHS,minusStaticVector<LHS,RHS> > EXPRESSION;
};

/*! Komplex minus trait */
template<typename LHS, typename RHS>                       //Komplex trait
class minusTrait<LHS,RHS,typeKomplex>
{
  public:
    typedef expressionKomplex<LHS,RHS,minusKomplex<LHS,RHS> > EXPRESSION;
};

/*! Difference operator */
template<typename LHS, typename RHS>                       //The difference operator
typename minusTrait<LHS,RHS,LHS::myType>::EXPRESSION
operator-(const LHS & left, const RHS & right)
{
  typedef typename minusTrait<LHS,RHS,LHS::myType>::EXPRESSION               EXPRESSION;
  assert(staticAssert<LHS::myType == RHS::myType>::returnValue);
  return(EXPRESSION(left,right));
}



//_________________________________________________________________________________________________
// MULT TRAITS AND OPERATOR : TYPE - TYPE
//-------------------------------------------------------------------------------------------------
/*! Empty multiplication trait */
template<typename LHS, morganaTypes>                       //Empty trait
class multTrait;

/*! Point 2d multiplication trait */
template<typename LHS>                                     //Point2d trait
class multTrait<LHS,typePoint2d>
{
  public:
    typedef expressionPoint2d<LHS,Real,multPoint2d<LHS> > EXPRESSION;
};

/*! Point 3d multiplication trait */
template<typename LHS>                                     //Point3d trait
class multTrait<LHS,typePoint3d>
{
  public:
    typedef expressionPoint3d<LHS,Real,multPoint3d<LHS> > EXPRESSION;
};

/*! Tensor 2d multiplication trait */
template<typename LHS>                                     //Tensor2d trait
class multTrait<LHS,typeTensor2d>
{
  public:
    typedef expressionTensor2d<LHS,Real,multTensor2d<LHS> > EXPRESSION;
};

/*! Tensor 3d multiplication trait */
template<typename LHS>                                     //Tensor3d trait
class multTrait<LHS,typeTensor3d>
{
  public:
    typedef expressionTensor3d<LHS,Real,multTensor3d<LHS> > EXPRESSION;
};

/*! State Vector multiplication trait */
template<typename LHS>                                     //StateVector trait
class multTrait<LHS,typeStateVector>
{
  public:
    typedef expressionStateVector<LHS,Real,multStateVector<LHS> > EXPRESSION;
};

/*! State Matrix multiplication trait */
template<typename LHS>                                     //StateMatrix trait
class multTrait<LHS,typeStateMatrix>
{
  public:
    typedef expressionStateMatrix<LHS,Real,multStateMatrix<LHS> > EXPRESSION;
};

/*! Static Vector multiplication trait */
template<typename LHS>                                     //StaticVector trait
class multTrait<LHS,typeStaticVector>
{
  public:
    typedef expressionStaticVector<LHS,Real,multStaticVector<LHS> > EXPRESSION;
};

/*! Komplex multiplication trait */
template<typename LHS>                                     //Komplex trait
class multTrait<LHS,typeKomplex>
{
  public:
    typedef expressionKomplex<LHS,Real,multKomplex<LHS> > EXPRESSION;
};

/*! Multiplication operator */
template<typename LHS>                                     //The multiplication operator
typename  multTrait<LHS,LHS::myType>::EXPRESSION
operator*(const LHS & left, const Real & right)
{
  typedef typename multTrait<LHS,LHS::myType>::EXPRESSION EXPRESSION;
  return(EXPRESSION(left,right));
}


//_________________________________________________________________________________________________
// DIV TRAITS AND OPERATOR : TYPE - TYPE
//-------------------------------------------------------------------------------------------------
/*! Empty division trait */
template<typename LHS, morganaTypes>                       //Empty trait
class divTrait;

/*! Point 2d division trait */
template<typename LHS>                                     //Point2d trait
class divTrait<LHS,typePoint2d>
{
  public:
    typedef expressionPoint2d<LHS,Real,divPoint2d<LHS> > EXPRESSION;
};

/*! Point 3d division trait */
template<typename LHS>                                     //Point3d trait
class divTrait<LHS,typePoint3d>
{
  public:
    typedef expressionPoint3d<LHS,Real,divPoint3d<LHS> > EXPRESSION;
};

/*! Tensor 2d division trait */
template<typename LHS>                                     //Tensor2d trait
class divTrait<LHS,typeTensor2d>
{
  public:
    typedef expressionTensor2d<LHS,Real,divTensor2d<LHS> > EXPRESSION;
};

/*! Tensor 3d division trait */
template<typename LHS>                                     //Tensor3d trait
class divTrait<LHS,typeTensor3d>
{
  public:
    typedef expressionTensor3d<LHS,Real,divTensor3d<LHS> > EXPRESSION;
};

/*! State Vector division trait */
template<typename LHS>                                     //StateVector trait
class divTrait<LHS,typeStateVector>
{
  public:
    typedef expressionStateVector<LHS,Real,divStateVector<LHS> > EXPRESSION;
};

/*! State Matrix division trait */
template<typename LHS>                                     //StateMatrix trait
class divTrait<LHS,typeStateMatrix>
{
  public:
    typedef expressionStateMatrix<LHS,Real,divStateMatrix<LHS> > EXPRESSION;
};

/*! Static Vector division trait */
template<typename LHS>                                     //StaticVector trait
class divTrait<LHS,typeStaticVector>
{
  public:
    typedef expressionStaticVector<LHS,Real,divStaticVector<LHS> > EXPRESSION;
};

/*! Komplex division trait */
template<typename LHS>                                     //Komplex trait
class divTrait<LHS,typeKomplex>
{
  public:
    typedef expressionKomplex<LHS,Real,divKomplex<LHS> > EXPRESSION;
};

/*! Division operator */
template<typename LHS>                                     //The division operator
typename  divTrait<LHS,LHS::myType>::EXPRESSION
operator/(const LHS & left, const Real & right)
{
  typedef typename divTrait<LHS,LHS::myType>::EXPRESSION EXPRESSION;
  return(EXPRESSION(left,right));
}


#endif
