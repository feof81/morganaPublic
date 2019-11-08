/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef TRAITSMULTIPLY_H
#define TRAITSMULTIPLY_H

#include "typesInterface.hpp"

/*! Traits for dof multiplications, the following operations are supported

Direct Product
<ol>
<li> Real         * Real         = Real
<li> Real         * point2d      = point2d
<li> Real         * point3d      = point3d
<li> Real         * tensor2d     = tensor2d
<li> Real         * tensor3d     = tensor3d
<li> Real         * stateVector  = stateVector
<li> Real         * stateMatrix  = stateMatrix
<li> Real         * staticVector = staticVector
<li> Real         * komplex      = komplex
<li> point2d      * Real         = point2d
<li> point2d      * point2d      = tensor2d
<li> point3d      * Real         = point3d
<li> point3d      * point3d      = tensor3d
<li> tensor2d     * Real         = tensor2d
<li> tensor3d     * Real         = tensor3d
<li> stateVector  * Real         = stateVector
<li> stateMatrix  * Real         = stateMatrix
<li> staticVector * Real         = staticVector
<li> komplex      * Real         = komplex
<li> komplex      * komplex      = komplex
</ol>

Cross product A
<ol>
<li> point2d  * point2d  = point3d  (E_{i j k} a_j     b_k     )
<li> point3d  * point3d  = point3d  (E_{i j k} a_j     b_k     )
</ol>

Scalar product A
<ol>
<li> Real      * Real      = Real
<li> point2d   * point2d   = Real   ( a_k     b_k)
<li> point3d   * point3d   = Real   ( a_k     b_k)
</ol>
*/
template<typename DOFA, typename DOFB>
class traitsMultiply;


//_________________________________________________________________________________________________
// REAL-REAL TRAIT
//-------------------------------------------------------------------------------------------------
/*! Real - Real */
template<> class traitsMultiply<Real,Real>
{
  public:
    typedef Real DIRECTTYPE;
    typedef Real SCALARTYPEA;
    
    traitsMultiply() { };
    static DIRECTTYPE  multiply(const Real & A, const Real & B);
    static SCALARTYPEA scalarProductA(const Real & A, const Real & B);
};


//_________________________________________________________________________________________________
// REAL-POINT2D TRAIT
//-------------------------------------------------------------------------------------------------
/*! Real - Point2d */
template<> class traitsMultiply<Real,point2d>
{
  public:
    typedef point2d DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const Real & A, const point2d & B);
};


//_________________________________________________________________________________________________
// REAL-POINT3D TRAIT
//-------------------------------------------------------------------------------------------------
/*! Real - Point3d */
template<> class traitsMultiply<Real,point3d>
{
  public:
    typedef point3d DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const Real & A, const point3d & B);
};


//_________________________________________________________________________________________________
// REAL-TENSOR2D TRAIT
//-------------------------------------------------------------------------------------------------
/*! Real - Tensor2d */
template<> class traitsMultiply<Real,tensor2d>
{
  public:
    typedef tensor2d DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const Real & A, const tensor2d & B);
};


//_________________________________________________________________________________________________
// REAL-TENSOR3D TRAIT
//-------------------------------------------------------------------------------------------------
/*! Real - Tensor3d */
template<> class traitsMultiply<Real,tensor3d>
{
  public:
    typedef tensor3d DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const Real & A, const tensor3d & B);
};


//_________________________________________________________________________________________________
// REAL-STATEVECTOR TRAIT
//-------------------------------------------------------------------------------------------------
/*! Real - StateVector */
template<> class traitsMultiply<Real,stateVector>
{
  public:
    typedef stateVector DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const Real & A, const stateVector & B);
};


//_________________________________________________________________________________________________
// REAL-STATEMATRIX TRAIT
//-------------------------------------------------------------------------------------------------
/*! Real - StateMatrix */
template<> class traitsMultiply<Real,stateMatrix>
{
  public:
    typedef stateMatrix DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const Real & A, const stateMatrix & B);
};


//_________________________________________________________________________________________________
// REAL-STATIC VECTOR
//-------------------------------------------------------------------------------------------------
/*! Real - StaticVector */
template<size_t N> class traitsMultiply<Real,staticVector<N> >
{
  public:
    typedef staticVector<N> DIRECTTYPE;
    
    traitsMultiply() { };
    static staticVector<N> multiply(const Real & A, const staticVector<N> & B);
};

template<size_t N>
staticVector<N>
traitsMultiply<Real,staticVector<N> >::
multiply(const Real & A, const staticVector<N> & B)
{
  return(B*A);
}


//_________________________________________________________________________________________________
// REAL-KOMPLEX TRAIT
//-------------------------------------------------------------------------------------------------
/*! Real - Point2d */
template<> class traitsMultiply<Real,komplex>
{
  public:
    typedef komplex DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const Real & A, const komplex & B);
};


//_________________________________________________________________________________________________
// POINT2D-REAL TRAIT
//-------------------------------------------------------------------------------------------------
/*! Point2d - Real */
template<> class traitsMultiply<point2d,Real>
{
  public:
    typedef point2d DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const point2d & A, const Real & B);
};


//_________________________________________________________________________________________________
// POINT2D-POINT2D TRAIT
//-------------------------------------------------------------------------------------------------
/*! Point2d - Point2d */
template<> class traitsMultiply<point2d,point2d>
{
  public:
    typedef tensor2d DIRECTTYPE;
    typedef point3d  CROSSTYPEA;
    typedef Real     SCALARTYPEA;
    
    traitsMultiply() { };
    static DIRECTTYPE  multiply(const point2d & A, const point2d & B);
    static CROSSTYPEA  crossProductA(const point2d & A, const point2d & B);
    static SCALARTYPEA scalarProductA(const point2d & A, const point2d & B);
};


//_________________________________________________________________________________________________
// POINT3D-REAL TRAIT
//-------------------------------------------------------------------------------------------------
/*! Point3d - Real */
template<> class traitsMultiply<point3d,Real>
{
  public:
    typedef point3d DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const point3d & A, const Real & B);
};


//_________________________________________________________________________________________________
// POINT3D-POINT3D TRAIT
//-------------------------------------------------------------------------------------------------
/*! Point3d - Point3d */
template<> class traitsMultiply<point3d,point3d>
{
  public:
    typedef tensor3d DIRECTTYPE;
    typedef point3d  CROSSTYPEA;
    typedef Real     SCALARTYPEA;
    
    traitsMultiply() { };
    static DIRECTTYPE  multiply(const point3d & A, const point3d & B);
    static CROSSTYPEA  crossProductA(const point3d & A, const point3d & B);
    static SCALARTYPEA scalarProductA(const point3d & A, const point3d & B);
};


//_________________________________________________________________________________________________
// TENSOR2D-REAL TRAIT
//-------------------------------------------------------------------------------------------------
/*! Tensor2d - Real */
template<> class traitsMultiply<tensor2d,Real>
{
  public:
    typedef tensor2d DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const tensor2d & A, const Real & B);
};



//_________________________________________________________________________________________________
// TENSOR3D-REAL TRAIT
//-------------------------------------------------------------------------------------------------
/*! Tensor3d - Real */
template<> class traitsMultiply<tensor3d,Real>
{
  public:
    typedef tensor3d DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const tensor3d & A, const Real & B);
};


//_________________________________________________________________________________________________
// STATEVECTOR-REAL TRAIT
//-------------------------------------------------------------------------------------------------
/*! StateVector - Real */
template<> class traitsMultiply<stateVector,Real>
{
  public:
    typedef stateVector DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const stateVector & A, const Real & B);
};


//_________________________________________________________________________________________________
// STATEMATRIX-REAL TRAIT
//-------------------------------------------------------------------------------------------------
/*! StateMatrix - Real */
template<> class traitsMultiply<stateMatrix,Real>
{
  public:
    typedef stateMatrix DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const stateMatrix & A, const Real & B);
};


//_________________________________________________________________________________________________
// STATICVECTOR-REAL TRAIT
//-------------------------------------------------------------------------------------------------
/*! StaticVector - Real */
template<size_t N> class traitsMultiply<staticVector<N>,Real>
{
  public:
    typedef staticVector<N> DIRECTTYPE;
    
    traitsMultiply() { };
    static staticVector<N> multiply(const staticVector<N> & A, const Real & B);
};

template<size_t N>
staticVector<N>
traitsMultiply<staticVector<N>,Real>::
multiply(const staticVector<N> & A, const Real & B)
{
  return(A*B);
}


//_________________________________________________________________________________________________
// KOMPLEX-REAL TRAIT
//-------------------------------------------------------------------------------------------------
/*! Point2d - Real */
template<> class traitsMultiply<komplex,Real>
{
  public:
    typedef komplex DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const komplex & A, const Real & B);
};


//_________________________________________________________________________________________________
// KOMPLEX-KOMPLEX TRAIT
//-------------------------------------------------------------------------------------------------
/*! Point2d - Real */
template<> class traitsMultiply<komplex,komplex>
{
  public:
    typedef komplex DIRECTTYPE;
    
    traitsMultiply() { };
    static DIRECTTYPE multiply(const komplex & A, const komplex & B);
};

#endif
