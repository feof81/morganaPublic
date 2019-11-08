/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef POLYDYNAMIC_H
#define POLYDYNAMIC_H

#include "polyDynamicCard.h"


/*! The dynamic evaluation class. The polynomial type is given setting the \c polyDynamicCard wich contains 
all the information that define the polynomial. */
class polyDynamic
{
    /*! @name Internal data */ //@{
  public:
    bool cardLoaded;
    polyDynamicCard card;
    //@}
    
    /*! @name Constructor */ //@{
  public:
    polyDynamic();
    //@}
    
    /*! @name Internal functions */ //@{
  public:
    void setPolyDynamicCard(const polyDynamicCard & Card);
    UInt factorial(const UInt & d, const UInt & p) const;
    //@}
    
    /*! @name Evaluate functions */ //@{
  public:
    Real evaluate(const point3d & X) const;
    //@}
    
    /*! @name Gradient evaluation functions */ //@{
  public:
    Real evaluateGradientX(const point3d & X) const;
    Real evaluateGradientY(const point3d & X) const;
    Real evaluateGradientZ(const point3d & X) const;
    point3d evaluateGradient(const point3d & X) const;
    //@}
    
    /*! @name Second derivative evaluation functions */ //@{
  public:
    Real evaluateHessianXX(const point3d & X) const;
    Real evaluateHessianXY(const point3d & X) const;
    Real evaluateHessianXZ(const point3d & X) const;
    Real evaluateHessianYX(const point3d & X) const;
    Real evaluateHessianYY(const point3d & X) const;
    Real evaluateHessianYZ(const point3d & X) const;
    Real evaluateHessianZX(const point3d & X) const;
    Real evaluateHessianZY(const point3d & X) const;
    Real evaluateHessianZZ(const point3d & X) const;
    tensor3d evaluateHessian(const point3d & X) const;
    //@}
    
    /*! @name Generic evaluation function */ //@{
  public:
    Real evaluateDerivative(const UInt & dx, const UInt & dy, const UInt & dz, const point3d & X) const;
    //@}
    
    /*! @name Other functions */ //@{
  public:
    UInt getDegree() const;
    //@}
};

#endif
