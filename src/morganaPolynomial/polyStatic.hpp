/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef POLYSTATIC_HPP
#define POLYSTATIC_HPP

#include "typesInterface.hpp"
#include "polyPow.hpp"

/*! The static evaluation class. Is template with respect to the type of polynomial used.
see \c polyCards.h for the list of polynomials availables */
template<typename POLY>
class polyStatic 
{
  public:
    typedef typename POLY::ROOT ROOT;
    
    /*! @name Constructors */ //@{
  public:
    polyStatic();
    //@}
    
    /*! @name Evaluation functions */ //@{
  public:
    Real evaluate(const point3d & X) const;
    static Real evaluateStatic(const point3d & X);
    //@}
    
    /*! @name Gradient evaluation functions */ //@{
  public:
    static Real evaluateGradientX(const point3d & X);
    static Real evaluateGradientY(const point3d & X);
    static Real evaluateGradientZ(const point3d & X);
    static point3d evaluateGradient(const point3d & X);
    //@}
    
    /*! @name Second derivative evaluation functions */ //@{
  public:
    static Real evaluateHessianXX(const point3d & X);
    static Real evaluateHessianXY(const point3d & X);
    static Real evaluateHessianXZ(const point3d & X);
    static Real evaluateHessianYX(const point3d & X);
    static Real evaluateHessianYY(const point3d & X);
    static Real evaluateHessianYZ(const point3d & X);
    static Real evaluateHessianZX(const point3d & X);
    static Real evaluateHessianZY(const point3d & X);
    static Real evaluateHessianZZ(const point3d & X);
    static tensor3d evaluateHessian(const point3d & X);
    //@}
    
    /*! @name Generic evaluation function */ //@{
  public:
    template<UInt dx, UInt dy, UInt dz>
    static Real evaluateDerivative(const point3d & X);
    //@}
   
    /*! @name Other functions */ //@{
  public:
    static UInt getDegree();
    //@}
    
};



//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename POLY>
polyStatic<POLY>::
polyStatic() 
{ }



//_________________________________________________________________________________________________
// EVALUATION FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename POLY>
Real
polyStatic<POLY>::
evaluate(const point3d & X) const
{
  return(polyEvalIter<ROOT,POLY::S>::eval(X));
}

template<typename POLY>
Real
polyStatic<POLY>::
evaluateStatic(const point3d & X)
{
  return(polyEvalIter<ROOT,POLY::S>::eval(X));
}


//_________________________________________________________________________________________________
// GRADIENT EVALUATION FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename POLY>
Real
polyStatic<POLY>::
evaluateGradientX(const point3d & X)
{
  return(polyDerIter<ROOT,1,0,0,POLY::S>::eval(X));
}

template<typename POLY>
Real
polyStatic<POLY>::
evaluateGradientY(const point3d & X)
{
  return(polyDerIter<ROOT,0,1,0,POLY::S>::eval(X));
}

template<typename POLY>
Real
polyStatic<POLY>::
evaluateGradientZ(const point3d & X)
{
  return(polyDerIter<ROOT,0,0,1,POLY::S>::eval(X));
}

template<typename POLY>
point3d
polyStatic<POLY>::
evaluateGradient(const point3d & X)
{
  return(
  point3d(
  polyDerIter<ROOT,1,0,0,POLY::S>::eval(X),
  polyDerIter<ROOT,0,1,0,POLY::S>::eval(X),
  polyDerIter<ROOT,0,0,1,POLY::S>::eval(X))
  );
}



//_________________________________________________________________________________________________
// HESSIAN EVALUATION FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename POLY>
Real
polyStatic<POLY>::
evaluateHessianXX(const point3d & X)
{
  return(polyDerIter<ROOT,2,0,0,POLY::S>::eval(X));
}

template<typename POLY>
Real
polyStatic<POLY>::
evaluateHessianXY(const point3d & X)
{
  return(polyDerIter<ROOT,1,1,0,POLY::S>::eval(X));
}

template<typename POLY>
Real
polyStatic<POLY>::
evaluateHessianXZ(const point3d & X)
{
  return(polyDerIter<ROOT,1,0,1,POLY::S>::eval(X));
}

template<typename POLY>
Real
polyStatic<POLY>::
evaluateHessianYX(const point3d & X)
{
  return(polyDerIter<ROOT,1,1,0,POLY::S>::eval(X));
}

template<typename POLY>
Real
polyStatic<POLY>::
evaluateHessianYY(const point3d & X)
{
  return(polyDerIter<ROOT,0,2,0,POLY::S>::eval(X));
}

template<typename POLY>
Real
polyStatic<POLY>::
evaluateHessianYZ(const point3d & X)
{
  return(polyDerIter<ROOT,0,1,1,POLY::S>::eval(X));
}

template<typename POLY>
Real
polyStatic<POLY>::
evaluateHessianZX(const point3d & X)
{
  return(polyDerIter<ROOT,1,0,1,POLY::S>::eval(X));
}

template<typename POLY>
Real
polyStatic<POLY>::
evaluateHessianZY(const point3d & X)
{
  return(polyDerIter<ROOT,0,1,1,POLY::S>::eval(X));
}

template<typename POLY>
Real
polyStatic<POLY>::
evaluateHessianZZ(const point3d & X)
{
  return(polyDerIter<ROOT,0,0,2,POLY::S>::eval(X));
}

template<typename POLY>
tensor3d
polyStatic<POLY>::
evaluateHessian(const point3d & X)
{
  return(tensor3d(
  polyDerIter<ROOT,2,0,0,POLY::S>::eval(X), polyDerIter<ROOT,1,1,0,POLY::S>::eval(X), polyDerIter<ROOT,1,0,1,POLY::S>::eval(X),
  polyDerIter<ROOT,1,1,0,POLY::S>::eval(X), polyDerIter<ROOT,0,2,0,POLY::S>::eval(X), polyDerIter<ROOT,0,1,1,POLY::S>::eval(X),
  polyDerIter<ROOT,1,0,1,POLY::S>::eval(X), polyDerIter<ROOT,0,1,1,POLY::S>::eval(X), polyDerIter<ROOT,0,0,2,POLY::S>::eval(X))
  );
}



//_________________________________________________________________________________________________
// DERIVATIVE
//-------------------------------------------------------------------------------------------------
template<typename POLY>
template<UInt dx, UInt dy, UInt dz>
Real
polyStatic<POLY>::
evaluateDerivative(const point3d & X)
{
  return(polyDerIter<ROOT,dx,dy,dz,POLY::S>::eval(X));
}


//_________________________________________________________________________________________________
// PRINTOUT
//-------------------------------------------------------------------------------------------------
template<typename POLY>
UInt
polyStatic<POLY>::
getDegree()
{
  return(polyDegree<ROOT,POLY::S>::maxDegree());
}
    
#endif
