/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef INTOSCSYMPLEXSUPPORT_H
#define INTOSCSYMPLEXSUPPORT_H

#include "typesInterface.hpp"

#include "geoShapes.h"
#include "pointElement.hpp"

#include "morganaIntegrator.hpp"
#include "morganaFiniteElements.hpp"

class intOscSymplexSupport
{
    /*! @name Typedefs */ //@{
  public:
    typedef pointElement<linearTetra>    TETRA;
    typedef pointElement<linearTriangle> TRIANGLE;
    //@}
    
     /*! @name Angle tolerance */ //@{
  public:
    static const Real ratioToll;
    //@}
  
    /*! @name Static functions */ //@{
  public:
    static komplex intS0_1d(const point3d & K,
                            const point3d & P1,
                            const point3d & P2);
    
    static komplex intS0_2d(const point3d & K,
                            const point3d & P1,
                            const point3d & P2,
                            const point3d & P3);
    
    static komplex intS0_3d(const point3d & K,
                            const point3d & P1,
                            const point3d & P2,
                            const point3d & P3,
                            const point3d & P4);
    
    static sVect<TETRA>       refineTetra(const UInt & nRef, const TETRA    & tetra);
    static sVect<TETRA>      refineTetra2(const UInt & nRef, const TETRA    & tetra);
    static sVect<TRIANGLE> refineTriangle(const UInt & nRef, const TRIANGLE & triangle);
    
    /*template<typename OPERATOR>
    static sVect<komplex> intS1_1d(OPERATOR & op,
                                   const point3d & K,
                                   const point3d & P1,
                                   const point3d & P2,
                                   const point3d & Y1,
                                   const point3d & Y2);*/
    //@}
};

#endif
