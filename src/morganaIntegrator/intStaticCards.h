/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef INTSTATICCARDS_H
#define INTSTATICCARDS_H

#include "typesInterface.hpp"
#include "geoShapes.h"
#include "geoSpecialShapes.h"
#include "morganaIntegrator.hpp"



/*! -------------------------------------------------------------------------------------------- */
/*!                                       ENUMS                                                  */
/*! -------------------------------------------------------------------------------------------- */

//! Forward declaration
template<typename GEOSHAPE, intTypes TYPE, UInt PRECISION>
struct intStaticCard;




/*! -------------------------------------------------------------------------------------------- */
/*!                                       INTEGRATION 1D                                         */
/*! -------------------------------------------------------------------------------------------- */

//! LinearLine, STANDARD, Precision 1----------------------------------------------------------------
template<> struct intStaticCard<linearLine,STANDARD,1>
{
  public:
    static const UInt N = 1;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearLine, STANDARD, Precision 2---------------------------------------------------------------.
template<> struct intStaticCard<linearLine,STANDARD,2>
{
  public:
    static const UInt N = 2;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearLine, STANDARD, Precision 3---------------------------------------------------------------
template<> struct intStaticCard<linearLine,STANDARD,3>
{
  public:
    static const UInt N = 3;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearLine, HYBRID, Precision 1----------------------------------------------------------------
template<> struct intStaticCard<linearLine,HYBRID,1>
{
  public:
    static const UInt N = 1;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearLine, HYBRID, Precision 2---------------------------------------------------------------.
template<> struct intStaticCard<linearLine,HYBRID,2>
{
  public:
    static const UInt N = 2;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearLine, HYBRID, Precision 3---------------------------------------------------------------
template<> struct intStaticCard<linearLine,HYBRID,3>
{
  public:
    static const UInt N = 3;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};



/*! -------------------------------------------------------------------------------------------- */
/*!                                       INTEGRATION 2D                                         */
/*! -------------------------------------------------------------------------------------------- */

//! LinearTriangle, STANDARD, Precision 1------------------------------------------------------------
template<> struct intStaticCard<linearTriangle,STANDARD,1>
{
  public:
    static const UInt N = 1;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearTriangle, STANDARD, Precision 2------------------------------------------------------------
template<> struct intStaticCard<linearTriangle,STANDARD,2>
{
  public:
    static const UInt N = 3;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearTriangle, STANDARD, Precision 3------------------------------------------------------------
template<> struct intStaticCard<linearTriangle,STANDARD,3>
{
  public:
    static const UInt N = 4;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearTriangle, P0DUAL, Precision 4------------------------------------------------------------
template<> struct intStaticCard<linearTriangle,STANDARD,4>
{
  public:
    static const UInt N = 6;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearTriangle, P0DUAL, Precision 1------------------------------------------------------------
template<> struct intStaticCard<linearTriangle,MINIELEMENT,1>
{
  public:
    static const UInt N = 3;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearTriangle, P0DUAL, Precision 1------------------------------------------------------------
template<> struct intStaticCard<linearTriangle,P0DUAL,1>
{
  public:
    static const UInt N = 3;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearTriangle, HYBRID, Precision 1------------------------------------------------------------
template<> struct intStaticCard<linearTriangle,HYBRID,1>
{
  public:
    static const UInt N = 1;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearTriangle, HYBRID, Precision 2------------------------------------------------------------
template<> struct intStaticCard<linearTriangle,HYBRID,2>
{
  public:
    static const UInt N = 4;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearTriangle, HYBRID, Precision 3------------------------------------------------------------
template<> struct intStaticCard<linearTriangle,HYBRID,3>
{
  public:
    static const UInt N = 6;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearQuad, STANDARD, Precision 1---------------------------------------------------------------
template<> struct intStaticCard<linearQuad,STANDARD,1>
{
  public:
    static const UInt N = 9;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearQuad, STANDARD, Precision 2---------------------------------------------------------------
template<> struct intStaticCard<linearQuad,STANDARD,2>
{
  public:
    static const UInt N = 16;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearQuad, STANDARD, Precision 3---------------------------------------------------------------
template<> struct intStaticCard<linearQuad,STANDARD,3>
{
  public:
    static const UInt N = 25;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearQuad, HYBRID, Precision 1---------------------------------------------------------------
template<> struct intStaticCard<linearQuad,HYBRID,1>
{
  public:
    static const UInt N = 1;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearQuad, HYBRID, Precision 2---------------------------------------------------------------
template<> struct intStaticCard<linearQuad,HYBRID,2>
{
  public:
    static const UInt N = 4;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearQuad, HYBRID, Precision 3---------------------------------------------------------------
template<> struct intStaticCard<linearQuad,HYBRID,3>
{
  public:
    static const UInt N = 9;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


/*! -------------------------------------------------------------------------------------------- */
/*!                                       INTEGRATION 3D                                         */
/*! -------------------------------------------------------------------------------------------- */

//! LinearTetra, STANDARD, Precision 1---------------------------------------------------------------
template<> struct intStaticCard<linearTetra,STANDARD,1>
{
  public:
    static const UInt N = 1;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearTetra, STANDARD, Precision 2---------------------------------------------------------------
template<> struct intStaticCard<linearTetra,STANDARD,2>
{
  public:
    static const UInt N = 4;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearTetra, STANDARD, Precision 3---------------------------------------------------------------
template<> struct intStaticCard<linearTetra,STANDARD,3>
{
  public:
    static const UInt N = 15;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearTetra, STANDARD, Precision 4---------------------------------------------------------------
template<> struct intStaticCard<linearTetra,STANDARD,4>
{
  public:
    static const UInt N = 64;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearTetra, MINIELEMENT, Precision 1---------------------------------------------------------------
template<> struct intStaticCard<linearTetra,MINIELEMENT,1>
{
  public:
    static const UInt N = 4;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearTetra, P0DUAL, Precision 1---------------------------------------------------------------
template<> struct intStaticCard<linearTetra,P0DUAL,1>
{
  public:
    static const UInt N = 4;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearHexa, STANDARD, Precision 1---------------------------------------------------------------
template<> struct intStaticCard<linearHexa,STANDARD,1>
{
  public:
    static const UInt N = 27;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearQuad, STANDARD, Precision 2---------------------------------------------------------------
template<> struct intStaticCard<linearHexa,STANDARD,2>
{
  public:
    static const UInt N = 64;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearQuad, STANDARD, Precision 3---------------------------------------------------------------
template<> struct intStaticCard<linearHexa,STANDARD,3>
{
  public:
    static const UInt N = 125;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearWedge, STANDARD, Precision 1---------------------------------------------------------------
template<> struct intStaticCard<linearWedge,STANDARD,1>
{
  public:
    static const UInt N = 1;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};


//! LinearWedge, STANDARD, Precision 2---------------------------------------------------------------
template<> struct intStaticCard<linearWedge,STANDARD,2>
{
  public:
    static const UInt N = 6;
    static point3d getYn(const UInt & i);
    static Real    getWn(const UInt & i);
};





#endif
