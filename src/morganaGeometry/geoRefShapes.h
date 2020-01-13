/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef GEOREFSHAPES_H
#define GEOREFSHAPES_H

#include "morganaTypes.hpp"
#include "morganaGeometry.hpp"

/*! Primary geometrical definition - point */
class geoPoint
{
  public:
    static const ReferenceShapes   Shape    = POINT; 
    static const ReferenceGeometry Geometry = VERTEX;
    static const UInt nDim        = 0;  
    static const UInt numFaces    = 0; 
    static const UInt numEdges    = 0;
    static const UInt numVertices = 1;
    static const UInt numPoints   = 1;
    static const Real refVol;
};


/*! Primary geometrical definition - line */
class geoLine
{
  public:
    static const ReferenceShapes   Shape    = LINE;
    static const ReferenceGeometry Geometry = EDGE;
    static const UInt nDim        = 1;
    static const UInt numFaces    = 0;
    static const UInt numEdges    = 1;
    static const UInt numVertices = 2;
    static const Real refVol;
};


/*! Primary geometrical definition - triangle */
class geoTriangle
{
  public:
    static const ReferenceShapes   Shape    = TRIANGLE;
    static const ReferenceGeometry Geometry = FACE;
    static const UInt nDim        = 2;
    static const UInt numVertices = 3;
    static const UInt numFaces    = 1;
    static const UInt numEdges    = 3;
    static const Real refVol;
};


/*! Primary geometrical definition - quad */
class geoQuad
{
  public:
    static const ReferenceShapes   Shape    = QUAD;
    static const ReferenceGeometry Geometry = FACE;
    static const UInt nDim        = 2;
    static const UInt numFaces    = 1;
    static const UInt numVertices = 4;
    static const UInt numEdges    = 4;
    static const Real refVol;
};


/*! Primary geometrical definition - tetra */
class geoTetra
{
  public:
    static const ReferenceShapes   Shape    = TETRA;
    static const ReferenceGeometry Geometry = VOLUME;
    static const UInt nDim        = 3;
    static const UInt numVertices = 4;
    static const UInt numFaces    = 4;
    static const UInt numEdges    = 6;
    static const Real refVol;
};


/*! Primary geometrical definition - hexa */
class geoHexa
{
  public:
    static const ReferenceShapes   Shape    = HEXA;
    static const ReferenceGeometry Geometry = VOLUME;
    static const UInt nDim        = 3;
    static const UInt numFaces    = 6;
    static const UInt numVertices = 8;
    static const UInt numEdges    = 12;
    static const Real refVol;
};

#endif
