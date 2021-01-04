/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef GEOSPECIALSHAPES_H
#define GEOSPECIALSHAPES_H

#include "geoShapes.h"

//_________________________________________________________________________________________________
// LINEAR WEDGE
//-------------------------------------------------------------------------------------------------

/*! Special geometrical definition - linearLine */
class linearWedge
{
  public:
    static const ReferenceShapes   Shape    = NONE;
    static const ReferenceGeometry Geometry = VOLUME;
    static const UInt nDim        = 3;
    static const UInt numFaces    = 5;
    static const UInt numVertices = 6;
    static const UInt numEdges    = 9;
    static const Real refVol;
  
  public:
    static const ShapeName geoName   = LINEARWEDGE;
    static const UInt numPoints      = 6;
    static const UInt nbPtsPerVertex = 1;
    static const UInt nbPtsPerEdge   = 0;
    static const UInt nbPtsPerFace   = 0;
    static const UInt nbPtsPerVolume = 0;
};

template<>  struct geoMap<linearWedge,1> {typedef W1_3d_A BASE; };
template<>  struct geoMap<linearWedge,2> {typedef W1_3d_B BASE; };
template<>  struct geoMap<linearWedge,3> {typedef W1_3d_C BASE; };
template<>  struct geoMap<linearWedge,4> {typedef W1_3d_D BASE; };
template<>  struct geoMap<linearWedge,5> {typedef W1_3d_E BASE; };
template<>  struct geoMap<linearWedge,6> {typedef W1_3d_F BASE; };



#endif
