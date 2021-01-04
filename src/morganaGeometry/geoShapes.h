/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef GEOSHAPES_H
#define GEOSHAPES_H

#include "morganaTypes.hpp"
#include "polyCards.h"
#include "geoRefShapes.h"


template<typename GEOSHAPE, UInt I> struct geoMap;


//_________________________________________________________________________________________________
// LINEAR LINE 
//-------------------------------------------------------------------------------------------------

/*! Secondary geometrical definition - linearLine */
class linearLine: public geoLine
{
  public:
    typedef geoLine   BASREFSHA;
    typedef geoPoint  GEOBSHAPE;
    
  public:
    static const ShapeName geoName   = LINEARLINE;
    static const UInt numPoints      = 2;
    static const UInt nbPtsPerVertex = 1;
    static const UInt nbPtsPerEdge   = 0;
    static const UInt nbPtsPerFace   = 0;
    static const UInt nbPtsPerVolume = 0;
    
  public:
    static point3d getRefNodes(const UInt & i);
    static UInt    getPointsOnEdge(const UInt & i);
    static UInt    getPointsOnFace(const UInt & i);
    static UInt    getEdgesOnFace(const UInt & i);
    static UInt    edgeToPoint(const UInt & localEdge, const UInt & pointId);
    static UInt    faceToPoint(const UInt & localFace, const UInt & pointId);
    static UInt    faceToEdge(const UInt & localFace, const UInt & edgeId);
    
  public:
    static bool isInside(const point3d & P);
    static bool projectInside(const sVect<point3d> & points, const point3d & P, point3d & Q);
    static point3d projection(const point3d & P);
    static point3d getBarycenter();
    static void boundingBox(const sVect<point3d> & points, point3d & Pmin, point3d & Pmax);
};

template<>  struct geoMap<linearLine,1> {typedef P1_1d_A BASE; };
template<>  struct geoMap<linearLine,2> {typedef P1_1d_B BASE; };



//_________________________________________________________________________________________________
// LINEAR TRIANGLE
//-------------------------------------------------------------------------------------------------

/*! Secondary geometrical definition - linearTriangle */
class linearTriangle: public geoTriangle
{
  public:
    typedef geoTriangle  BASREFSHA;
    typedef linearLine   GEOBSHAPE;
    
  public:
    static const ShapeName geoName   = LINEARTRIANGLE;
    static const UInt numPoints      = 3;
    static const UInt nbPtsPerVertex = 1;
    static const UInt nbPtsPerEdge   = 0;
    static const UInt nbPtsPerFace   = 0;
    static const UInt nbPtsPerVolume = 0;
    
  public:
    static point3d getRefNodes(const UInt & i);
    static UInt    getPointsOnEdge(const UInt & i);
    static UInt    getPointsOnFace(const UInt & i);
    static UInt    getEdgesOnFace(const UInt & i);
    static UInt    edgeToPoint(const UInt & localEdge, const UInt & pointId);
    static UInt    faceToPoint(const UInt & localFace, const UInt & pointId);
    static UInt    faceToEdge(const UInt & localFace, const UInt & edgeId);
    
  public:
    static bool isInside(const point3d & P);
    static bool projectInside(const sVect<point3d> & points, const point3d & P, point3d & Q);
    static point3d projection(const point3d & P);
    static point3d getBarycenter();
    static void boundingBox(const sVect<point3d> & points, point3d & Pmin, point3d & Pmax);
};

template<>  struct geoMap<linearTriangle,1> {typedef P1_2d_A BASE; };
template<>  struct geoMap<linearTriangle,2> {typedef P1_2d_B BASE; };
template<>  struct geoMap<linearTriangle,3> {typedef P1_2d_C BASE; };



//_________________________________________________________________________________________________
// LINEAR QUAD
//-------------------------------------------------------------------------------------------------

/*! Secondary geometrical definition - linearQuad */
class linearQuad: public geoQuad
{
  public:
    typedef geoQuad    BASREFSHA;
    typedef linearLine GEOBSHAPE;
    
  public:
    static const ShapeName geoName   = LINEARQUAD;
    static const UInt numPoints      = 4;
    static const UInt nbPtsPerVertex = 1;
    static const UInt nbPtsPerEdge   = 0;
    static const UInt nbPtsPerFace   = 0;
    
  public:
    static point3d getRefNodes(const UInt & i);
    static UInt    getPointsOnEdge(const UInt & i);
    static UInt    getPointsOnFace(const UInt & i);
    static UInt    getEdgesOnFace(const UInt & i);
    static UInt    edgeToPoint(const UInt & localEdge, const UInt & pointId);
    static UInt    faceToPoint(const UInt & localFace, const UInt & pointId);
    static UInt    faceToEdge(const UInt & localFace, const UInt & edgeId);
    
  public:
    static bool isInside(const point3d & P);
    static bool projectInside(const sVect<point3d> & points, const point3d & P, point3d & Q);
    static point3d projection(const point3d & P);
    static point3d getBarycenter();
    static void boundingBox(const sVect<point3d> & points, point3d & Pmin, point3d & Pmax);
};

template<>  struct geoMap<linearQuad,1> {typedef Q1_2d_A BASE; };
template<>  struct geoMap<linearQuad,2> {typedef Q1_2d_B BASE; };
template<>  struct geoMap<linearQuad,3> {typedef Q1_2d_C BASE; };
template<>  struct geoMap<linearQuad,4> {typedef Q1_2d_D BASE; };



//_________________________________________________________________________________________________
// LINEAR TETRA
//-------------------------------------------------------------------------------------------------

/*! Secondary geometrical definition - linearTetra */
class linearTetra: public geoTetra
{
  public:
    typedef geoTetra       BASREFSHA;
    typedef linearTriangle GEOBSHAPE;
    
  public:
    static const ShapeName geoName   = LINEARTETRA;
    static const UInt numPoints      = 4;
    static const UInt nbPtsPerVertex = 1;
    static const UInt nbPtsPerEdge   = 0;
    static const UInt nbPtsPerFace   = 0;
    static const UInt nbPtsPerVolume = 0;
    
  public:
    static point3d getRefNodes(const UInt & i);
    static UInt    getPointsOnEdge(const UInt & i);
    static UInt    getPointsOnFace(const UInt & i);
    static UInt    getEdgesOnFace(const UInt & i);
    static UInt    edgeToPoint(const UInt & localEdge, const UInt & pointId);
    static UInt    faceToPoint(const UInt & localFace, const UInt & pointId);
    static UInt    faceToEdge(const UInt & localFace, const UInt & edgeId);
    static UInt    faceOppositeToNode(const UInt & localNode);
    
  public:
    static bool isInside(const point3d & P);
    static bool projectInside(const sVect<point3d> & points, const point3d & P, point3d & Q);
    static point3d projection(const point3d & P);
    static point3d getBarycenter();
    static void boundingBox(const sVect<point3d> & points, point3d & Pmin, point3d & Pmax);
};

template<>  struct geoMap<linearTetra,1> {typedef P1_3d_A BASE; };
template<>  struct geoMap<linearTetra,2> {typedef P1_3d_B BASE; };
template<>  struct geoMap<linearTetra,3> {typedef P1_3d_C BASE; };
template<>  struct geoMap<linearTetra,4> {typedef P1_3d_D BASE; };



//_________________________________________________________________________________________________
// LINEAR HEXA
//-------------------------------------------------------------------------------------------------

/*! Secondary geometrical definition - linearHexa */
class linearHexa: public geoHexa
{
  public:
    typedef geoHexa    BASREFSHA;
    typedef linearQuad GEOBSHAPE;
    
  public:
    static const ShapeName geoName   = LINEARHEXA;
    static const UInt numPoints      = 8;
    static const UInt nbPtsPerVertex = 1;
    static const UInt nbPtsPerEdge   = 0;
    static const UInt nbPtsPerFace   = 0;
    static const UInt nbPtsPerVolume = 0;
    
  public:
    static point3d getRefNodes(const UInt & i);
    static UInt    getPointsOnEdge(const UInt & i);
    static UInt    getPointsOnFace(const UInt & i);
    static UInt    getEdgesOnFace(const UInt & i);
    static UInt    edgeToPoint(const UInt & localEdge, const UInt & pointId);
    static UInt    faceToPoint(const UInt & localFace, const UInt & pointId);
    static UInt    faceToEdge(const UInt & localFace, const UInt & edgeId);
    
  public:
    static bool isInside(const point3d & P);
    static bool projectInside(const sVect<point3d> & points, const point3d & P, point3d & Q);
    static point3d projection(const point3d & P);
    static point3d getBarycenter();
    static void boundingBox(const sVect<point3d> & points, point3d & Pmin, point3d & Pmax);
};

template<>  struct geoMap<linearHexa,1> {typedef Q1_3d_A BASE; };
template<>  struct geoMap<linearHexa,2> {typedef Q1_3d_B BASE; };
template<>  struct geoMap<linearHexa,3> {typedef Q1_3d_C BASE; };
template<>  struct geoMap<linearHexa,4> {typedef Q1_3d_D BASE; };
template<>  struct geoMap<linearHexa,5> {typedef Q1_3d_E BASE; };
template<>  struct geoMap<linearHexa,6> {typedef Q1_3d_F BASE; };
template<>  struct geoMap<linearHexa,7> {typedef Q1_3d_G BASE; };
template<>  struct geoMap<linearHexa,8> {typedef Q1_3d_H BASE; };

#endif
