/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "geoShapes.h"
#include "../morganaDofs/morganaTypes.hpp"

using namespace std;

//_________________________________________________________________________________________________
// LINEAR LINE 
//-------------------------------------------------------------------------------------------------
point3d
linearLine::
getRefNodes(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 2);
  
  //Data definition
  static const point3d refNodes[2] =  { 
    point3d(0.0, 0.0, 0.0),
    point3d(1.0, 0.0, 0.0) };
    
  return(refNodes[i-1]);
}

UInt
linearLine::
getPointsOnEdge(const UInt & i)
{
  assert(i==i);
  
  static const UInt p2e = 2;
  return(p2e);
}

UInt
linearLine::
getPointsOnFace(const UInt & i)
{
  assert(i==i);
  
  static const UInt p2f = 0;
  return(p2f);
}

UInt
linearLine::
getEdgesOnFace(const UInt & i)
{
  assert(i==i);
  
  static const UInt e2f = 0;
  return(e2f);
}

UInt
linearLine::
edgeToPoint(const UInt & localEdge, const UInt & pointId)
{
  assert(1==2);
  assert(localEdge == localEdge);
  assert(pointId   == pointId);
  
  return(0);
}

UInt
linearLine::
faceToPoint(const UInt & localFace, const UInt & pointId)
{
  assert(1==2);
  assert(localFace == localFace);
  assert(pointId   == pointId);
  
  return(0);
}

UInt
linearLine::
faceToEdge(const UInt & localFace, const UInt & edgeId)
{
  assert(1==2);
  assert(localFace == localFace);
  assert(edgeId    == edgeId);
  
  return(0);
}

bool
linearLine::
isInside(const point3d & P)
{
  return(
     ( (P.getX() >= -geoToll) && (P.getX() <= (1.0 + geoToll)) )
  && ( (P.getY() >= -geoToll) && (P.getY() <= geoToll)         )
  && ( (P.getZ() >= -geoToll) && (P.getZ() <= geoToll)         )  
  );
}

bool
linearLine::
projectInside(const sVect<point3d> & points, const point3d & P, point3d & Q)
{
  assert(points.size() == 2);
  
  point3d D = points(2) - points(1);
  Real  out = max(0.0,min(1.0, (point3d::dot(P - points(1), D)) / point3d::dot(D,D)) );
          Q = points(1) * (1.0-out) + points(2) * out;
  
  return(true);
}

point3d
linearLine::
projection(const point3d & P)
{
  point3d Y;
  
  Y.setX( max(0.0, P.getX()) ); Y.setX( min(1.0, Y.getX()) );
  
  return(Y);
}

point3d
linearLine::
getBarycenter()
{
  static point3d Pb = point3d(0.5,0.0,0.);
  return(Pb);
}

void
linearLine::
boundingBox(const sVect<point3d> & points, point3d & Pmin, point3d & Pmax)
{
  assert(points.size() == 2);
  
  Pmin = points(1);
  Pmax = points(1);
  
  for(UInt i=1; i <= 2; ++i)
  {
    Pmin.setX( std::min(Pmin.getX(), points(i).getX()) );
    Pmin.setY( std::min(Pmin.getY(), points(i).getY()) );
    Pmin.setZ( std::min(Pmin.getZ(), points(i).getZ()) );
    
    Pmax.setX( std::max(Pmax.getX(), points(i).getX()) );
    Pmax.setY( std::max(Pmax.getY(), points(i).getY()) );
    Pmax.setZ( std::max(Pmax.getZ(), points(i).getZ()) );
  }
}



//_________________________________________________________________________________________________
// LINEAR TRIANGLE
//-------------------------------------------------------------------------------------------------
point3d
linearTriangle::
getRefNodes(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 3);
  
  //Data definition
  static const point3d refNodes[3] =  { 
    point3d(0.0, 0.0, 0.0),
    point3d(1.0, 0.0, 0.0),
    point3d(0.0, 1.0, 0.0)};
    
  return(refNodes[i-1]);
}

UInt
linearTriangle::
getPointsOnEdge(const UInt & i)
{
  assert(i==i);
  
  static const UInt p2e = 2;
  return(p2e);
}

UInt
linearTriangle::
getPointsOnFace(const UInt & i)
{
  assert(i==i);
  
  static const UInt p2f = 0;
  return(p2f);
}

UInt
linearTriangle::
getEdgesOnFace(const UInt & i)
{
  assert(i==i);
  
  static const UInt e2f = 0;
  return(e2f);
}

UInt
linearTriangle::
edgeToPoint(const UInt & localEdge, const UInt & pointId)
{
  static const UInt eToP[ 2 * numEdges ] =  { 1, 2, 2, 3, 3, 1 };
    
  assert( pointId > 0 && pointId < 3 ) ;
  assert( localEdge > 0 && localEdge <= numEdges ) ;
  return eToP[ 2 * localEdge + pointId - 3 ];
}

UInt
linearTriangle::
faceToPoint(const UInt & localFace, const UInt & pointId)
{
  assert(1==2);
  assert(localFace == localFace);
  assert(pointId   == pointId);
  
  return(0);
}

UInt
linearTriangle::
faceToEdge(const UInt & localFace, const UInt & edgeId)
{
  assert(1==2);
  assert(localFace == localFace);
  assert(edgeId    == edgeId);
  
  return(0);
}

bool
linearTriangle::
isInside(const point3d & P)
{
  return(
  (P.getX() >= -geoToll) &&
  (P.getY() >= -geoToll) &&
  ( (P.getX() + P.getY()) <= (1.0 + geoToll) ) &&
  (P.getZ() >= -geoToll) && (P.getZ() <= geoToll) );
}

bool
linearTriangle::
projectInside(const sVect<point3d> & points, const point3d & P, point3d & Q)
{
  /*assert(points.size() == 3);
  
  //Servo-alloc
  point3d P1 = points(1);
  point3d P2 = points(2);
  point3d P3 = points(3);
  
  point3d Q1 = P - P1;
  point3d Q2 = P - P2;
  point3d Q3 = P - P3;
  
  point3d P12 = P2 - P1;
  point3d P23 = P3 - P2;
  point3d P31 = P1 - P3;
  
  //Volume
  Real a   =   point3d::dot(P12,P12);
  Real b   = - point3d::dot(P12,P31);
  Real c   =   point3d::dot(P31,P31);
  Real t1  =   point3d::dot(Q1,P12);
  Real t2  = - point3d::dot(Q1,P31);
  Real det = a * c - b*b;
  Real csi = (c * t1 - b * t2) / det;
  Real eta = (a * t2 - b * t1) / det;
  
  point3d Pv = P1 + (P12 * csi) - (P31 * eta); 
  bool volume = (csi >= 0.0) && (eta >= 0.0) && ((csi + eta) <= 1.0);
  
  //Nodes
  Real l12 = P12.norm2();
  Real l23 = P23.norm2();
  Real l31 = P31.norm2();
  
  P12 /= l12;
  P23 /= l23;
  P31 /= l31;
  
  bool node1 = (point3d::dot(Q1, P31 * (-1.0)) <= 0.0) && (point3d::dot(Q1, P12) <= 0.0);
  bool node2 = (point3d::dot(Q2, P12 * (-1.0)) <= 0.0) && (point3d::dot(Q2, P23) <= 0.0);
  bool node3 = (point3d::dot(Q3, P23 * (-1.0)) <= 0.0) && (point3d::dot(Q3, P31) <= 0.0);
  bool node  = (node1 || node2 || node3);
  assert( (UInt(node1) + UInt(node2) + UInt(node3)) <= 1 );
  
  //Faces
  point3d D1 = Q1 - P12 * std::min( abs(point3d::dot(Q1, P12)) , l12);
  point3d D2 = Q2 - P23 * std::min( abs(point3d::dot(Q2, P23)) , l23);
  point3d D3 = Q3 - P31 * std::min( abs(point3d::dot(Q3, P31)) , l31);

  Real minDist = std::min(D1.norm2(), std::min(D2.norm2(), D3.norm2()));
  bool face1   = (minDist == D1.norm2());
  bool face3   = (minDist == D2.norm2());
  bool face2   = (minDist == D3.norm2());
  bool face    = (!node) && (!volume);
  
  D1 = P1 + P12 * std::min( abs(point3d::dot(Q1, P12)) , l12);
  D2 = P2 + P23 * std::min( abs(point3d::dot(Q2, P23)) , l23);
  D3 = P3 + P31 * std::min( abs(point3d::dot(Q3, P31)) , l31);
  
  Q = ( P1 * Real(node1) + P2 * Real(node2) + P3 * Real(node3) ) * Real(node) +
      ( D1 * Real(face1) + D2 * Real(face2) + D3 * Real(face3) ) * Real(face) +
        Pv * Real(volume);
  
  return(true);*/
  
  
  //Servo-alloc
  point3d P1 = points(1);
  point3d P2 = points(2);
  point3d P3 = points(3);
  
  point3d Q1 = P - P1;
  point3d Q2 = P - P2;
  point3d Q3 = P - P3;
  
  point3d P12 = P2 - P1;
  point3d P23 = P3 - P2;
  point3d P31 = P1 - P3;
  
  //Volume
  Real a   =   point3d::dot(P12,P12);
  Real b   = - point3d::dot(P12,P31);
  Real c   =   point3d::dot(P31,P31);
  Real t1  =   point3d::dot(Q1,P12);
  Real t2  = - point3d::dot(Q1,P31);
  Real det = a * c - b*b;
  Real csi = (c * t1 - b * t2) / det;
  Real eta = (a * t2 - b * t1) / det;
  
  point3d Pv = P1 + (P12 * csi) - (P31 * eta); 
  bool volume = (csi >= 0.0) && (eta >= 0.0) && ((csi + eta) <= 1.0);
  
  //Faces
  sVect<point3d> facePoints(2);
  point3d D1, D2, D3;

  facePoints(1) = P1; facePoints(2) = P2; linearLine::projectInside(facePoints,P,D1);
  facePoints(1) = P2; facePoints(2) = P3; linearLine::projectInside(facePoints,P,D2);
  facePoints(1) = P3; facePoints(2) = P1; linearLine::projectInside(facePoints,P,D3);

  Real d1 = point3d::norm2(P-D1);
  Real d2 = point3d::norm2(P-D2);
  Real d3 = point3d::norm2(P-D3);

  Real minDist = std::min(d1,std::min(d2,d3));
  bool face1   = (minDist == d1);
  bool face2   = (minDist == d2);
  bool face3   = (minDist == d3);
  bool face    = (!volume);

  if( (UInt(face1) + UInt(face2) + UInt(face3)) >=1 )
  {
    if(face1) { face2 = false; face3 = false;}
    if(face2) { face3 = false;} 
  }

  assert( (UInt(face1) + UInt(face2) + UInt(face3)) == 1);
  
  Q = ( D1 * Real(face1) + D2 * Real(face2) + D3 * Real(face3) ) * Real(face) +
        Pv * Real(volume);
  
  return(true);
  
}

point3d
linearTriangle::
projection(const point3d & P)
{  
  //Plane projection
  point3d Y = P;
  Y.setZ(0.0);
  
  if( isInside(Y) )
  { return(Y); }
  else
  {
    point3d Y1, Y2, Y3;
    sVect<point3d> points(2);
    
    points(1) = point3d(1.0, 0.0, 0.0); points(2) = point3d(0.0, 1.0, 0.0); linearLine::projectInside(points,Y,Y1);
    points(1) = point3d(0.0, 0.0, 0.0), points(2) = point3d(0.0, 1.0, 0.0); linearLine::projectInside(points,Y,Y2);
    points(1) = point3d(0.0, 0.0, 0.0), points(2) = point3d(1.0, 0.0, 0.0); linearLine::projectInside(points,Y,Y3);
    
    Real dist1 = point3d::norm2(Y1-Y);
    Real dist2 = point3d::norm2(Y2-Y);
    Real dist3 = point3d::norm2(Y3-Y);
    
    if( (dist1 <= dist2) && (dist1 <= dist3) )
    { return(Y1); }
    
    if( (dist2 <= dist1) && (dist2 <= dist3) )
    { return(Y2); }
    
    if( (dist3 <= dist2) && (dist3 <= dist1) )
    { return(Y3); }
  }
  
  assert(2==1);
  return(P);
}

point3d
linearTriangle::
getBarycenter()
{
  static point3d Pb = point3d(1.0/3.0 , 1.0/3.0 ,0.);
  return(Pb);
}

void
linearTriangle::
boundingBox(const sVect<point3d> & points, point3d & Pmin, point3d & Pmax)
{
  assert(points.size() == 3);
  
  Pmin = points(1);
  Pmax = points(1);
  
  for(UInt i=1; i <= 3; ++i)
  {
    Pmin.setX( std::min(Pmin.getX(), points(i).getX()) );
    Pmin.setY( std::min(Pmin.getY(), points(i).getY()) );
    Pmin.setZ( std::min(Pmin.getZ(), points(i).getZ()) );
    
    Pmax.setX( std::max(Pmax.getX(), points(i).getX()) );
    Pmax.setY( std::max(Pmax.getY(), points(i).getY()) );
    Pmax.setZ( std::max(Pmax.getZ(), points(i).getZ()) );
  }
}



//_________________________________________________________________________________________________
// LINEAR QUAD
//-------------------------------------------------------------------------------------------------
point3d
linearQuad::
getRefNodes(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 4);
  
  //Data definition
  static const point3d refNodes[4] =  { 
    point3d(0.0, 0.0, 0.0),
    point3d(1.0, 0.0, 0.0),
    point3d(1.0, 1.0, 0.0),
    point3d(0.0, 1.0, 0.0)};
    
  return(refNodes[i-1]);
}

UInt
linearQuad::
getPointsOnEdge(const UInt & i)
{
  assert(i==i);
  
  static const UInt p2e = 2;
  return(p2e);
}

UInt
linearQuad::
getPointsOnFace(const UInt & i)
{
  assert(i==i);
  
  static const UInt p2f = 0;
  return(p2f);
}

UInt
linearQuad::
getEdgesOnFace(const UInt & i)
{
  assert(i==i);
  
  static const UInt e2f = 0;
  return(e2f);
}

UInt
linearQuad::
edgeToPoint(const UInt & localEdge, const UInt & pointId)
{
  static const UInt eToP[ 2 * numEdges ] = { 1, 2, 2, 3, 3, 4, 4, 1 };
  
  assert( pointId > 0 && pointId < 3 ) ;
  assert( localEdge > 0 && localEdge <= numEdges ) ;
  return eToP[ 2 * localEdge + pointId - 3 ];
}

UInt
linearQuad::
faceToPoint(const UInt & localFace, const UInt & pointId)
{
  assert(1==2);
  assert(localFace == localFace);
  assert(pointId   == pointId);
  
  return(0);
}

UInt
linearQuad::
faceToEdge(const UInt & localFace, const UInt & edgeId)
{
  assert(1==2);
  assert(localFace == localFace);
  assert(edgeId    == edgeId);
  
  return(0);
}
    
bool
linearQuad::
isInside(const point3d & P)
{
  return(
  (P.getX() >= -geoToll)        &&
  (P.getY() >= -geoToll)        &&
  (P.getX() <= (1.0 + geoToll)) &&
  (P.getY() <= (1.0 + geoToll))  );
}

bool
linearQuad::
projectInside(const sVect<point3d> & points, const point3d & P, point3d & Q)
{
  assert(points.size() == 4);
  
  point3d Q1, Q2;
  sVect<point3d> subPoints(3);
  
  subPoints(1) = points(1);
  subPoints(2) = points(2);
  subPoints(3) = points(3);
  linearTriangle::projectInside(subPoints,P,Q1);
  
  subPoints(1) = points(1);
  subPoints(2) = points(3);
  subPoints(3) = points(4);
  linearTriangle::projectInside(subPoints,P,Q2);
  
  Real d1 = point3d::norm2(P-Q1);
  Real d2 = point3d::norm2(P-Q2);
  
  Q = Q1 * (d1 <= d2) + Q2 * (d2 < d1); 
  
  return(true);
}

point3d
linearQuad::
projection(const point3d & P)
{
  point3d Y;
  
  Y.setX( max(0.0, P.getX()) ); Y.setX( min(1.0, Y.getX()) );
  Y.setY( max(0.0, P.getY()) ); Y.setY( min(1.0, Y.getY()) );
  
  return(Y);
}

point3d
linearQuad::
getBarycenter()
{
  return(point3d(0.5,0.5,0.0));
}

void
linearQuad::
boundingBox(const sVect<point3d> & points, point3d & Pmin, point3d & Pmax)
{
  assert(points.size() == 4);
  
  Pmin = points(1);
  Pmax = points(1);
  
  for(UInt i=1; i <= 4; ++i)
  {
    Pmin.setX( std::min(Pmin.getX(), points(i).getX()) );
    Pmin.setY( std::min(Pmin.getY(), points(i).getY()) );
    Pmin.setZ( std::min(Pmin.getZ(), points(i).getZ()) );
    
    Pmax.setX( std::max(Pmax.getX(), points(i).getX()) );
    Pmax.setY( std::max(Pmax.getY(), points(i).getY()) );
    Pmax.setZ( std::max(Pmax.getZ(), points(i).getZ()) );
  }
}



//_________________________________________________________________________________________________
// LINEAR TETRA
//-------------------------------------------------------------------------------------------------
point3d
linearTetra::
getRefNodes(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 4);
  
  //Data definition
  static const point3d refNodes[4] =  {
  point3d(0.0, 0.0, 0.0),
  point3d(1.0, 0.0, 0.0),
  point3d(0.0, 1.0, 0.0),
  point3d(0.0, 0.0, 1.0) };
  
  return(refNodes[i-1]);
}

UInt
linearTetra::
getPointsOnEdge(const UInt & i)
{
  assert(i==i);
  
  static const UInt p2e = 2;
  return(p2e);
}

UInt
linearTetra::
getPointsOnFace(const UInt & i)
{
  assert(i==i);
  
  static const UInt p2f = 3;
  return(p2f);
}

UInt
linearTetra::
getEdgesOnFace(const UInt & i)
{
  assert(i==i);
  
  static const UInt e2f = 3;
  return(e2f);
}

UInt
linearTetra::
edgeToPoint(const UInt & localEdge, const UInt & pointId)
{
  static const UInt eToP[ 2 * numEdges ] = { 1, 2, 2, 3, 3, 1, 1, 4, 2, 4, 3, 4 };
    
  assert( pointId > 0 && pointId < 3 ) ;
  assert( localEdge > 0 && localEdge <= numEdges ) ;
  return eToP[ 2 * localEdge + pointId - 3 ];
}

UInt
linearTetra::
faceToPoint(const UInt & localFace, const UInt & pointId)
{
  static const UInt fToP[ 3 * numFaces ] = { 1, 3, 2, 1, 2, 4, 2, 3, 4, 1, 4, 3 }; 

  assert( pointId > 0 && pointId < 4 ) ;
  assert( localFace > 0 && localFace <= numFaces ) ;
  return fToP[ 3 * localFace + pointId - 4 ];
}

UInt
linearTetra::
faceToEdge(const UInt & localFace, const UInt & edgeId)
{
  static const UInt fToE[ 3 * numFaces ] = { 3, 2, 1, 1, 5, 4, 2, 6, 5, 4, 6, 3 };
  
  assert( edgeId > 0 && edgeId < 4 ) ;
  assert( localFace > 0 && localFace <= numFaces ) ;
  return fToE[ 3 * localFace + edgeId - 4 ];
}

UInt
linearTetra::
faceOppositeToNode(const UInt & localNode)
{
  assert(localNode >= 1);
  assert(localNode <= 4);
  
  static const UInt fOpN[ 4 ] = {3, 4, 2, 1};
  
  return( fOpN[ localNode-1 ] );
}

bool
linearTetra::
isInside(const point3d & P)
{
  return(
  (P.getX() >= -geoToll) &&
  (P.getY() >= -geoToll) &&
  (P.getZ() >= -geoToll) &&
  ( (P.getX() + P.getY() + P.getZ()) <= (1.0 + geoToll) )
  );
}

bool
linearTetra::
projectInside(const sVect<point3d> & points, const point3d & P, point3d & Q)
{
  assert(points.size() == 4);
  
  //Download data
  point3d P1 = points(1);
  point3d P2 = points(2);
  point3d P3 = points(3);
  point3d P4 = points(4);
  
  //Volume prj
  point3d P12 = P2 - P1;
  point3d P13 = P3 - P1;
  point3d P14 = P4 - P1;
  
  Real a11 = point3d::dot(P12,P12);
  Real a12 = point3d::dot(P12,P13);
  Real a13 = point3d::dot(P12,P14);
  Real a22 = point3d::dot(P13,P13);
  Real a23 = point3d::dot(P13,P14);
  Real a33 = point3d::dot(P14,P14);
  
  Real t1 = point3d::dot(P-P1, P12);
  Real t2 = point3d::dot(P-P1, P13);
  Real t3 = point3d::dot(P-P1, P14);
  
  Real det = a33 * a12 * a12 -
       2.0 * a12 * a13 * a23 +
             a22 * a13 * a13 +
             a11 * a23 * a23 -
             a11 * a22 * a33;

  Real b11 = a23 * a23 - a22 * a33;
  Real b12 = a12 * a33 - a13 * a23;
  Real b13 = a13 * a22 - a12 * a23;
  Real b22 = a13 * a13 - a11 * a33;
  Real b23 = a11 * a23 - a12 * a13;
  Real b33 = a12 * a12 - a11 * a22;
  
  Real x1 = (b11 * t1 + b12 * t2 + b13 * t3) / det;
  Real x2 = (b12 * t1 + b22 * t2 + b23 * t3) / det;
  Real x3 = (b13 * t1 + b23 * t2 + b33 * t3) / det;
  
  bool volume = (x1 >= 0.0) && (x2 >= 0.0) && (x3 >= 0.0) && ((x1+x2+x3) <= 1.0);
  
  //Face projection
  point3d Q1, Q2, Q3, Q4;  
  sVect<point3d> points2d(3);
  
  points2d(1) = P2;
  points2d(2) = P3;
  points2d(3) = P4;
  linearTriangle::projectInside(points2d,P,Q1);
  
  points2d(1) = P1;
  points2d(2) = P4;
  points2d(3) = P3;
  linearTriangle::projectInside(points2d,P,Q2);
  
  points2d(1) = P1;
  points2d(2) = P2;
  points2d(3) = P4;
  linearTriangle::projectInside(points2d,P,Q3);
  
  points2d(1) = P1;
  points2d(2) = P3;
  points2d(3) = P2;
  linearTriangle::projectInside(points2d,P,Q4);
  
  Real  d1 = point3d::norm2(P-Q1);
  Real  d2 = point3d::norm2(P-Q2);
  Real  d3 = point3d::norm2(P-Q3);
  Real  d4 = point3d::norm2(P-Q4);
  Real   d = std::min(std::min(d1,d2), std::min(d3,d4));
  
  bool flag1 = (d1 == d);
  bool flag2 = (d2 == d);
  bool flag3 = (d3 == d);
  bool flag4 = (d4 == d);
  Real num = Real(flag1) + Real(flag2) + Real(flag3) + Real(flag4);
  assert((num != 0) || (volume));
  
  if(num != 1)
  {    
    if(flag1) { flag2 = false; flag3 = false; flag4 = false; }
    if(flag2) { flag3 = false; flag4 = false; }
    if(flag3) { flag4 = false; }
  }

  Q = P * Real(volume) +
    (Q1 * Real(flag1) + Q2 * Real(flag2) + Q3 * Real(flag3) + Q4 * Real(flag4)) * Real(!volume);
  
  return(true);
}

point3d
linearTetra::
projection(const point3d & P)
{ 
  if( isInside(P) )
  { 
    return(P);
  }
  else
  {
    point3d Y1, Y2, Y3, Y4;
    sVect<point3d> points(3);
    
    points(1) = point3d(0.0, 0.0, 0.0); points(2) = point3d(0.0, 1.0, 0.0); points(3) = point3d(1.0, 0.0, 0.0); linearTriangle::projectInside(points,P,Y1);
    points(1) = point3d(0.0, 0.0, 0.0); points(2) = point3d(1.0, 0.0, 0.0); points(3) = point3d(0.0, 0.0, 1.0); linearTriangle::projectInside(points,P,Y2);
    points(1) = point3d(1.0, 0.0, 0.0); points(2) = point3d(0.0, 1.0, 0.0); points(3) = point3d(0.0, 0.0, 1.0); linearTriangle::projectInside(points,P,Y3);
    points(1) = point3d(0.0, 0.0, 0.0); points(2) = point3d(0.0, 0.0, 1.0); points(3) = point3d(0.0, 1.0, 0.0); linearTriangle::projectInside(points,P,Y4);
    
    Real dist1 = point3d::norm2(Y1 - P);
    Real dist2 = point3d::norm2(Y2 - P);
    Real dist3 = point3d::norm2(Y3 - P);
    Real dist4 = point3d::norm2(Y4 - P);
    
    Real minDist = min(min(dist1, dist2) , min(dist3 , dist4));
    
    if(dist1 <= minDist + geoToll)
    { return(Y1); }
    
    if(dist2 <= minDist + geoToll)
    { return(Y2); }
    
    if(dist3 <= minDist + geoToll)
    { return(Y3); }
    
    if(dist4 <= minDist + geoToll)
    { return(Y4); }
  }
  
  assert(2==1);
  return(P);
}

point3d
linearTetra::
getBarycenter()
{
  static point3d Pb = point3d(1.0/4.0 , 1.0/4.0, 1.0/4.0);
  return(Pb);
}

void
linearTetra::
boundingBox(const sVect<point3d> & points, point3d & Pmin, point3d & Pmax)
{
  assert(points.size() == 4);
  
  Pmin = points(1);
  Pmax = points(1);
  
  for(UInt i=1; i <= 4; ++i)
  {
    Pmin.setX( std::min(Pmin.getX(), points(i).getX()) );
    Pmin.setY( std::min(Pmin.getY(), points(i).getY()) );
    Pmin.setZ( std::min(Pmin.getZ(), points(i).getZ()) );
    
    Pmax.setX( std::max(Pmax.getX(), points(i).getX()) );
    Pmax.setY( std::max(Pmax.getY(), points(i).getY()) );
    Pmax.setZ( std::max(Pmax.getZ(), points(i).getZ()) );
  }
}



//_________________________________________________________________________________________________
// LINEAR HEXA
//-------------------------------------------------------------------------------------------------
point3d
linearHexa::
getRefNodes(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 8);
  
  //Data definition
  static const point3d refNodes[8] =  {
    point3d(0.0, 0.0, 0.0),
    point3d(1.0, 0.0, 0.0),
    point3d(1.0, 1.0, 0.0),
    point3d(0.0, 1.0, 0.0),
    point3d(0.0, 0.0, 1.0),
    point3d(1.0, 0.0, 1.0),
    point3d(1.0, 1.0, 1.0),
    point3d(0.0, 1.0, 1.0) };
  
  return(refNodes[i-1]);
}

UInt
linearHexa::
getPointsOnEdge(const UInt & i)
{
  assert(i==i);
  
  static const UInt p2f = 3;
  return(p2f);
}

UInt
linearHexa::
getPointsOnFace(const UInt & i)
{
  assert(i==i);
  
  static const UInt p2f = 4;
  return(p2f);
}

UInt
linearHexa::
getEdgesOnFace(const UInt & i)
{
  assert(i==i);
  
  static const UInt e2f = 4;
  return(e2f);
}

UInt
linearHexa::
edgeToPoint(const UInt & localEdge, const UInt & pointId)
{
  static const UInt eToP[ 2 * numEdges ] = {
            1, 2, 2, 3, 3, 4, 4, 1,
            1, 5, 2, 6, 3, 7, 4, 8,
            5, 6, 6, 7, 7, 8, 8, 5 };
	
    assert( pointId > 0 && pointId < 3 ) ;
    assert( localEdge > 0 && localEdge <= numEdges ) ;
    return eToP[ 2 * localEdge + pointId - 3 ];
}

UInt
linearHexa::
faceToPoint(const UInt & localFace, const UInt & pointId)
{
  static const UInt fToP[ 4 * numFaces ] = {
            1, 4, 3, 2,
            1, 5, 8, 4,
            1, 2, 6, 5,
            2, 3, 7, 6,
            3, 4, 8, 7,
            5, 6, 7, 8 };
	    
    assert( pointId > 0 && pointId < 5 ) ;
    assert( localFace > 0 && localFace <= numFaces ) ;
    return fToP[ 4 * localFace + pointId - 5 ];
}

UInt
linearHexa::
faceToEdge(const UInt & localFace, const UInt & edgeId)
{
  static const UInt fToE[ 4 * numFaces ] = {
            4, 3, 2, 1, 5, 12, 8, 4, 1, 6, 9, 5,
            2, 7, 10, 6, 3, 8, 11, 7, 9, 10, 11, 12 };

  assert( edgeId > 0 && edgeId < 5 ) ;
  assert( localFace > 0 && localFace <= numFaces ) ;
  return  fToE[ 4 * localFace + edgeId - 5 ];
}
    
bool
linearHexa::
isInside(const point3d & P)
{
  return(
  (P.getX() >= -geoToll)        &&
  (P.getY() >= -geoToll)        &&
  (P.getZ() >= -geoToll)        &&
  (P.getX() <= (1.0 + geoToll)) &&
  (P.getY() <= (1.0 + geoToll)) &&
  (P.getZ() <= (1.0 + geoToll))  
  );
}

bool
linearHexa::
projectInside(const sVect<point3d> & points, const point3d & P, point3d & Q)
{
  assert(points.size() == 8);
  
  point3d Q1, Q2, Q3, Q4, Q5, Q6;
  sVect<point3d> subPoints(4);
  
  subPoints(1) = points(1);
  subPoints(2) = points(4);
  subPoints(3) = points(3);
  subPoints(4) = points(2);
  linearQuad::projectInside(subPoints,P,Q1);
  
  subPoints(1) = points(1);
  subPoints(2) = points(5);
  subPoints(3) = points(8);
  subPoints(4) = points(4);
  linearQuad::projectInside(subPoints,P,Q2);
  
  subPoints(1) = points(1);
  subPoints(2) = points(2);
  subPoints(3) = points(6);
  subPoints(4) = points(5);
  linearQuad::projectInside(subPoints,P,Q3);
  
  subPoints(1) = points(2);
  subPoints(2) = points(3);
  subPoints(3) = points(7);
  subPoints(4) = points(6);
  linearQuad::projectInside(subPoints,P,Q4);
  
  subPoints(1) = points(3);
  subPoints(2) = points(4);
  subPoints(3) = points(8);
  subPoints(4) = points(7);
  linearQuad::projectInside(subPoints,P,Q5);
  
  subPoints(1) = points(5);
  subPoints(2) = points(6);
  subPoints(3) = points(7);
  subPoints(4) = points(8);
  linearQuad::projectInside(subPoints,P,Q6);
  
  Real d1 = point3d::norm2(P-Q1);
  Real d2 = point3d::norm2(P-Q2);
  Real d3 = point3d::norm2(P-Q3);
  Real d4 = point3d::norm2(P-Q4);
  Real d5 = point3d::norm2(P-Q5);
  Real d6 = point3d::norm2(P-Q6);
  Real d  =  std::min(std::min(std::min(d1,d2), std::min(d3,d4)), std::min(d5,d6));
  
  bool flag1 = (d1 == d);
  bool flag2 = (d2 == d);
  bool flag3 = (d3 == d);
  bool flag4 = (d4 == d);
  bool flag5 = (d5 == d);
  bool flag6 = (d6 == d);
  
  UInt num = UInt(flag1) + UInt(flag2) + UInt(flag3) + UInt(flag4) + UInt(flag5) + UInt(flag6);
  assert(num != 0);
  
  if(num != 1)
  {
    if(flag1) { flag2 = false; flag3 = false; flag4 = false; flag5 = false; flag6 = false; }
    if(flag2) { flag3 = false; flag4 = false; flag5 = false; flag6 = false; }
    if(flag3) { flag4 = false; flag5 = false; flag6 = false; }
    if(flag4) { flag5 = false; flag6 = false; }
    if(flag5) { flag6 = false; }
  }
  
  Q = Q1 * Real(flag1) +
      Q2 * Real(flag2) +
      Q3 * Real(flag3) +
      Q4 * Real(flag4) +
      Q5 * Real(flag5) +
      Q6 * Real(flag6);
  
  return(true);
}

point3d
linearHexa::
projection(const point3d & P)
{
  point3d Y;
  
  Y.setX( max(0.0, P.getX()) ); Y.setX( min(1.0, Y.getX()) );
  Y.setY( max(0.0, P.getY()) ); Y.setY( min(1.0, Y.getY()) );
  Y.setZ( max(0.0, P.getZ()) ); Y.setZ( min(1.0, Y.getZ()) );
  
  return(Y);
}

point3d
linearHexa::
getBarycenter()
{
  static point3d Pb = point3d(0.5,0.5,0.5);
  return(Pb);
}

void
linearHexa::
boundingBox(const sVect<point3d> & points, point3d & Pmin, point3d & Pmax)
{
  assert(points.size() == 8);
  
  Pmin = points(1);
  Pmax = points(1);
  
  for(UInt i=1; i <= 8; ++i)
  {
    Pmin.setX( std::min(Pmin.getX(), points(i).getX()) );
    Pmin.setY( std::min(Pmin.getY(), points(i).getY()) );
    Pmin.setZ( std::min(Pmin.getZ(), points(i).getZ()) );
    
    Pmax.setX( std::max(Pmax.getX(), points(i).getX()) );
    Pmax.setY( std::max(Pmax.getY(), points(i).getY()) );
    Pmax.setZ( std::max(Pmax.getZ(), points(i).getZ()) );
  }
}
