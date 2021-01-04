/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MORGANAGEOMETRY_HPP
#define MORGANAGEOMETRY_HPP

enum ShapeName         {GEOPOINT = 1, 
LINEARLINE     = 2, QUADRATICLINE = 3,
LINEARTRIANGLE = 4, QUADRATICTRIANGLE = 5, LINEARQUAD = 6,  QUADRATICQUAD = 7,
LINEARTETRA    = 8, QUADRATICTETRA    = 9, LINEARHEXA = 10, QUADRATICHEXA = 11,
LINEARWEDGE    = 12
};

enum ReferenceShapes   {NONE, POINT, LINE, TRIANGLE, QUAD, HEXA, TETRA};

enum ReferenceGeometry {VERTEX = 1, EDGE = 2, FACE = 3, VOLUME = 4};

/*! STDN no standard, STDA standard distribution A (with overlap), STDB standard distribution B (without overlap), STDL standard distribution load (the mesh is mized on the processes),
STDU standard distribution unique (everything on a process)*/
enum MeshStandards     {STDN = 0, STDA = 1, STDB = 2, STDL = 3, STDU = 4};

#endif
