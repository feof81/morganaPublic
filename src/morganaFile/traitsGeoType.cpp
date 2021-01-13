/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "traitsGeoType.h"


//-------------------------------------------------------------------------------------------------
/*! LinearLine trait */

traitsGeoType<linearLine>::
traitsGeoType()
{ }

string
traitsGeoType<linearLine>::
parariewString()
{
  return("line");
}

void
traitsGeoType<linearLine>::
hdft5String(ofstream & out)
{
  out << "TopologyType=\"Polyline\" NodesPerElement=\"2\"";
}



//-------------------------------------------------------------------------------------------------
/*! LinearTriangle trait */

traitsGeoType<linearTriangle>::
traitsGeoType()
{ }

string
traitsGeoType<linearTriangle>::
parariewString()
{
  return("tri");
}

void
traitsGeoType<linearTriangle>::
hdft5String(ofstream & out)
{
  out << "TopologyType=\"Triangle\"";
}



//-------------------------------------------------------------------------------------------------
/*! LinearQuad trait */

traitsGeoType<linearQuad>::
traitsGeoType()
{ }

string
traitsGeoType<linearQuad>::
parariewString()
{
  return("quad");
}

void
traitsGeoType<linearQuad>::
hdft5String(ofstream & out)
{
  out << "TopologyType=\"Quadrilateral\"";
}



//-------------------------------------------------------------------------------------------------
/*! LinearTetra trait */

traitsGeoType<linearTetra>::
traitsGeoType()
{ }

string
traitsGeoType<linearTetra>::
parariewString()
{
  return("tet");
}

void
traitsGeoType<linearTetra>::
hdft5String(ofstream & out)
{
  out << "TopologyType=\"Tetrahedron\"";
}



//-------------------------------------------------------------------------------------------------
/*! LinearHexa trait */

traitsGeoType<linearHexa>::
traitsGeoType()
{ }

string
traitsGeoType<linearHexa>::
parariewString()
{
  return("hex");
}

void
traitsGeoType<linearHexa>::
hdft5String(ofstream & out)
{
  out << "TopologyType=\"Hexahedron\"";
}
