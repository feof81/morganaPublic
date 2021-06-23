/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef GEOTYPETRAITS_HPP
#define GEOTYPETRAITS_HPP

#include "geoShapes.h"


//-------------------------------------------------------------------------------------------------
/*! Generic traits - empty */
template<typename GEOTYPE>
class geotypeTraits;


//-------------------------------------------------------------------------------------------------
/*! LinearLine trait */
template<> class geotypeTraits<linearLine>
{
  public:
    geotypeTraits();
    static string parariewString();
};

geotypeTraits<linearLine>::
geotypeTraits()
{ }

string
geotypeTraits<linearLine>::
parariewString()
{
  return("line");
}


//-------------------------------------------------------------------------------------------------
/*! LinearTriangle trait */
template<> class geotypeTraits<linearTriangle>
{
  public:
    geotypeTraits();
    static string parariewString();
};

geotypeTraits<linearTriangle>::
geotypeTraits()
{ }

string
geotypeTraits<linearTriangle>::
parariewString()
{
  return("tri");
}


//-------------------------------------------------------------------------------------------------
/*! LinearQuad trait */
template<> class geotypeTraits<linearQuad>
{
  public:
    geotypeTraits();
    static string parariewString();
};

geotypeTraits<linearQuad>::
geotypeTraits()
{ }

string
geotypeTraits<linearQuad>::
parariewString()
{
  return("quad");
}


//-------------------------------------------------------------------------------------------------
/*! LinearTetra trait */
template<> class geotypeTraits<linearTetra>
{
  public:
    geotypeTraits();
    static string parariewString();
};

geotypeTraits<linearTetra>::
geotypeTraits()
{ }

string
geotypeTraits<linearTetra>::
parariewString()
{
  return("tet");
}


//-------------------------------------------------------------------------------------------------
/*! LinearTetra trait */
template<> class geotypeTraits<linearHexa>
{
  public:
    geotypeTraits();
    static string parariewString();
};

geotypeTraits<linearHexa>::
geotypeTraits()
{ }

string
geotypeTraits<linearHexa>::
parariewString()
{
  return("hex");
}


#endif
