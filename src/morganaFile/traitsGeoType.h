/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef TRAITSGEOTYPE_H
#define TRAITSGEOTYPE_H

#include "geoShapes.h"


//-------------------------------------------------------------------------------------------------
/*! Generic traits - empty */
template<typename GEOTYPE>
class traitsGeoType;


//-------------------------------------------------------------------------------------------------
/*! LinearLine trait */
template<> class traitsGeoType<linearLine>
{
  public:
    traitsGeoType();
    static string parariewString();
    void hdft5String(ofstream & out);
};



//-------------------------------------------------------------------------------------------------
/*! LinearTriangle trait */
template<> class traitsGeoType<linearTriangle>
{
  public:
    traitsGeoType();
    static string parariewString();
    void hdft5String(ofstream & out);
};



//-------------------------------------------------------------------------------------------------
/*! LinearQuad trait */
template<> class traitsGeoType<linearQuad>
{
  public:
    traitsGeoType();
    static string parariewString();
    void hdft5String(ofstream & out);
};



//-------------------------------------------------------------------------------------------------
/*! LinearTetra trait */
template<> class traitsGeoType<linearTetra>
{
  public:
    traitsGeoType();
    static string parariewString();
    void hdft5String(ofstream & out);
};



//-------------------------------------------------------------------------------------------------
/*! LinearHexa trait */
template<> class traitsGeoType<linearHexa>
{
  public:
    traitsGeoType();
    static string parariewString();
    void hdft5String(ofstream & out);
};

#endif
