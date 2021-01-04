/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef SEARCHCARD_H
#define SEARCHCARD_H

#include "typesInterface.hpp"

using namespace std;

class searchCard
{
  public:
    point3d Pmax, Pmin;
    UInt numConnected;
    
  public:
    searchCard();
    void setPmin(const point3d & P);
    void setPmax(const point3d & P);
    void setConnected(const UInt & num);
    const point3d & getPmin() const;
    const point3d & getPmax() const;
    point3d & getPmin();
    point3d & getPmax();
    const UInt & getConnected() const;
    friend ostream & operator<<(ostream & f, const searchCard & C);
};

#endif
