/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <cmath>
#include <iostream>

#include "typesInterface.hpp"
#include "simpleFormats.hpp"


int main(int argc, char *argv[])
{
  UInt N = 3000000000;
  time_t start, end;
  
  boost::array<point3d,3> statVector;
  vector<point3d>         dinVector(3);
  
  
  //Static vector
  time(&start);
  
  for(UInt i=1; i <= N; ++i)
  {
    statVector[0] = statVector[0] + statVector[2];
    statVector[1] = statVector[1] + statVector[0];
    statVector[2] = statVector[2] + statVector[1];
  }
  
  time(&end);
  cout << "Static Vector: " << difftime(end,start) << " val " << statVector[2] << endl;
  
  
  //Dynamic vector
  time(&start);
  
  for(UInt i=1; i <= N; ++i)
  {
    dinVector[0] = dinVector[0] + dinVector[2];
    dinVector[1] = dinVector[1] + dinVector[0];
    dinVector[2] = dinVector[2] + dinVector[1];
  }
  
  time(&end);
  cout << "Din Vector: " << difftime(end,start) << " val " << dinVector[2] << endl;
}
