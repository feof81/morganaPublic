/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "time.h"
#include <cmath>
#include <iostream>

#include "typesInterface.hpp"
#include "polyCards.h"
#include "polyStatic.hpp"


int main(int argc, char *argv[])
{
  sVect<point3d> points(27);
  
  points(1) = point3d(0.0,0.0,0.0);
  points(2) = point3d(1.0,0.0,0.0);
  points(3) = point3d(1.0,1.0,0.0);
  points(4) = point3d(0.0,1.0,0.0);
  points(5) = point3d(0.5,0.0,0.0);
  points(6) = point3d(1.0,0.5,0.0);
  points(7) = point3d(0.5,1.0,0.0);
  points(8) = point3d(0.0,0.5,0.0);
  points(9) = point3d(0.5,0.5,0.0);
  
  points(10) = point3d(0.0,0.0,0.5);
  points(11) = point3d(1.0,0.0,0.5);
  points(12) = point3d(1.0,1.0,0.5);
  points(13) = point3d(0.0,1.0,0.5);
  points(14) = point3d(0.5,0.0,0.5);
  points(15) = point3d(1.0,0.5,0.5);
  points(16) = point3d(0.5,1.0,0.5);
  points(17) = point3d(0.0,0.5,0.5);
  points(18) = point3d(0.5,0.5,0.5);
  
  points(19) = point3d(0.0,0.0,1.0);
  points(20) = point3d(1.0,0.0,1.0);
  points(21) = point3d(1.0,1.0,1.0);
  points(22) = point3d(0.0,1.0,1.0);
  points(23) = point3d(0.5,0.0,1.0);
  points(24) = point3d(1.0,0.5,1.0);
  points(25) = point3d(0.5,1.0,1.0);
  points(26) = point3d(0.0,0.5,1.0);
  points(27) = point3d(0.5,0.5,1.0);
  
  polyStatic<Q2_3d_BI> polyEval;
  
  for(UInt i=1; i <= 27; ++i)
  {
    cout << polyEval.evaluateStatic(points(i)) << endl;
  }
}
