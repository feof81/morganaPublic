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

#include "feStaticEvalIterators.hpp"

int main(int argc, char *argv[])
{
  point3d Y(0.25, 0.25, 0.25);
  
  //Testing RT0 eval
  const UInt numBasisRT0 = 4;
  sVect<point3d> listRT0(numBasisRT0);
  
  feStaticPoint3dEval<FE_RT0LT_3d, numBasisRT0>::eval(listRT0,Y);
  
  cout << "List RT0" << endl;
  cout << listRT0 << endl;
  
  
  //Testing RT0 gradX
  sVect<point3d> listRT0gradX(numBasisRT0);
  
  feStaticPoint3dDx<FE_RT0LT_3d, numBasisRT0>::eval(listRT0gradX,Y);
  
  cout << "List RT0 - Grad X" << endl;
  cout << listRT0gradX << endl;
  
  
  //Testing derivative 
  sVect<point3d> listRT0d(numBasisRT0);
  
  feStaticPoint3dDerivative<FE_RT0LT_3d, numBasisRT0, 0,0,1>::eval(listRT0d,Y);
  
  cout << "List RT0 - Derivative" << endl;
  cout << listRT0d << endl;
}
