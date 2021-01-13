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
  
  //Testing P0
  const UInt numBasisP0 = 1;
  sVect<Real> listP0(numBasisP0);
  
  feStaticRealEval<FE_P0_3d, numBasisP0>::eval(listP0,Y);
  
  cout << "List P0" << endl;
  cout << listP0 << endl;
  
  
  //Testing P1
  const UInt numBasisP1 = 4;
  sVect<Real> listP1(numBasisP1);
  
  feStaticRealEval<FE_P1_3d, numBasisP1>::eval(listP1,Y);
  
  cout << "List P1" << endl;
  cout << listP1 << endl;
  
  
  //Testing Dx P1
  sVect<Real> listDxP1(numBasisP1);
  
  feStaticRealDx<FE_P1_3d, numBasisP1>::eval(listDxP1,Y);
  
  cout << "List Dx P1" << endl;
  cout << listDxP1 << endl;
  
  
  //Testing generic derivative
  sVect<Real> listDDP1(numBasisP1);
  
  feStaticRealDerivative<FE_P1_3d, numBasisP1, 1,0,0>::eval(listDDP1,Y);
  
  cout << "List DD P1" << endl;
  cout << listDDP1 << endl;
}
