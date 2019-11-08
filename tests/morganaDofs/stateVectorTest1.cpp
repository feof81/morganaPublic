/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "typesInterface.hpp"
#include <cmath>

using namespace std;

int main(int argc, char *argv[])
{
  stateVector A(4), B(4), C(4);
  
  A(1) = 1.0;
  A(2) = 2.0;
  A(3) = 3.0;
  A(4) = 4.0;
  
  B(1) = 1.0;
  B(2) = 2.0;
  B(3) = 1.0;
  B(4) = 2.0;
  
  C = A + B;
  C.cancel();
  
  C += A;
  cout << C;
  
  C = A-B;
  cout << C;
  
  C = A / 2.0;
  C = A * 2.0;
  cout << C;
}
