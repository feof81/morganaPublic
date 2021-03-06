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

using namespace boost;
using namespace std;

int main(int argc, char *argv[])
{
  staticVector<5> vectA, vectB, vectC;
  
  vectA(1) = 1.0;
  vectA(3) = 4.0;
  
  vectB(1) = 1.0;
  vectB(5) = 2.0;
  
  vectC = vectA + vectB;
  cout << vectC << endl;
  
  vectC = vectA - vectB;
  cout << vectC << endl;
  
  vectC = vectA * 2.0;
  cout << vectC << endl;
  
  vectC = vectA / 2.0;
  cout << vectC << endl;
}

