/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "komplex.h"
#include "typesInterface.hpp"
#include <cmath>

using namespace std;

int main(int argc, char *argv[])
{
  //Alloc
  komplex A(1.0,1.5), B(0.5,0.75), C;
  
  //Algebraic tests
  C = A * 2.0 + B / 3.0;
  cout << C << endl;
  
  C = A * B;
  cout << C << endl;
  
  C = A / B;
  cout << C << endl;
  
  C = (A - A*2) * (B + B/2);
  cout << C << endl;
  
  //Compare
  cout << (A <  B) << endl;
  cout << (B <  A) << endl;
  cout << (A != B) << endl;
  cout << (A == B) << endl << endl;
  
  C = A;
  
  cout << (A <  C) << endl;
  cout << (C <  A) << endl;
  cout << (A != C) << endl;
  cout << (A == C) << endl << endl;
  
  //Components
  cout << A.getReal() << endl;
  cout << A.getImag() << endl;
  cout << A.norm()    << endl;
  cout << A.phase()   << endl << endl;
  
  cout << komplex::cReal(A) << endl;
  cout << komplex::cImag(A) << endl;
  cout << komplex::norm(A)  << endl;
  cout << komplex::phase(A) << endl << endl;
}
