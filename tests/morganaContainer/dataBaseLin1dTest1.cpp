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

#include "dataBaseLin1d.hpp"

using namespace std;
using namespace boost::mpi;

int main(int argc, char *argv[])
{
  typedef dataBaseLin1d<Real,Real> DATABASE;
  
  sVect<Real> x(5);
  sVect<Real> y(5);
  
  x(1) = 1.0; y(1) = 2.0;
  x(2) = 2.0; y(2) = 4.0;
  x(3) = 3.0; y(3) = 6.0;
  x(4) = 4.0; y(4) = 8.0;
  x(5) = 5.0; y(5) = 10.0;
  
  DATABASE dataBase(x,y); //Add a list of 5 evaluations of the function
  
  cout << "change: " << dataBase.replacePoint(5.0,12.0) << endl; //Replace the value at point x=5.0 : 12.0 replaces 10.0
  cout << "database" << endl;                                    //Prinotout the database
  cout << dataBase << endl;
  
  cout << "min:        " << dataBase.getMinX() << endl;
  cout << "max:        " << dataBase.getMaxX() << endl;
  cout << "isInternal: " << dataBase.isInternal(3.0) << endl;
  cout << "interp:     " << dataBase.interp(4.5) << endl; //The correct value is (8+12) / 2 = 10 since 4.5 is just in the middle of the x-interval [4.0, 5.0]
}
