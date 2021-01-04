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


#include "polyPow.hpp"

int main(int argc, char *argv[])
{
  time_t start, end;
  long int top = 10e9;
  Real r = 2.0, x = 0.0;
  
  cout << top << endl;
  
  time(&start);
  for(long int i=1; i<=top; ++i)
  {
    x += polyPow<4>::eval(r);
    //x += polyPosPow<-1>::eval(r);
  }
  time(&end);
  cout << difftime(end,start) << endl;
  cout << x << endl;
  
  x = 0.0;
  time(&start);
  for(long int i=1; i<=top; ++i)
  {
    x += pow(r,4);
  }
  time(&end);
  cout << difftime(end,start) << endl;
  cout << x << endl;
  
  x = 0.0;
  time(&start);
  for(long int i=1; i<=top; ++i)
  {
    x += r * r * r * r;
  }
  time(&end);
  cout << difftime(end,start) << endl;
  cout << x << endl;
}
