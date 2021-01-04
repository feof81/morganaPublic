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


int main(int argc, char *argv[])
{
  time_t start, end;
  long int top = 1e9;
  
  //point3d----------------------------------------------------------------------------------------
  point3d P, Q;
  
  
  cout << "Expression point3d" << endl; time(&start);
  for(long int i=1; i<=top; ++i)
  {
    P = (P + Q) / 2.0 - (P * 2.0);
  }
  time(&end);
  cout << difftime(end,start) << endl;
  
  
  cout << "Bench point3d" << endl; time(&start);
  for(long int i=1; i<=top; ++i)
  {
    P.X[0] = (P.X[0] + Q.X[0]) / 2.0 - (P.X[0] * 2.0);
    P.X[1] = (P.X[1] + Q.X[1]) / 2.0 - (P.X[1] * 2.0);
    P.X[2] = (P.X[2] + Q.X[2]) / 2.0 - (P.X[2] * 2.0);
  }
  time(&end);
  cout << difftime(end,start) << endl;
}
