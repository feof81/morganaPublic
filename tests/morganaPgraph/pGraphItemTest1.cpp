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

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "pGraphItem.h"


using namespace std;
using namespace boost::mpi;


/*! Run with one processor */
int main(int argc, char *argv[])
{
  pGraphItem P(2), Qa, Qb;
  
  P(1) = 5;
  P(2) = 6;
  
  P.push_back(3);
  cout << "P: " << P.getSorted(1) << " " << P.getSorted(2) << " " << P.getSorted(3) << endl;
  P.push_back(2);
  cout << "P: " << P.getSorted(1) << " " << P.getSorted(2) << " " << P.getSorted(3) << " " << P.getSorted(4) << endl;
  P.push_back(1);
  cout << "P: " << P.getSorted(1) << " " << P.getSorted(2) << " " << P.getSorted(3) << " " << P.getSorted(4) << " " << P.getSorted(5) << endl;
  
  Qa.push_back(3,false);
  Qa.push_back(2,false);
  Qa.push_back(1,false);
  Qa.updateSorting();
  
  Qb.push_back(6,false);
  Qb.push_back(4,false);
  Qb.push_back(2,false);
  Qb.updateSorting();
  
  cout << "Qa: " << Qa.getSorted(1) << " " << Qa.getSorted(2) << " " << Qa.getSorted(3) << endl;
  
  cout << "logic1: " << (Qa < Qb) << endl;
  cout << "logic2: " << (Qb < Qa) << endl;
}

/*
P: 3 5 6
P: 2 3 5 6
P: 1 2 3 5 6
Qa: 1 2 3
logic1: 1
logic2: 0
*/
