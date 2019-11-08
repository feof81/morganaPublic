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

#include "pVect.hpp"
#include "pVectManip.hpp"
#include "pMapItem.h"
#include "pVectComm.hpp"
#include "pVectGlobalManip.hpp"


//! Run with one processor
int main(int argc, char *argv[])
{
  point3d P1(0.0,1.0,0.0);
  point3d P2(1.0,1.0,1.0);
  
  point3d Pmin = std::min(P1,P2);
  point3d Pmax = std::max(P1,P2);
  cout << Pmax << endl;
  
  UInt n = 2;
  dataSegmentationUtility<point3d> segmenter(n,Pmin,Pmax);
  
  point3d P(0.5,1.0,0.0);
  cout << segmenter.index(P2) << endl;
}

/*
1 1 1

2
*/
