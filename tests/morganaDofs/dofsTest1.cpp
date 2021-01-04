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

#include "typesInterface.hpp"
#include "traitsSegmentationUtility.hpp"


int main(int argc, char *argv[])
{
  //Allocations
  point2d P2A(1.0,2.0), P2B(3.0,4.0), P2;
  point3d P3A(1.0,2.0,3.0), P3B(3.0,2.0,1.0), P3;
  
  tensor2d T2A, T2B, T2;
  tensor3d T3A, T3B, T3;
  
  T2A.setIJ(1,1,1.0); T2A.setIJ(1,2,2.0);
  T2A.setIJ(2,1,3.0); T2A.setIJ(2,2,4.0);
  
  T2B.setIJ(1,1,4.0); T2B.setIJ(1,2,3.0);
  T2B.setIJ(2,1,2.0); T2B.setIJ(2,2,1.0);
  
  T3A.setIJ(1,1,1.0); T3A.setIJ(1,2,2.0); T3A.setIJ(1,3,3.0);
  T3A.setIJ(2,1,4.0); T3A.setIJ(2,2,5.0); T3A.setIJ(2,3,6.0);
  T3A.setIJ(3,1,7.0); T3A.setIJ(3,2,8.0); T3A.setIJ(3,3,9.0);
  
  T3B.setIJ(1,1,9.0); T3B.setIJ(1,2,8.0); T3B.setIJ(1,3,7.0);
  T3B.setIJ(2,1,6.0); T3B.setIJ(2,2,5.0); T3B.setIJ(2,3,4.0);
  T3B.setIJ(3,1,3.0); T3B.setIJ(3,2,2.0); T3B.setIJ(3,3,1.0);
  
  stateVector VA(2), VB(2), V(2), Vt(2);
  
  VA(1) = 1.0; VA(2) = 2.0;
  VB(1) = 2.0; VB(2) = 1.0;
  
  stateMatrix MA(2,2), MB(2,2), M(2,2), Mt(2,2);
  
  MA(1,1) = 1.0; MA(1,2) = 2.0;
  MA(2,1) = 3.0; MA(2,2) = 4.0;
  
  staticVector<2> SA, SB, S, St;
  SA(1) = 1.0; SA(2) = 2.0;
  SB(1) = 2.0; SB(2) = 1.0;
  
  
  //Testing
  Real a = 1.0/3.0;
  
  P2 = linearCombination(P2A,P2B,a);
  P3 = linearCombination(P3A,P3B,a);
  T2 = linearCombination(T2A,T2B,a);
  T3 = linearCombination(T3A,T3B,a);
  
  V = linearCombination(VA,VB,a);
  Vt = VA * (1.0-a) + VB * a;
  
  M = linearCombination(MA,MB,a);
  Mt = MA * (1.0-a) + MB * a;
  
  S = linearCombination(SA,SB,a);
  St = SA * (1.0-a) + SB * a;
  
  cout << "point2d: " << !(P2 != (P2A * (1.0-a) + P2B * a)) << endl;
  cout << "point3d: " << !(P3 != (P3A * (1.0-a) + P3B * a)) << endl;
  cout << "tensor2d: " << !(T2 != (T2A * (1.0-a) + T2B * a)) << endl;
  cout << "tensor3d: " << !(T3 != (T3A * (1.0-a) + T3B * a)) << endl;
  
  cout << "stateVector: " << endl;
  cout << V(1) << " " << Vt(1) << endl;
  cout << V(2) << " " << Vt(2) << endl;
  
  cout << "stateMatrix: " << endl;
  cout << M(1,1) << " " << Mt(1,1) << endl;
  cout << M(1,2) << " " << Mt(1,2) << endl;
  cout << M(2,1) << " " << Mt(2,1) << endl;
  cout << M(2,2) << " " << Mt(2,2) << endl;
  
  cout << "staticVector: " << endl;
  cout << S(1) << " " << St(1) << endl;
  cout << S(2) << " " << St(2) << endl;
}
