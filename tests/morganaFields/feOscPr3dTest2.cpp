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

#include "feOscPr3d.hpp"
#include "pMapItem.h"

int main(int argc, char *argv[])
{
  //Typedefs
  typedef pMapItem                 PMAPTYPE;
  typedef feOscPr3d<1,1,PMAPTYPE>  FETYPE;
  typedef typename FETYPE::FECARD  FECARD;
  typedef typename FETYPE::ELCARD  ELCARD;
  
  //FeCard
  const Real pi = 3.1415926535897;
  
  point3d K1( 4* pi,   0.0,   0.0);
  point3d K2(-4* pi,   0.0,   0.0);
  point3d K3(   0.0, 4* pi,   0.0);
  point3d K4(   0.0,-4* pi,   0.0);
  point3d K5(   0.0,   0.0, 4* pi);
  point3d K6(   0.0,   0.0,-4* pi);
  
  point3d Y1(0.0, 0.0, 0.0);
  point3d Y2(1.0, 0.0, 0.0);
  point3d Y3(0.0, 1.0, 0.0);
  point3d Y4(0.0, 0.0, 1.0);
  
  FECARD FeCard(24);
  FeCard.set(Y1,K1,1);
  FeCard.set(Y1,K2,2);
  FeCard.set(Y1,K3,3);
  FeCard.set(Y1,K4,4);
  FeCard.set(Y1,K5,5);
  FeCard.set(Y1,K6,6);
  
  FeCard.set(Y2,K1,7);
  FeCard.set(Y2,K2,8);
  FeCard.set(Y2,K3,9);
  FeCard.set(Y2,K4,10);
  FeCard.set(Y2,K5,11);
  FeCard.set(Y2,K6,12);
  
  FeCard.set(Y3,K1,13);
  FeCard.set(Y3,K2,14);
  FeCard.set(Y3,K3,15);
  FeCard.set(Y3,K4,16);
  FeCard.set(Y3,K5,17);
  FeCard.set(Y3,K6,18);
  
  FeCard.set(Y4,K1,19);
  FeCard.set(Y4,K2,20);
  FeCard.set(Y4,K3,21);
  FeCard.set(Y4,K4,22);
  FeCard.set(Y4,K5,23);
  FeCard.set(Y4,K6,24);
  
  //Cards
  ELCARD ElCard;
  
  sVect<point3d> nodes(4);
  nodes(1) = point3d(0.0, 0.0, 0.0);
  nodes(2) = point3d(1.0, 0.0, 0.0);
  nodes(3) = point3d(0.0, 1.0, 0.0);
  nodes(4) = point3d(0.0, 0.0, 1.0);
  
  ElCard.setNodes(nodes);
  
  //Finite element
  FETYPE finiteElement(FeCard,ElCard);
  
  cout << "Dof Cards ------------------" << endl;
  for(UInt i=1; i<= FETYPE::numBasis; ++i)
  {
    cout << i << endl;
    cout << finiteElement.getDofCard(i) << endl;
  }
  
  //Evaluation
  UInt numBasis = FETYPE::numBasis;
  point3d Y(0.0, 0.0, 0.0);
  sVect<komplex> locEval(numBasis), locGradX(numBasis), locGradY(numBasis), locGradZ(numBasis);
  sVect<komplex> globalEval(numBasis), globalGradX(numBasis), globalGradY(numBasis), globalGradZ(numBasis);
  
  finiteElement.localEval(Y, locEval);
  finiteElement.localGrad(Y, locGradX, locGradY, locGradZ);
  finiteElement.globalEval(Y, globalEval);
  finiteElement.globalGrad(Y, globalGradX, globalGradY, globalGradZ);
  
  
  cout << "LocEval" << endl;
  cout << locEval << endl;
  
  cout << "LocGradX" << endl;
  cout << locGradX << endl;
  
  cout << "LocGradY" << endl;
  cout << locGradY << endl;
  
  cout << "LocGradZ" << endl;
  cout << locGradZ << endl << endl;
  
  
  cout << "GlobalEval" << endl;
  cout << globalEval << endl;
  
  cout << "GlobalGradX" << endl;
  cout << globalGradX << endl;
  
  cout << "GlobalGradY" << endl;
  cout << globalGradY << endl;
  
  cout << "GlobalGradZ" << endl;
  cout << globalGradZ << endl;
}
