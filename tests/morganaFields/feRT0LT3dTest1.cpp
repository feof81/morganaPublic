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

#include "feRt0LT3d.hpp"
#include "pMapItem.h"


int main(int argc, char *argv[])
{
  typedef pMapItem                 PMAPTYPE;
  typedef feRt0LT3d<PMAPTYPE>      FETYPE;
  typedef typename FETYPE::FECARD  FECARD;
  typedef typename FETYPE::ELCARD  ELCARD;
  
  
  //Cards
  FECARD FeCard;
  ELCARD ElCard;
  
  sVect<point3d> nodes(4);
  nodes(1) = point3d(0.0, 0.0, 0.0);
  nodes(2) = point3d(1.0, 0.0, 0.0);
  nodes(3) = point3d(0.0, 1.0, 0.0);
  nodes(4) = point3d(0.0, 0.0, 0.5);
  
  sVect<bool> orientation(4);
  orientation(1) = true;
  orientation(2) = true;
  orientation(3) = true;
  orientation(4) = true;
  
  ElCard.setNodes(nodes);
  FeCard.setFacesOrientation(orientation);
  
  
  //Finite element
  FETYPE finiteElement(FeCard,ElCard);  
  
  //Evaluation
  UInt numBasis = FETYPE::numBasis;
  point3d Y(0.5, 0.0, 0.0);
  sVect<point3d> locEval(numBasis), locGradX(numBasis), locGradY(numBasis), locGradZ(numBasis);
  sVect<point3d> globalEval(numBasis), globalGradX(numBasis), globalGradY(numBasis), globalGradZ(numBasis);
  
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
