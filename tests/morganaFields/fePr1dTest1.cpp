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

#include "fePr1d.hpp"
#include "pMapItem.h"

int main(int argc, char *argv[])
{
  typedef pMapItem                 PMAPTYPE;
  typedef fePr1d<0,PMAPTYPE>       FETYPE;
  typedef typename FETYPE::FECARD  FECARD;
  typedef typename FETYPE::ELCARD  ELCARD;
  
  
  //Cards
  FECARD FeCard;
  ELCARD ElCard;
  
  sVect<point3d> nodes(2);
  nodes(1) = point3d(0.0, 0.0, 0.0);
  nodes(2) = point3d(1.0, 0.0, 0.0);
  
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
  point3d Y(0.5, 0.0, 0.0);
  sVect<Real> locEval(numBasis), locGradX(numBasis), locGradY(numBasis), locGradZ(numBasis);
  sVect<Real> globalEval(numBasis), globalGradX(numBasis), globalGradY(numBasis), globalGradZ(numBasis);
  
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
