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

#include "pMapItem.h"
#include "feXFEM_LT3d.hpp"

int main(int argc, char *argv[])
{
  typedef pMapItem                  PMAPTYPE;
  typedef feXFEM_LT3d<2,PMAPTYPE>   FETYPE;
  typedef typename FETYPE::FECARD   FECARD;
  typedef typename FETYPE::ELCARD   ELCARD;
  typedef typename FETYPE::BASETYPE BASETYPE;
  
  //ElCard
  ELCARD ElCard;
  
  sVect<point3d> nodes(4);
  nodes(1) = point3d(0.0, 1.0, 0.0);
  nodes(2) = point3d(0.0, 0.0, 0.0);
  nodes(3) = point3d(1.0, 1.0, 0.0);
  nodes(4) = point3d(0.0, 1.0, 1.0);
  
  
  //FeCard
  FECARD FeCard;
  ElCard.setNodes(nodes);
  
  FeCard.setActive(true);
  FeCard.setXtoll(0.1);
  
  FeCard.getPhi(1)  = 1.0;
  FeCard.getPhi(2)  = 1.0;
  FeCard.getPhi(3)  = 1.0;
  FeCard.getPhi(4)  = 0.0;
  
  FeCard.resizeEnrichedDof(4);
  FeCard.getEnrichedDof(1)  = 1;
  FeCard.getEnrichedDof(2)  = 2;
  FeCard.getEnrichedDof(3)  = 3;
  FeCard.getEnrichedDof(4)  = 4;
  
  //Finite element
  FETYPE finiteElement;
  finiteElement.setCards(FeCard,ElCard);
  
  UInt numBasis = finiteElement.getNumBasis();
  
  cout << "NumBasis : " << numBasis << endl;
  cout << "DofCards" << endl;
  cout << finiteElement.getDofCards() << endl;
  
  sVect<BASETYPE> val(numBasis);
  sVect<BASETYPE> gradX(numBasis), gradY(numBasis), gradZ(numBasis);
  point3d Y(0.0, 0.0, 1.0);
  
  finiteElement.localEval(Y,val);
  finiteElement.localGrad(Y,gradX,gradY,gradZ);
  
  cout << "Local eval  : " << endl << val << endl;
  cout << "Local gradX : " << endl << gradX << endl;
  cout << "Local gradY : " << endl << gradY << endl;
  cout << "Local gradZ : " << endl << gradZ << endl;
  
  finiteElement.globalEval(Y,val);
  finiteElement.globalGrad(Y,gradX,gradY,gradZ);
  
  cout << "Global eval  : " << endl << val << endl;
  cout << "Global gradX : " << endl << gradX << endl;
  cout << "Global gradY : " << endl << gradY << endl;
  cout << "Global gradZ : " << endl << gradZ << endl;
}
