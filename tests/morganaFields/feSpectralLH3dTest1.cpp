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

#include "feSpectralLH3d.hpp"
#include "pMapItem.h"

int main(int argc, char *argv[])
{
  typedef pMapItem                  PMAPTYPE;
  typedef feSpectralLH3d<PMAPTYPE>  FETYPE;
  typedef typename FETYPE::FECARD   FECARD;
  typedef typename FETYPE::ELCARD   ELCARD;
  
  //Element cards
  ELCARD ElCard;
  
  sVect<point3d> nodes(8);
  nodes(1) = point3d(1.0, 1.0, 0.0);
  nodes(2) = point3d(0.0, 1.0, 0.0);
  nodes(3) = point3d(0.0, 0.0, 0.0);
  nodes(4) = point3d(1.0, 0.0, 0.0);
  nodes(5) = point3d(1.0, 1.0, 1.0);
  nodes(6) = point3d(0.0, 1.0, 1.0);
  nodes(7) = point3d(0.0, 0.0, 1.0);
  nodes(8) = point3d(1.0, 0.0, 1.0);
  
  ElCard.setNodes(nodes);
  
  //Finite element cards
  FECARD FeCard;
  FeCard.setIsActive(true);
  FeCard.setR(1,1,1);
  
  //Finite element
  FETYPE finiteElement(FeCard,ElCard);
  UInt numBasis = finiteElement.getNumBasis();
  
  cout << "Dof Cards ------------------" << endl;
  cout << finiteElement.getDofCards() << endl;
  
  
  //Evaluation
  point3d Y(0.0, 0.0, 0.0);
  
  sVect<Real> locEval(numBasis), locGradX(numBasis), locGradY(numBasis), locGradZ(numBasis);
  sVect<Real> globalEval(numBasis), globalGradX(numBasis), globalGradY(numBasis), globalGradZ(numBasis);
  
  finiteElement.localEval(Y, locEval);
  finiteElement.localGrad(Y, locGradX, locGradY, locGradZ);
  finiteElement.globalEval(Y, globalEval);
  finiteElement.globalGrad(Y, globalGradX, globalGradY, globalGradZ);
  
  Real  totLocEval = 0.0,  totLocGradX = 0.0,  totLocGradY = 0.0,  totLocGradZ = 0.0;
  Real totGlobEval = 0.0, totGlobGradX = 0.0, totGlobGradY = 0.0, totGlobGradZ = 0.0;
  
  for(UInt i=1; i <= locEval.size(); ++i)
  {
    totLocEval  += locEval(i);
    totLocGradX += locGradX(i);
    totLocGradY += locGradY(i);
    totLocGradZ += locGradZ(i);
    
    totGlobEval  += globalEval(i);
    totGlobGradX += globalGradX(i);
    totGlobGradY += globalGradY(i);
    totGlobGradZ += globalGradZ(i);
  }
  
  cout << "local  : " << totLocEval << " " <<  totLocGradX << " " << totLocGradY << " " << totLocGradZ << endl;
  cout << "global : " << totGlobEval << " " << totGlobGradX << " " << totGlobGradY << " " << totGlobGradZ << endl;
}

