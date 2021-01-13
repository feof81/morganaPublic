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
  typedef pMapItem                 PMAPTYPE;
  typedef feXFEM_LT3d<1,PMAPTYPE>  FETYPE;
  typedef typename FETYPE::FECARD  FECARD;
  typedef typename FETYPE::ELCARD  ELCARD;
  
  
  //Cards
  FECARD FeCard;
  ELCARD ElCard;
  
  sVect<point3d> nodes(4);
  nodes(1) = point3d(0.0, 1.0, 0.0);
  nodes(2) = point3d(0.0, 0.0, 0.0);
  nodes(3) = point3d(1.0, 1.0, 0.0);
  nodes(4) = point3d(0.0, 1.0, 1.0);
  
  ElCard.setNodes(nodes);
  
  FeCard.setActive(true);
  FeCard.setXtoll(0.1);
  
  FeCard.getPhi(1)  = 1.0;
  FeCard.getPhi(2)  = 0.0;
  FeCard.getPhi(3)  = 0.5;
  FeCard.getPhi(4)  = 1.0;
  
  //Testing the supporter
  sVect<Xvolumes>         volLabels;
  sVect< sVect<point3d> > points;
  
  xfemSupportLT3d<PMAPTYPE> supporter;
  sVect<Real> volumes = supporter.volumes(FeCard,ElCard,volLabels,points);
  
  cout << "Labels" << endl;
  cout << volLabels << endl << endl;
  
  cout << "Points" << endl;
  cout << points << endl;
  
  cout << "volumes: " << endl << volumes << endl;
}
