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

#include "typesInterface.hpp"
#include "pMapItemShare.h"
#include "pVect.hpp"
#include "pVectManip.hpp"
#include "pVectComm.hpp"
#include "pVectGlobalManip.hpp"

using namespace Teuchos;


//! Run with two processors
int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  RCP<communicator> world(new communicator);
  
  typedef pVect<point3d,pMapItemShare> PVECT;
  PVECT sVector, rVector;
  sVect<bool> isLocal;
  pMapItemShare item;
  
  assert(world->size() == 2);
  
  
  if(world->rank() == 0)
  {
    isLocal.push_back(false);
    item.setPid(0); item.setLid(1); item.setGid(1); sVector.push_back(item,point3d( 1.0, 0.0, 0.0),false);
    
    isLocal.push_back(true);
    item.setPid(0); item.setLid(2); item.setGid(1); sVector.push_back(item,point3d( 0.0, 1.0, 0.0),false);
    
    isLocal.push_back(false);
    item.setPid(0); item.setLid(3); item.setGid(1); sVector.push_back(item,point3d( 0.0, 0.0, 1.0),false);
    
    isLocal.push_back(false);
    item.setPid(0); item.setLid(4); item.setGid(1); sVector.push_back(item,point3d( 2.0, 0.0, 0.0),false);
    
    isLocal.push_back(false);
    item.setPid(0); item.setLid(5); item.setGid(1); sVector.push_back(item,point3d( 0.0, 2.0, 2.0),false);
    
    isLocal.push_back(false);
    item.setPid(0); item.setLid(6); item.setGid(1); sVector.push_back(item,point3d( 3.0, 0.0, 0.0),false);
    
    sVector.updateFinder();
  }
  
  if(world->rank() == 1)
  {
    isLocal.push_back(false);
    item.setPid(1); item.setLid(1); item.setGid(1); sVector.push_back(item,point3d( 0.0, 0.0, 1.0),false);
    
    isLocal.push_back(false);
    item.setPid(1); item.setLid(2); item.setGid(1); sVector.push_back(item,point3d( 2.0, 0.0, 0.0),false);
    
    isLocal.push_back(false);
    item.setPid(1); item.setLid(3); item.setGid(1); sVector.push_back(item,point3d( 0.0, 2.0, 2.0),false);
    
    isLocal.push_back(false);
    item.setPid(1); item.setLid(4); item.setGid(1); sVector.push_back(item,point3d( 3.0, 0.0, 0.0),false);
    
    isLocal.push_back(true);
    item.setPid(1); item.setLid(5); item.setGid(1); sVector.push_back(item,point3d( 0.0, 3.0, 0.0),false);
    
    isLocal.push_back(false);
    item.setPid(1); item.setLid(6); item.setGid(1); sVector.push_back(item,point3d( 0.0, 0.0, 3.0),false);
    
    sVector.updateFinder();
  }
  
  //Global numbering
  pVectGlobalManip<point3d,pMapItemShare> manipulator(world);
  manipulator.buildGlobalNumbering(sVector,isLocal);
  
  
  world->barrier();
  if(world->rank() == 0)
  { cout << sVector << endl; }
  
  world->barrier();
  if(world->rank() == 1)
  { cout << sVector << endl; }
}

