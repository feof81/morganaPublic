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

//! Run with one processor - serial
int main(int argc, char *argv[])
{
  typedef pVect<Real,pMapItem> PVECT;
  Teuchos::RCP<PVECT> pVector(new PVECT);    
  pMapItem item;
  
  item.setPid(0); item.setLid(1); item.setGid(4); pVector->push_back(item,1.0,false);
  item.setPid(0); item.setLid(2); item.setGid(5); pVector->push_back(item,2.0,false);
  item.setPid(1); item.setLid(3); item.setGid(6); pVector->push_back(item,3.0,false);
  item.setPid(1); item.setLid(4); item.setGid(7); pVector->push_back(item,4.0,false);
  
  item.setPid(0); item.setLid(5); item.setGid(4); pVector->push_back(item,1.0,false);  
  item.setPid(0); item.setLid(6); item.setGid(5); pVector->push_back(item,2.0,false);
  item.setPid(1); item.setLid(7); item.setGid(6); pVector->push_back(item,3.0,false);
  item.setPid(1); item.setLid(8); item.setGid(7); pVector->push_back(item,4.0,false);
  
  pVectManip<Real,pMapItem> manipulator(pVector);
  manipulator.unionExclusive();
  
  cout << *pVector << endl;
}


/*
1
 map:  pid: 0 lid: 1 gid: 4
 data: 1
2
 map:  pid: 0 lid: 2 gid: 5
 data: 2
3
 map:  pid: 1 lid: 3 gid: 6
 data: 3
4
 map:  pid: 1 lid: 4 gid: 7
 data: 4
*/

