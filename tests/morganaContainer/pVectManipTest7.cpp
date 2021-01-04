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
#include "traitsSegmentationUtility.hpp"


//! Run with one processor - serial
int main(int argc, char *argv[])
{
  typedef pVect<Real,pMapItem> PVECT;
  Teuchos::RCP<PVECT> pVector(new PVECT);    
  pMapItem item;
  sVect<PVECT> targetVects;
  
  item.setPid(0); item.setLid(3); item.setGid(5); pVector->push_back(item,5.0);
  item.setPid(1); item.setLid(2); item.setGid(6); pVector->push_back(item,6.0);
  item.setPid(0); item.setLid(4); item.setGid(4); pVector->push_back(item,4.0);
  item.setPid(1); item.setLid(1); item.setGid(7); pVector->push_back(item,7.0);
  pVector->updateFinder();
  
   
  pVectManip<Real,pMapItem> manipulator(pVector);
  manipulator.segmentationData(targetVects,3, 4.0,7.0);
  
  
  for(UInt i=1; i <= targetVects.size(); ++i)
  {
    cout << "-------------------------------" << endl;
    cout << targetVects(i) << endl << endl;
  }
  
  
}

/*
-------------------------------
1
 map:  pid: 0 lid: 3 gid: 4
 data: 4


-------------------------------
1
 map:  pid: 0 lid: 1 gid: 5
 data: 5


-------------------------------
1
 map:  pid: 1 lid: 2 gid: 6
 data: 6
2
 map:  pid: 1 lid: 4 gid: 7
 data: 7
*/
