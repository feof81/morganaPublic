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

#include "pMapItem.h"
#include "pMapItemShare.h"
#include "pMap.hpp"
#include "pMapGlobalManip.h"
#include "pVectGlobalManip.hpp"

using namespace std;
using namespace boost::mpi;

 
int main(int argc, char *argv[])
{
  //Comm dev---------------------------------------------------------
  boost::mpi::environment  env(argc, argv);
  boost::mpi::communicator oldCommDev;
  
  bool commActive = (oldCommDev.rank() < 2);
  boost::mpi::communicator newCommDev = oldCommDev.split(commActive);
  
  assert(oldCommDev.size() == 4);
  
  //Build map--------------------------------------------------------
  pVect<Real,pMapItem> oldVect;
  
  if(oldCommDev.rank() == 0)
  {
    oldVect.push_back(pMapItem(1,1,oldCommDev.rank()),1);
    oldVect.push_back(pMapItem(2,2,oldCommDev.rank()),2);
  }
  
  if(oldCommDev.rank() == 1)
  {
    oldVect.push_back(pMapItem(1,3,oldCommDev.rank()),3);
    oldVect.push_back(pMapItem(2,4,oldCommDev.rank()),4);
  }
  
  if(oldCommDev.rank() == 2)
  {
    oldVect.push_back(pMapItem(1,5,oldCommDev.rank()),5);
    oldVect.push_back(pMapItem(2,6,oldCommDev.rank()),6);
  }
  
  if(oldCommDev.rank() == 3)
  {
    oldVect.push_back(pMapItem(1,7,oldCommDev.rank()),7);
    oldVect.push_back(pMapItem(2,8,oldCommDev.rank()),8);
  }
  
  //Reduce map-------------------------------------------------------
  if(commActive)
  {    
    pVect<Real,pMapItem> newVect;
    
    pVectGlobalManip<Real,pMapItem> pManip(oldCommDev);
    pManip.reduceCommunicator(commActive,
                              oldCommDev,
                              oldVect,
                              newCommDev,
                              newVect);
    
    for(UInt pid=0; pid < newCommDev.size(); pid++)
    {
      if(pid == newCommDev.rank())
      { cout << newVect << endl; }
      
      newCommDev.barrier();
    }
  }
}
