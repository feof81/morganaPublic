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

#include "geoShapes.h"
#include "meshInit3d.hpp"

#include "traitsGeometry.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;

int main(int argc, char *argv[])
{
  typedef linearTetra            GEOSHAPE3D;
  typedef geoElement<GEOSHAPE3D> GEOELEMENT3D;
  typedef pMapItemShare          ELMAP;
  typedef pMapItemShare          NODEMAP;
  
  
  pVect<point3d,NODEMAP> nodes3d;
    
  nodes3d.reserve(4);
  nodes3d.push_back(point3d(0.0, 0.0, 0.0),pMapItemShare(1,1,true,false));
  nodes3d.push_back(point3d(1.0, 0.0, 0.0),pMapItemShare(2,2,true,false));
  nodes3d.push_back(point3d(1.0, 1.0, 0.0),pMapItemShare(3,3,true,false));
  nodes3d.push_back(point3d(0.5, 0.5, 1.0),pMapItemShare(4,4,true,false));
  nodes3d.updateFinder();
  
  //Elements 3d
  GEOELEMENT3D tet(true);
  pGraph<GEOELEMENT3D,ELMAP,NODEMAP> elList3d;
   
  elList3d.reserve(1);
  tet.setGeoId(1); tet(1) = 1; tet(2) = 2; tet(3) = 3; tet(4) = 4;
  elList3d.push_back(tet, pMapItemShare(1,1,false,false));
  
  elList3d.setColMap(nodes3d.getMapRef());
  
  
  auxiliaryGraphCheck3d<GEOSHAPE3D,ELMAP,NODEMAP> checker;
  cout << "Check: " << checker.check(elList3d) << endl;
}