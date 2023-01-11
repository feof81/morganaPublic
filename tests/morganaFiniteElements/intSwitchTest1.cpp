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

#include "feQr3d.hpp"
#include "feSpectralLH3d.hpp"

#include "dofMapStatic3d.hpp"
#include "elCardFeeder3d.hpp"

#include "feDynamicField3d.hpp"
#include "feStaticField3d.hpp"
#include "feStaticField3dGlobalManip.hpp"

#include "opEL3dA.hpp"

#include "intTableEL3d.hpp"
#include "intSwitch.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;

  assert(world.size() == 2);
  
  typedef pMapItemShare             PMAPTYPE;
  typedef feSpectralLH3d<PMAPTYPE>  FIELD_FETYPE;
  typedef Real                      FIELD_DOFTYPE;
  typedef feSpectralLH3d<PMAPTYPE>  TEST_FETYPE;
  typedef Real                      TEST_DOFTYPE;
  typedef feQr3d<0,PMAPTYPE>        COEFF_FETYPE;
  typedef Real                      COEFF_DOFTYPE;
  
  typedef typename FIELD_FETYPE::FECARD FIELD_FECARD;
  typedef typename TEST_FETYPE::FECARD  TEST_FECARD;

  typedef typename FIELD_FETYPE::GEOSHAPE GEOSHAPE;

  
  //Testing
  static const intClass currentIntFlag = intDefault;
  //static const INTclass currentIntFlag = intSTD3d;
  
  static const intClass defaultIntFlag = intTableEL3d<FIELD_FETYPE::feBaseLabel, TEST_FETYPE::feBaseLabel>::intFlag;
  static const intClass finalIntFlag   = intSwitch<currentIntFlag,defaultIntFlag>::intFlag;
  
  cout << finalIntFlag << endl;
}
