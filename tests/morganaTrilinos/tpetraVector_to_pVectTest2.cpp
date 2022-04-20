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

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map_decl.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"

#include "typesInterface.hpp"
#include "pMapItem.h"
#include "pMapItemShare.h"
#include "pVect.hpp"
#include "pVectManip.hpp"
#include "pVectComm.hpp"
#include "pVectGlobalManip.hpp"

#include "tpetraBlockVectorAssemble.h"
#include "tpetraVector_to_pVect.h"


using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


//! Run with two processors
int main(int argc, char *argv[])
{
  //Typedefs-------------------------------------------------------------------
  typedef tpetraBlockVectorAssemble::RCP_MAP                RCP_MAP;
  typedef tpetraBlockVectorAssemble::TPETRA_GLOBAL_TYPE     TPETRA_GLOBAL_TYPE;
  typedef tpetraBlockVectorAssemble::TPETRA_GLOBAL_ORDINAL  TPETRA_GLOBAL_ORDINAL;
  typedef tpetraBlockVectorAssemble::TPETRA_SCALAR          TPETRA_SCALAR;  
    
  //Startup--------------------------------------------------------------------
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  
  RCP<const Comm<int> > teuchosComm = rcp(new MpiComm<int> (world));
  
  //Maps-----------------------------------------------------------------------
  typedef Tpetra::Map<>         TPETRA_MAP;
  typedef Tpetra::MultiVector<> TPETRA_VECTOR;
  typedef pMapItemShare MAP;
  typedef pMap<MAP>    PMAP;
  
  PMAP pMapRef;
  
  //Epetra map
  TPETRA_GLOBAL_ORDINAL numMyElements;
  
  if(world.rank() == 0) { numMyElements = 3; }
  else                  { numMyElements = 2; }
  
  TPETRA_GLOBAL_ORDINAL myGlobalElements[numMyElements];
  
  if(world.rank() == 0)
  {
    myGlobalElements[0] = 2;
    myGlobalElements[1] = 3;
    myGlobalElements[2] = 4;
    
    pMapRef.resize(3);
    pMapRef(1) = MAP(pMapItemShare(1,3,0,false,true));
    pMapRef(2) = MAP(pMapItemShare(2,4,0,false,true));
    pMapRef(3) = MAP(pMapItemShare(3,5,0,false,true));
  }
  else
  {
    myGlobalElements[0] = 0;
    myGlobalElements[1] = 1;
    
    pMapRef.resize(2);
    pMapRef(1) = MAP(pMapItemShare(1,1,1,false,true));
    pMapRef(2) = MAP(pMapItemShare(2,2,1,false,true));
  }

  Teuchos::RCP<TPETRA_MAP> tpetraMap(new TPETRA_MAP(-1,myGlobalElements,numMyElements,0,teuchosComm));
  
  
  //Epetra vector
  Teuchos::RCP<TPETRA_VECTOR> tpetraVector(new TPETRA_VECTOR(tpetraMap,1));
  
  TPETRA_GLOBAL_ORDINAL I;
  TPETRA_SCALAR       tot;
  
  if(world.rank() == 0)
  {
    I = 2; tot = 2;  tpetraVector->sumIntoGlobalValue(I,0,tot);
    I = 3; tot = 3;  tpetraVector->sumIntoGlobalValue(I,0,tot);
    I = 4; tot = 4;  tpetraVector->sumIntoGlobalValue(I,0,tot);
  }
  else
  {
    I = 0; tot = 0;  tpetraVector->sumIntoGlobalValue(I,0,tot);
    I = 1; tot = 1;  tpetraVector->sumIntoGlobalValue(I,0,tot);
  }
  
  //Estract the map 
  tpetraVector_to_pVect<MAP> converter(world);
  pVect<Real,MAP> mvect = converter.convert(*tpetraVector,pMapRef);
  
  if(world.rank() == 0)
  { cout << mvect << endl; }
}
