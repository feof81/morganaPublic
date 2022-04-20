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
  typedef Tpetra::Map<> MAP;
  
  TPETRA_GLOBAL_ORDINAL numRowA, numRowB;
  
  if(world.rank() == 0) { numRowA = 3; numRowB = 2; }
  else                  { numRowA = 2; numRowB = 2; }
  
  std::vector<TPETRA_GLOBAL_ORDINAL> elementsRowA(numRowA);
  std::vector<TPETRA_GLOBAL_ORDINAL> elementsRowB(numRowB);
  
  if(world.rank() == 0)
  {
    elementsRowA[0] = 2;
    elementsRowA[1] = 3;
    elementsRowA[2] = 4;
    
    elementsRowB[0] = 2;
    elementsRowB[1] = 3;
  }
  else
  {
    elementsRowA[0] = 0;
    elementsRowA[1] = 1;
    
    elementsRowB[0] = 0;
    elementsRowB[1] = 1;
  }
  
  TPETRA_GLOBAL_TYPE globSize = Teuchos::OrdinalTraits<TPETRA_GLOBAL_TYPE>::invalid();
  Teuchos::ArrayView<TPETRA_GLOBAL_ORDINAL> teuchosRowA(elementsRowA);
  Teuchos::ArrayView<TPETRA_GLOBAL_ORDINAL> teuchosRowB(elementsRowB);
  
  Teuchos::RCP<MAP> rowMapA(new MAP(globSize,teuchosRowA,0,teuchosComm));
  Teuchos::RCP<MAP> rowMapB(new MAP(globSize,teuchosRowB,0,teuchosComm));

  //Vectors--------------------------------------------------------------------
  typedef Tpetra::MultiVector<> VECTOR;
  
  TPETRA_GLOBAL_ORDINAL I;
  TPETRA_SCALAR       tot;
  
  Teuchos::RCP<VECTOR> vectorA(new VECTOR(rowMapA,1));
  Teuchos::RCP<VECTOR> vectorB(new VECTOR(rowMapB,1));
  
  if(world.rank() == 0)
  {
    I = 2; tot = 3;
    vectorA->sumIntoGlobalValue(I,0,tot);
  
    I = 3; tot = 2;
    vectorA->sumIntoGlobalValue(I,0,tot);
  
    I = 4; tot = 1;
    vectorA->sumIntoGlobalValue(I,0,tot);
    
    I = 2; tot = 6;
    vectorB->sumIntoGlobalValue(I,0,tot);
  
    I = 3; tot = 5;
    vectorB->sumIntoGlobalValue(I,0,tot);
  }
  else
  {
    I = 0; tot = 5;
    vectorA->sumIntoGlobalValue(I,0,tot);
  
    I = 1; tot = 4;
    vectorA->sumIntoGlobalValue(I,0,tot);
    
    I = 0; tot = 3;
    vectorB->sumIntoGlobalValue(I,0,tot);
  
    I = 1; tot = 8;
    vectorB->sumIntoGlobalValue(I,0,tot);
  }
  
  //Assembler
  tpetraBlockVectorAssemble assembler(2);
  assembler.setTpetraComm(teuchosComm);
  
  assembler.setRowMap(1,rowMapA);
  assembler.setRowMap(2,rowMapB);
  assembler.startup();
  
  assembler.setVector(vectorA,1);
  assembler.setVector(vectorB,2);
  
  auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  
  Teuchos::RCP<MAP> mapT = assembler.getRowMap();
  mapT->describe(*out,Teuchos::VERB_EXTREME);
  
  Teuchos::RCP<VECTOR> vectorT = assembler.getVector();
  vectorT->describe(*out,Teuchos::VERB_EXTREME);
}
