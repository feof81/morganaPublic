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

#include "tpetraBlockMatrixAssemble.h"

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  //Typedefs
  typedef tpetraBlockMatrixAssemble::RCP_MAP                RCP_MAP;
  typedef tpetraBlockMatrixAssemble::TPETRA_GLOBAL_TYPE     TPETRA_GLOBAL_TYPE;
  typedef tpetraBlockMatrixAssemble::TPETRA_GLOBAL_ORDINAL  TPETRA_GLOBAL_ORDINAL;
  typedef tpetraBlockMatrixAssemble::TPETRA_SCALAR          TPETRA_SCALAR;  
  
  
  //Startup
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  
  RCP<MpiComm<int> > teuchosComm = rcp(new MpiComm<int> (world));
  
  
  //Maps
  TPETRA_GLOBAL_ORDINAL numRowA, numRowB;
  TPETRA_GLOBAL_ORDINAL numColA, numColB;
  
  if(world.rank() == 0) { numRowA = 3; numRowB = 2; }
  else                  { numRowA = 2; numRowB = 2; }
  
  if(world.rank() == 0) { numColA = 2; numColB = 1; }
  else                  { numColA = 1; numColB = 1; }
  
  std::vector<TPETRA_GLOBAL_ORDINAL> elementsRowA(numRowA);
  std::vector<TPETRA_GLOBAL_ORDINAL> elementsRowB(numRowB);
  std::vector<TPETRA_GLOBAL_ORDINAL> elementsColA(numColA);
  std::vector<TPETRA_GLOBAL_ORDINAL> elementsColB(numColB);
  
  if(world.rank() == 0)
  {
    elementsRowA[0] = 2;
    elementsRowA[1] = 3;
    elementsRowA[2] = 4;
    
    elementsRowB[0] = 2;
    elementsRowB[1] = 3;
    
    elementsColA[0] = 0;
    elementsColA[1] = 1;
    
    elementsColB[0] = 0;
  }
  else
  {
    elementsRowA[0] = 0;
    elementsRowA[1] = 1;
    
    elementsRowB[0] = 0;
    elementsRowB[1] = 1;
    
    elementsColA[0] = 2;
    
    elementsColB[0] = 1;
  }
  
  Teuchos::ArrayView<TPETRA_GLOBAL_ORDINAL> teuchosRowA(elementsRowA);
  Teuchos::ArrayView<TPETRA_GLOBAL_ORDINAL> teuchosRowB(elementsRowB);
  Teuchos::ArrayView<TPETRA_GLOBAL_ORDINAL> teuchosColA(elementsColA);
  Teuchos::ArrayView<TPETRA_GLOBAL_ORDINAL> teuchosColB(elementsColB);
  
  TPETRA_GLOBAL_TYPE globSize = Teuchos::OrdinalTraits<TPETRA_GLOBAL_TYPE>::invalid();
  Tpetra::Map<> rowMapA(globSize,teuchosRowA,0,teuchosComm);
  Tpetra::Map<> rowMapB(globSize,teuchosRowB,0,teuchosComm);
  
  Tpetra::Map<> colMapA(globSize,teuchosColA,0,teuchosComm);
  Tpetra::Map<> colMapB(globSize,teuchosColB,0,teuchosComm);
  
  
  //EpetraBlockMatrixAssemble
  tpetraBlockMatrixAssemble assembler(2,2);
  assembler.setTpetraComm(teuchosComm);
  
  assembler.setRowMap(1,rowMapA);
  assembler.setRowMap(2,rowMapB);
  
  assembler.setColMap(1,colMapA);
  assembler.setColMap(2,colMapB);
  
  assembler.startup();
  
  RCP_MAP rowMap = assembler.getRowMap();
  RCP_MAP colMap = assembler.getColMap();
  
  
  auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  
  rowMap->describe(*out,Teuchos::VERB_EXTREME);
  colMap->describe(*out,Teuchos::VERB_EXTREME);
}
