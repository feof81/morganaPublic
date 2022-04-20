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
#include "Teuchos_ScalarTraits.hpp"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_FEVector.h"

#include "typesInterface.hpp"
#include "pMapItem.h"
#include "pMapItemShare.h"
#include "pVect.hpp"
#include "pVectManip.hpp"
#include "pVectComm.hpp"
#include "pVectGlobalManip.hpp"

#include "epetraVector_to_pVect.h"
#include "epetraBlockVectorAssemble.h"
#include "epetraVectorStrip.h"



//! Run with two processors
int main(int argc, char *argv[])
{
  //Startup--------------------------------------------------------------------
  environment env(argc,argv);
  
  communicator world;
  Epetra_MpiComm epetraComm(MPI_COMM_WORLD);
  
  assert(world.size() == 2);
  
  
  //Maps-----------------------------------------------------------------------
  int numRowA, numRowB;
  
  if(world.rank() == 0) { numRowA = 3; numRowB = 2; }
  else                  { numRowA = 2; numRowB = 2; }
  
  int elementsRowA[numRowA];
  int elementsRowB[numRowB];

  
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

  Epetra_Map rowMapA(-1,numRowA,elementsRowA,0,epetraComm);
  Epetra_Map rowMapB(-1,numRowB,elementsRowB,0,epetraComm);
  
  
  //Vectors--------------------------------------------------------------------
  int I;
  double tot;
  
  Epetra_MultiVector vectorA(rowMapA,1);
  Epetra_MultiVector vectorB(rowMapB,1);
  
  if(world.rank() == 0)
  {
    I = 2; tot = 3;
    vectorA.SumIntoGlobalValue(I,0,tot);
  
    I = 3; tot = 2;
    vectorA.SumIntoGlobalValue(I,0,tot);
  
    I = 4; tot = 1;
    vectorA.SumIntoGlobalValue(I,0,tot);
    
    I = 2; tot = 6;
    vectorB.SumIntoGlobalValue(I,0,tot);
  
    I = 3; tot = 5;
    vectorB.SumIntoGlobalValue(I,0,tot);
  }
  else
  {
    I = 0; tot = 5;
    vectorA.SumIntoGlobalValue(I,0,tot);
  
    I = 1; tot = 4;
    vectorA.SumIntoGlobalValue(I,0,tot);
    
    I = 0; tot = 3;
    vectorB.SumIntoGlobalValue(I,0,tot);
  
    I = 1; tot = 8;
    vectorB.SumIntoGlobalValue(I,0,tot);
  }
  
  
  //Assembler
  epetraBlockVectorAssemble assembler(2);
  assembler.setEpetraComm(epetraComm);
  
  assembler.setRowMap(1,rowMapA);
  assembler.setRowMap(2,rowMapB);
  assembler.startup();
  
  assembler.setVector(vectorA,1);
  assembler.setVector(vectorB,2);
  
  typedef epetraBlockVectorAssemble::RCP_VECTOR RCP_VECTOR;
  RCP_VECTOR globVector = assembler.getVector();
  
  if(world.rank() == 0)
  { cout << "The merged vector" << endl; }
  //cout << *globVector << endl;
  
  
  //Stripper
  epetraVectorStrip epetraVectorStrip(2); 
  epetraVectorStrip.setEpetraComm(epetraComm);
  
  epetraVectorStrip.setRowMap(1,rowMapA);
  epetraVectorStrip.setRowMap(2,rowMapB);
  epetraVectorStrip.startup();
  
  cout << "Part of the vector" << endl;
  cout << epetraVectorStrip.getVectorRef(*globVector,1) << endl;
}
