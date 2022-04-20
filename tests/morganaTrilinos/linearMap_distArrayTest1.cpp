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

#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "EpetraExt_DistArray.h"
#include "EpetraExt_HDF5.h"
#include "Epetra_Import.h"
#include "linearMap_distArray.hpp"


int main(int argc, char *argv[])
{ 
  MPI_Init(&argc,&argv);
  Epetra_MpiComm epetraComm(MPI_COMM_WORLD);
  int numMyElements;
  int MyGlobalElements[4];
  
  if(epetraComm.MyPID() == 0)
  {
     numMyElements = 4;
     MyGlobalElements[0] = 0;
     MyGlobalElements[1] = 1;
     MyGlobalElements[2] = 2;
     MyGlobalElements[3] = 3;
  }
  
  if(epetraComm.MyPID() == 1)
  {
     numMyElements = 4;
     MyGlobalElements[0] = 4;
     MyGlobalElements[1] = 5;
     MyGlobalElements[2] = 6;
     MyGlobalElements[3] = 7;
  }
  
  if(epetraComm.MyPID() == 2)
  {
     numMyElements = 4;
     MyGlobalElements[0] = 8;
     MyGlobalElements[1] = 9;
     MyGlobalElements[2] = 10;
     MyGlobalElements[3] = 11;
  }
  
  if(epetraComm.MyPID() == 3)
  {
     numMyElements = 2;
     MyGlobalElements[0] = 12;
     MyGlobalElements[1] = 13;
  }
  
  std::string GroupName = "nodes";
  
  Epetra_Map mappa(-1,numMyElements,MyGlobalElements,0,epetraComm);
  EpetraExt::DistArray<double> nodes_epetra(mappa,3);
  
  EpetraExt::DistArray<double> nodesLinear = linearMapDistArray(nodes_epetra,epetraComm);
  
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Create("err.h5");
  hdf5.Write("alias",nodesLinear);
  hdf5.Close();
  
  
  MPI_Finalize();
}
