/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "resumeHDF5.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
resumeHDF5::
resumeHDF5()
{
  commDevLoaded = false;
}

resumeHDF5::
resumeHDF5(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

resumeHDF5::
resumeHDF5(communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

void
resumeHDF5::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

void
resumeHDF5::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}


//_________________________________________________________________________________________________
// LOAD AND WRITE FILES TRILINOS
//-------------------------------------------------------------------------------------------------
void
resumeHDF5::
printToFile(const string & s,
            const Epetra_CrsMatrix & A)
{
  //Assert
  assert(commDevLoaded);
  
  //Write
  Epetra_MpiComm epetraComm(*commDev);
  
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Create(s + ".h5");
  hdf5.Write("domainMap",A.DomainMap());
  hdf5.Write("rangeMap",A.RangeMap());
  hdf5.Write("matrix",A);
  hdf5.Close();
}
    
void
resumeHDF5::
loadFromFile(const string & s,
                   Epetra_CrsMatrix & A)
{
  //Assert
  assert(commDevLoaded);
  
  //Alloc
  Epetra_Map * domainMap;
  Epetra_Map * rangeMap;
  Epetra_CrsMatrix * tempA;
  
  //Read
  Epetra_MpiComm epetraComm(*commDev);
  
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Open(s + ".h5");
  hdf5.Read("domainMap",domainMap);
  hdf5.Read("rangeMap",rangeMap);
  hdf5.Read("matrix",*domainMap,*rangeMap,tempA);
  hdf5.Close();
  
  //Trnsfer data
  A = *tempA;
}

