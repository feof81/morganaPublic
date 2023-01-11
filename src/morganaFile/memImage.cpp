/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "memImage.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS AND FUNCTIONS
//-------------------------------------------------------------------------------------------------
memImage::
memImage()
{
  commDevLoaded = false;
}

memImage::
memImage(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

memImage::
memImage(communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

void
memImage::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

void
memImage::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}


//_________________________________________________________________________________________________
// PRINT TO FILE - DOFS
//-------------------------------------------------------------------------------------------------
void
memImage::
printMem(fstream & f, const bool & A) const
{
  f.write((char*)&A, sizeof(bool));
}

void
memImage::
printMem(fstream & f, const Int & A) const
{
  f.write((char*)&A, sizeof(Int));
}

void
memImage::
printMem(fstream & f, const UInt & A) const
{
  f.write((char*)&A, sizeof(UInt));
}

void
memImage::
printMem(fstream & f, const Real & A) const
{
  f.write((char*)&A, sizeof(Real));
}

void
memImage::
printMem(fstream & f, const komplex & A) const
{
  f.write((char*)&A.getReal(), sizeof(Real));
  f.write((char*)&A.getImag(), sizeof(Real));
}

void
memImage::
printMem(fstream & f, const point2d & A) const
{
  f.write((char*)&A.getX(), sizeof(Real));
  f.write((char*)&A.getY(), sizeof(Real));
  
  UInt id = A.getId();
  f.write((char*)&id, sizeof(UInt));
}

void
memImage::
printMem(fstream & f, const point3d & A) const
{
  f.write((char*)&A.getX(), sizeof(Real));
  f.write((char*)&A.getY(), sizeof(Real));
  f.write((char*)&A.getZ(), sizeof(Real));
  
  UInt id = A.getId();
  f.write((char*)&id, sizeof(UInt));
}

void
memImage::
printMem(fstream & f, const tensor2d & A) const
{
  f.write((char*)&A.getIJ(1,1), sizeof(Real));
  f.write((char*)&A.getIJ(1,2), sizeof(Real));
  f.write((char*)&A.getIJ(2,1), sizeof(Real));
  f.write((char*)&A.getIJ(2,2), sizeof(Real));
  
  UInt id = A.getId();
  f.write((char*)&id, sizeof(UInt));
}

void
memImage::
printMem(fstream & f, const tensor3d & A) const
{
  f.write((char*)&A.getIJ(1,1), sizeof(Real));
  f.write((char*)&A.getIJ(1,2), sizeof(Real));
  f.write((char*)&A.getIJ(1,3), sizeof(Real));
  
  f.write((char*)&A.getIJ(2,1), sizeof(Real));
  f.write((char*)&A.getIJ(2,2), sizeof(Real));
  f.write((char*)&A.getIJ(2,3), sizeof(Real));
  
  f.write((char*)&A.getIJ(3,1), sizeof(Real));
  f.write((char*)&A.getIJ(3,2), sizeof(Real));
  f.write((char*)&A.getIJ(3,3), sizeof(Real));
  
  UInt id = A.getId();
  f.write((char*)&id, sizeof(UInt));
}

void
memImage::
printMem(fstream & f, const stateVector & A) const
{
  UInt length = A.Length();
  f.write((char*)&length, sizeof(UInt));
  
  for(UInt i=1; i <= length; ++i)
  { f.write((char*)&A(i), sizeof(Real)); }
}

void
memImage::
printMem(fstream & f, const stateMatrix & A) const
{
  UInt numRows = A.RowDim();
  UInt numCols = A.ColDim();
  
  f.write((char*)&numRows, sizeof(UInt));
  f.write((char*)&numCols, sizeof(UInt));
  
  for(UInt i=1; i <= numRows; ++i)
  {
    for(UInt j=1; j <= numCols; ++j)
    { f.write((char*)&A(i,j), sizeof(Real)); }
  }
}


//_________________________________________________________________________________________________
// LOAD FROM FILE - DOFS
//-------------------------------------------------------------------------------------------------
void
memImage::
loadMem(fstream & f, bool & A) const
{
  f.read((char*)&A, sizeof(bool));
}

void
memImage::
loadMem(fstream & f, Int & A) const
{
  f.read((char*)&A, sizeof(Int));
}

void
memImage::
loadMem(fstream & f, UInt & A) const
{
  f.read((char*)&A, sizeof(UInt));
}

void
memImage::
loadMem(fstream & f, Real & A) const
{
  f.read((char*)&A, sizeof(Real));
}

void
memImage::
loadMem(fstream & f, komplex & A) const
{
  Real compReal, compImag;
  
  f.read((char*)&compReal, sizeof(Real));
  f.read((char*)&compImag, sizeof(Real));
  
  A.setReal(compReal);
  A.setImag(compImag);
}

void
memImage::
loadMem(fstream & f, point2d & A) const
{
  UInt id;
  Real X, Y;
  
  f.read((char*)&X,  sizeof(Real));
  f.read((char*)&Y,  sizeof(Real));
  f.read((char*)&id, sizeof(UInt));
  
  A.setX(X);
  A.setY(Y);
  A.setId(id);
}

void
memImage::
loadMem(fstream & f, point3d & A) const
{
  UInt id;
  Real X, Y, Z;
  
  f.read((char*)&X,  sizeof(Real));
  f.read((char*)&Y,  sizeof(Real));
  f.read((char*)&Z,  sizeof(Real));
  f.read((char*)&id, sizeof(UInt));
  
  A.setX(X);
  A.setY(Y);
  A.setZ(Z);
  A.setId(id);
}

void
memImage::
loadMem(fstream & f, tensor2d & A) const
{
  UInt id;
  Real T11, T12, T21, T22;
  
  f.read((char*)&T11, sizeof(Real));
  f.read((char*)&T12, sizeof(Real));
  f.read((char*)&T21, sizeof(Real));
  f.read((char*)&T22, sizeof(Real));
  f.read((char*)&id,  sizeof(UInt));
  
  A.setIJ(1,1,T11);
  A.setIJ(1,2,T12);
  A.setIJ(2,1,T21);
  A.setIJ(2,2,T22);
  A.setId(id);
}

void
memImage::
loadMem(fstream & f, tensor3d & A) const
{
  UInt id;
  Real T11, T12, T13, T21, T22, T23, T31, T32, T33;
  
  f.read((char*)&T11, sizeof(Real));
  f.read((char*)&T12, sizeof(Real));
  f.read((char*)&T13, sizeof(Real));
  
  f.read((char*)&T21, sizeof(Real));
  f.read((char*)&T22, sizeof(Real));
  f.read((char*)&T23, sizeof(Real));
  
  f.read((char*)&T31, sizeof(Real));
  f.read((char*)&T32, sizeof(Real));
  f.read((char*)&T33, sizeof(Real));
  
  A.setId(id);
  
  A.setIJ(1,1,T11);
  A.setIJ(1,2,T12);
  A.setIJ(1,3,T13);
  
  A.setIJ(2,1,T21);
  A.setIJ(2,2,T22);
  A.setIJ(2,3,T23);
  
  A.setIJ(3,1,T31);
  A.setIJ(3,2,T32);
  A.setIJ(3,3,T33);
}

void
memImage::
loadMem(fstream & f, stateVector & A) const
{
  UInt length;
  f.read((char*)&length, sizeof(UInt));
  
  stateVector V(length);
  for(UInt i=1; i <= length; ++i)
  { f.read((char*)&V(i), sizeof(Real)); }
  
  A = V;
}

void
memImage::
loadMem(fstream & f, stateMatrix & A) const
{
  UInt numRows;
  f.read((char*)&numRows, sizeof(UInt));
  
  UInt numCols;
  f.read((char*)&numCols, sizeof(UInt));
  
  stateMatrix V(numRows,numCols);
  for(UInt i=1; i <= numRows; ++i)
  {
    for(UInt j=1; j <= numCols; ++j)
    { f.read((char*)&V(i,j), sizeof(Real)); }
  }
  
  A = V;
}


//_________________________________________________________________________________________________
// PRINT TO FILE - CONTAINER
//-------------------------------------------------------------------------------------------------
void
memImage::
printMem(fstream & f, const pMapItem & A) const
{
  f.write((char*)&A.getLid(),    sizeof(UInt));
  f.write((char*)&A.getGid(),    sizeof(UInt));
  f.write((char*)&A.getPid(),    sizeof(UInt));
  f.write((char*)&A.getBufLid(), sizeof(UInt));
}

void
memImage::
printMem(fstream & f, const pMapItemShare & A) const
{
  f.write((char*)&A.getLid(),    sizeof(UInt));
  f.write((char*)&A.getGid(),    sizeof(UInt));
  f.write((char*)&A.getPid(),    sizeof(UInt));
  f.write((char*)&A.getBufLid(), sizeof(UInt));
  
  f.write((char*)&A.getShared(), sizeof(bool));
  f.write((char*)&A.getOwned(),  sizeof(bool));
}

void
memImage::
printMem(fstream & f, const pMapItemSendRecv & A) const
{
  f.write((char*)&A.getLid(),    sizeof(UInt));
  f.write((char*)&A.getGid(),    sizeof(UInt));
  f.write((char*)&A.getPid(),    sizeof(UInt));
  f.write((char*)&A.getBufLid(), sizeof(UInt));
  
  f.write((char*)&A.getSid(), sizeof(UInt));
  f.write((char*)&A.getRid(), sizeof(UInt));
}

void
memImage::
printMem(fstream & f, const sOctTreeItem & A) const
{
  f.write((char*)&A.getFather(), sizeof(UInt));
  f.write((char*)&A.getIsLeaf(), sizeof(bool));
  f.write((char*)&A.getId(),     sizeof(UInt));
  
  for(UInt i=1; i <= 8; ++i)
  { f.write((char*)&A.getSons(i), sizeof(UInt)); }
}


//_________________________________________________________________________________________________
// LOAD FROM FILE - CONTAINER
//-------------------------------------------------------------------------------------------------
void
memImage::
loadMem(fstream & f, pMapItem & A) const
{
  UInt lid;    f.read((char*)&lid,    sizeof(UInt));
  UInt gid;    f.read((char*)&gid,    sizeof(UInt));
  UInt pid;    f.read((char*)&pid,    sizeof(UInt));
  UInt bufLid; f.read((char*)&bufLid, sizeof(UInt));
  
  A.setLid(lid);
  A.setGid(gid);
  A.setPid(pid);
  A.setBufLid(bufLid);
}

void
memImage::
loadMem(fstream & f, pMapItemShare & A) const
{
  UInt lid;    f.read((char*)&lid,    sizeof(UInt));
  UInt gid;    f.read((char*)&gid,    sizeof(UInt));
  UInt pid;    f.read((char*)&pid,    sizeof(UInt));
  UInt bufLid; f.read((char*)&bufLid, sizeof(UInt));
  bool shared; f.read((char*)&shared, sizeof(bool));
  bool owned;  f.read((char*)&owned,  sizeof(bool));
  
  A.setLid(lid);
  A.setGid(gid);
  A.setPid(pid);
  A.setBufLid(bufLid);
  A.setShared(shared);
  A.setOwned(owned);
}

void
memImage::
loadMem(fstream & f, pMapItemSendRecv & A) const
{
  UInt lid;    f.read((char*)&lid,    sizeof(UInt));
  UInt gid;    f.read((char*)&gid,    sizeof(UInt));
  UInt pid;    f.read((char*)&pid,    sizeof(UInt));
  UInt bufLid; f.read((char*)&bufLid, sizeof(UInt));
  UInt sid;    f.read((char*)&sid,    sizeof(UInt));
  UInt rid;    f.read((char*)&rid,    sizeof(UInt));
  
  A.setLid(lid);
  A.setGid(gid);
  A.setPid(pid);
  A.setBufLid(bufLid);
  A.setSid(sid);
  A.setRid(rid);
}

void
memImage::
loadMem(fstream & f, sOctTreeItem & A) const
{
  f.read((char*)&A.father, sizeof(UInt));
  f.read((char*)&A.isLeaf, sizeof(bool));
  f.read((char*)&A.id,     sizeof(UInt));
  
  for(UInt i=1; i <= 8; ++i)
  { f.read((char*)&A.sons[i-1], sizeof(UInt)); }
}


//_________________________________________________________________________________________________
// PRINT TO FILE - PGRAPH
//-------------------------------------------------------------------------------------------------
void
memImage::
printMem(fstream & f, const pGraphItem & A) const
{
  sVect<UInt> connected    = A.connected;
  sVect<UInt> orderedNodes = A.orderedNodes;
  
  printMem(f, connected);
  printMem(f, orderedNodes);
  printMem(f, A.nodesOrdered);
}

void
memImage::
printMem(fstream & f, const pGraphItemOriented & A) const
{
  sVect<UInt> connected    = A.connected;
  sVect<UInt> orderedNodes = A.orderedNodes;
  sVect<bool> orientation  = A.orientation;
  
  printMem(f, connected);
  printMem(f, orderedNodes);
  printMem(f, orientation);
  printMem(f, A.nodesOrdered);
}

void
memImage::
printMem(fstream & f, const pGraphItemSubLoc & A) const
{
  sVect<UInt> connected    = A.connected;
  sVect<UInt> orderedNodes = A.orderedNodes;
  sVect<UInt> subLocId     = A.subLocId;
  
  printMem(f, connected);
  printMem(f, orderedNodes);
  printMem(f, subLocId);
  printMem(f, A.nodesOrdered);
}


//_________________________________________________________________________________________________
// LOAD FROM FILE - PGRAPH
//-------------------------------------------------------------------------------------------------
void
memImage::
loadMem(fstream & f, pGraphItem & A) const
{
  sVect<UInt> connected;
  sVect<UInt> orderedNodes;
  bool        nodesOrdered;
  
  loadMem(f, connected);
  loadMem(f, orderedNodes);
  loadMem(f, nodesOrdered);
  
  A.connected    = connected;
  A.orderedNodes = orderedNodes;
  A.nodesOrdered = nodesOrdered;
}

void
memImage::
loadMem(fstream & f, pGraphItemOriented & A) const
{
  sVect<UInt> connected;
  sVect<UInt> orderedNodes;
  sVect<bool> orientation;
  bool        nodesOrdered;
  
  loadMem(f, connected);
  loadMem(f, orderedNodes);
  loadMem(f, orientation);
  loadMem(f, nodesOrdered);
  
  A.connected    = connected;
  A.orderedNodes = orderedNodes;
  A.orientation  = orientation;
  A.nodesOrdered = nodesOrdered;
}

void
memImage::
loadMem(fstream & f, pGraphItemSubLoc & A) const
{
  sVect<UInt> connected;
  sVect<UInt> orderedNodes;
  sVect<UInt> subLocId;
  bool        nodesOrdered;
  
  loadMem(f, connected);
  loadMem(f, orderedNodes);
  loadMem(f, subLocId);
  loadMem(f, nodesOrdered);
  
  A.connected    = connected;
  A.orderedNodes = orderedNodes;
  A.subLocId     = subLocId;
  A.nodesOrdered = nodesOrdered;
}

