/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "branch1d_bif.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
branch1d_bif::
branch1d_bif() : nodeGid(0), extNodeGid(0), extElementGid(0)
{
}

branch1d_bif::
branch1d_bif(const branch1d_bif & card)
{
  nodeGid        = card.nodeGid;
  intElementsGid = card.intElementsGid;
  intNormals     = card.intNormals;
  extNodeGid     = card.extNodeGid;
  extElementGid  = card.extElementGid;
  extNormal      = card.extNormal;
}

branch1d_bif &
branch1d_bif::
operator=(const branch1d_bif & card)
{
  nodeGid        = card.nodeGid;
  intElementsGid = card.intElementsGid;
  intNormals     = card.intNormals;
  extNodeGid     = card.extNodeGid;
  extElementGid  = card.extElementGid;
  extNormal      = card.extNormal;
  
  return(*this);
}


//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
branch1d_bif::
setNodeGid(const UInt & NodeGid)
{
  nodeGid = NodeGid;
}

void
branch1d_bif::
setIntElementsGid(const sVect<UInt> & IntElementsGid)
{
  intElementsGid = IntElementsGid;
}

void
branch1d_bif::
setIntNormal(const sVect<point3d> & IntNormal)
{
  intNormals = IntNormal;
}

void
branch1d_bif::
setExtNodeGid(const UInt & ExtNodeGid)
{
  extNodeGid = ExtNodeGid;
}

void
branch1d_bif::
setExtElementGid(const UInt & ExtElementGid)
{
  extElementGid = ExtElementGid;
}

void
branch1d_bif::
setExtNormal(const point3d & ExtNormal)
{
  extNormal = ExtNormal;
}


//_________________________________________________________________________________________________
// ADD FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
branch1d_bif::
addIntElementsGid(const UInt & gid)
{
  intElementsGid.push_back(gid);
}

void
branch1d_bif::
addIntNormals(const point3d & P)
{
  intNormals.push_back(P);
}


//_________________________________________________________________________________________________
// GET FUNCTIONS - CONST
//-------------------------------------------------------------------------------------------------
const UInt &
branch1d_bif::
getNodeGid() const
{
  return(nodeGid);
}

UInt
branch1d_bif::
getIntElementsGidNum() const
{
  return(intElementsGid.size());
}

const sVect<UInt> &
branch1d_bif::
getIntElementsGid() const
{
  return(intElementsGid);
}

UInt
branch1d_bif::
getIntNormalsNum() const
{
  return(intNormals.size());
}

const sVect<point3d> &
branch1d_bif::
getIntNormals() const
{
  return(intNormals);
}

const UInt &
branch1d_bif::
getExtNodeGid() const
{
  return(extNodeGid);
}

const UInt &
branch1d_bif::
getExtElementGid() const
{
  return(extElementGid);
}

const point3d &
branch1d_bif::
getExtNormal() const
{
  return(extNormal);
}


//_________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
UInt &
branch1d_bif::
getNodeGid()
{
  return(nodeGid);
}

sVect<UInt> &
branch1d_bif::
getIntElementsGid()
{
  return(intElementsGid);
}

sVect<point3d> &
branch1d_bif::
getIntNormals()
{
  return(intNormals);
}

UInt &
branch1d_bif::
getExtNodeGid()
{
  return(extNodeGid);
}

UInt &
branch1d_bif::
getExtElementGid()
{
  return(extElementGid);
}

point3d &
branch1d_bif::
getExtNormal()
{
  return(extNormal);
}


//_________________________________________________________________________________________________
// OUTSTREAM OPERATOR
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const branch1d_bif & card)
{
  f << "intNode gid     : " << card.nodeGid << endl;
  f << "intElements gid : " << card.intElementsGid << endl;
  f << "intNormals      : " << card.intNormals << endl;
  f << "extNode gid     : " << card.extNodeGid << endl;
  f << "extElements gid : " << card.extElementGid << endl;
  f << "extNormal       : " << card.extNormal << endl;
  
  return(f);
}
