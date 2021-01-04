/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "branch1d_card.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
branch1d_card::
branch1d_card() : nodeGid(0), intElementGid(0), extNodeGid(0)
{
}

branch1d_card::
branch1d_card(const branch1d_card & card)
{
  nodeGid        = card.nodeGid;
  intElementGid  = card.intElementGid;
  intNormal      = card.intNormal;
  extNodeGid     = card.extNodeGid;
  extElementsGid = card.extElementsGid;
  extNormals     = card.extNormals;
}

branch1d_card &
branch1d_card::
operator=(const branch1d_card & card)
{
  nodeGid        = card.nodeGid;
  intElementGid  = card.intElementGid;
  intNormal      = card.intNormal;
  extNodeGid     = card.extNodeGid;
  extElementsGid = card.extElementsGid;
  extNormals     = card.extNormals;
  
  return(*this);
}


//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
branch1d_card::
setNodeGid(const UInt & NodeGid)
{
  nodeGid = NodeGid;
}

void
branch1d_card::
setIntElementGid(const UInt & IntElementGid)
{
  intElementGid = IntElementGid;
}

void
branch1d_card::
setIntNormal(const point3d & IntNormal)
{
  intNormal = IntNormal;
}

void
branch1d_card::
setExtNodeGid(const UInt & ExtNodeGid)
{
  extNodeGid = ExtNodeGid;
}

void
branch1d_card::
setExtElementsGid(const sVect<UInt> & ExtElementsGid)
{
  extElementsGid = ExtElementsGid;
}

void
branch1d_card::
setExtNormals(const sVect<point3d> & ExtNormals)
{
  extNormals = ExtNormals;
}


//_________________________________________________________________________________________________
// ADD FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
branch1d_card::
addExtElementsGid(const UInt & gid)
{
  extElementsGid.push_back(gid);
}

void
branch1d_card::
addExtNormals(const point3d & P)
{
  extNormals.push_back(P);
}


//_________________________________________________________________________________________________
// GET FUNCTIONS - CONST
//-------------------------------------------------------------------------------------------------
const UInt &
branch1d_card::
getNodeGid() const
{
  return(nodeGid);
}

const UInt &
branch1d_card::
getIntElementGid() const
{
  return(intElementGid);
}

const point3d &
branch1d_card::
getIntNormal() const
{
  return(intNormal);
}

const UInt &
branch1d_card::
getExtNodeGid() const
{
  return(extNodeGid);
}
    
UInt
branch1d_card::
getExtElementsNum() const
{
  return(extElementsGid.size());
}

const UInt &
branch1d_card::
getExtElementsGid(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= extElementsGid.size());
    
  return(extElementsGid(i));
}

const sVect<UInt> &
branch1d_card::
getExtElementsGid() const
{
  return(extElementsGid);
}
    
UInt
branch1d_card::
getExtNormalsNum() const
{
  return(extNormals.size());
}

const point3d &
branch1d_card::
getExtNormals(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= extNormals.size());
    
  return(extNormals(i));
}

const sVect<point3d> &
branch1d_card::
getExtNormals() const
{
  return(extNormals);
}


//_________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
UInt &
branch1d_card::
getNodeGid()
{
  return(nodeGid);
}

UInt &
branch1d_card::
getIntElementGid()
{
  return(intElementGid);
}

point3d &
branch1d_card::
getIntNormal()
{
  return(intNormal);
}

UInt &
branch1d_card::
getExtNodeGid()
{
  return(extNodeGid);
}
    
UInt &
branch1d_card::
getExtElementsGid(const UInt & i)
{
  assert(i >= 1);
  assert(i <= extElementsGid.size());
    
  return(extElementsGid(i));
}

sVect<UInt> &
branch1d_card::
getExtElementsGid()
{
  return(extElementsGid);
}

point3d &
branch1d_card::
getExtNormals(const UInt & i)
{
  assert(i >= 1);
  assert(i <= extNormals.size());
    
  return(extNormals(i));
}

sVect<point3d> &
branch1d_card::
getExtNormals()
{
  return(extNormals);
}


//_________________________________________________________________________________________________
// OUTSTREAM OPERATOR
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const branch1d_card & card)
{
  f << "intNode gid      : " << card.nodeGid << endl;
  f << "intElement gid   : " << card.intElementGid << endl;
  f << "intNormal        : " << card.intNormal << endl;
  f << "extNode gid      : " << card.extNodeGid << endl;
  f << "extElements gids : " << card.extElementsGid << endl;
  f << "extNormals       : " << card.extNormals << endl;
  return(f);
}
