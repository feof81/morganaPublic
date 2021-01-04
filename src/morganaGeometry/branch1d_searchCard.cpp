/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "branch1d_searchCard.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
branch1d_searchCard::
branch1d_searchCard()
{
}

branch1d_searchCard::
branch1d_searchCard(const branch1d_searchCard & card)
{
  tgtNodeGid = card.tgtNodeGid;
  tgtP       = card.tgtP;
  tgtElGid   = card.tgtElGid;
  tgtNormal  = card.tgtNormal;
  
  srcNodeGid      = card.srcNodeGid;
  srcElementsGids = card.srcElementsGids;
  srcNormals      = card.srcNormals;
}

branch1d_searchCard &
branch1d_searchCard::
operator=(const branch1d_searchCard & card)
{
  tgtNodeGid = card.tgtNodeGid;
  tgtP       = card.tgtP;
  tgtElGid   = card.tgtElGid;
  tgtNormal  = card.tgtNormal;
  
  srcNodeGid      = card.srcNodeGid;
  srcElementsGids = card.srcElementsGids;
  srcNormals      = card.srcNormals;
  
  return(*this);
}


//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
branch1d_searchCard::
setTgtNodeGid(const UInt & TgtNodeGid)
{
  tgtNodeGid = TgtNodeGid;
}

void
branch1d_searchCard::
setTgtNode(const point3d & Q)
{
  tgtP = Q;
}

void
branch1d_searchCard::
setTgtElGid(const UInt & TgtElGid)
{
  tgtElGid = TgtElGid;
}

void
branch1d_searchCard::
setTgtNormal(const point3d & TgtNormal)
{
  tgtNormal = TgtNormal;
}

void
branch1d_searchCard::
setSrcNodeGid(const UInt & SrcNodeGid)
{
  srcNodeGid = SrcNodeGid;
}

void
branch1d_searchCard::
setSrcElementsGids(const sVect<UInt> & SrcElementsGids)
{
  srcElementsGids = SrcElementsGids;
}

void
branch1d_searchCard::
setSrcNormals(const sVect<point3d> & SrcNormals)
{
  srcNormals = SrcNormals;
}


//_________________________________________________________________________________________________
// ADD FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
branch1d_searchCard::
addSrcElementsGids(const UInt & gid)
{
  srcElementsGids.push_back(gid);
}

void
branch1d_searchCard::
addSrcNormals(const point3d & P)
{
  srcNormals.push_back(P);
}


//_________________________________________________________________________________________________
// GET FUNCTIONS - CONST
//-------------------------------------------------------------------------------------------------
const UInt &
branch1d_searchCard::
getTgtNodeGid() const
{
  return(tgtNodeGid);
}

const point3d &
branch1d_searchCard::
getTgtNode() const
{
  return(tgtP);
}

const UInt &
branch1d_searchCard::
getTgtElGid() const
{
  return(tgtElGid);
}

const point3d &
branch1d_searchCard::
getTgtNormal() const
{
  return(tgtNormal);
}

const UInt &
branch1d_searchCard::
getSrcNodeGid() const
{
  return(srcNodeGid);
}

UInt
branch1d_searchCard::
getSrcElementsNum() const
{
  return(srcElementsGids.size());
}

const UInt &
branch1d_searchCard::
getSrcElementsGids(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= srcElementsGids.size());
  
  return(srcElementsGids(i));
}

const sVect<UInt> &
branch1d_searchCard::
getSrcElementsGids() const
{
  return(srcElementsGids);
}


UInt
branch1d_searchCard::
getSrcNormalsNum() const
{
  return(srcNormals.size());
}

const point3d &
branch1d_searchCard::
getSrcNormals(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= srcNormals.size());
  
  return(srcNormals(i));
}

const sVect<point3d> &
branch1d_searchCard::
getSrcNormals() const
{
  return(srcNormals);
}


//_________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
UInt &
branch1d_searchCard::
getTgtNodeGid()
{
  return(tgtNodeGid);
}

point3d &
branch1d_searchCard::
getTgtNode()
{
  return(tgtP);
}

UInt &
branch1d_searchCard::
getTgtElGid()
{
  return(tgtElGid);
}

point3d &
branch1d_searchCard::
getTgtNormal()
{
  return(tgtNormal);
}

UInt &
branch1d_searchCard::
getSrcNodeGid()
{
  return(srcNodeGid);
}

sVect<UInt> &
branch1d_searchCard::
getSrcElementsGids()
{
  return(srcElementsGids);
}

UInt &
branch1d_searchCard::
getSrcElementsGids(const UInt & i)
{
  assert(i >= 1);
  assert(i <= srcElementsGids.size());
    
  return(srcElementsGids(i));
}

sVect<point3d> &
branch1d_searchCard::
getSrcNormals()
{
  return(srcNormals);
}

point3d &
branch1d_searchCard::
getSrcNormals(const UInt & i)
{
  assert(i >= 1);
  assert(i <= srcNormals.size());
    
  return(srcNormals(i));
}


//_________________________________________________________________________________________________
// OUTSTREAM OPERATOR
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const branch1d_searchCard & card)
{
  f << "tgtNode   : " << card.tgtNodeGid << endl;
  f << "tgtElGid  : " << card.tgtElGid   << endl;
  f << "tgtNormal : " << card.tgtNormal  << endl;
  f << "srcNode   : " << card.srcNodeGid << endl;
  f << "srcElementsGids : " << card.srcElementsGids;
  f << "srcNormals      : " << card.srcNormals << endl;
  return(f);
}
