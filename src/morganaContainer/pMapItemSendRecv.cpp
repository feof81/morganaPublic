/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "pMapItemSendRecv.h"

pMapItemSendRecv::
pMapItemSendRecv() : pMapItem(), sid(0), rid(0)
{
}

pMapItemSendRecv::
pMapItemSendRecv(const UInt & Lid, const UInt & Gid) : pMapItem(Lid,Gid), sid(0), rid(0)
{
}

pMapItemSendRecv::
pMapItemSendRecv(const UInt & Lid, const UInt & Gid, const UInt & Sid, const UInt & Rid) : pMapItem(Lid,Gid), sid(Sid), rid(Rid)
{
}

pMapItemSendRecv::
pMapItemSendRecv(const pMapItemSendRecv & M) : pMapItem(M), sid(M.sid), rid(M.rid)
{
}

pMapItemSendRecv &
pMapItemSendRecv::
operator=(const pMapItemSendRecv & M)
{
  lid = M.lid;
  gid = M.gid;
  pid = M.pid;
  sid = M.sid;
  rid = M.rid;
  
  bufLid = M.bufLid;
  
  return(*this);
}

bool
pMapItemSendRecv::
operator<(const pMapItemSendRecv & E) const
{
  if(gid == E.gid)
  {
    if(sid == E.sid) { return(rid < E.rid); }
    return(sid < E.sid);
  }
  
  return(gid < E.gid);
}

bool
pMapItemSendRecv::
operator!=(const pMapItemSendRecv & E) const
{
  return((gid != E.gid) || (sid != E.sid) || (rid != E.rid));
}

void
pMapItemSendRecv::
setSid(const UInt & Sid)
{
  sid = Sid;
}

void
pMapItemSendRecv::
setRid(const UInt & Rid)
{
  rid = Rid;
}

UInt &
pMapItemSendRecv::
getSid()
{
  return(sid);
}

UInt &
pMapItemSendRecv::
getRid()
{
  return(rid);
}

const UInt &
pMapItemSendRecv::
getSid() const
{
  return(sid);
}

const UInt &
pMapItemSendRecv::
getRid() const
{
  return(rid);
}

ostream & operator<<(ostream & f, const pMapItemSendRecv & M)
{
  f << "pid: " << M.getPid() << " lid: " << M.getLid() << " gid: " << M.getGid() << " sid: " << M.getSid() << " rid " << M.getRid() << endl;
  return(f);
}

size_t
pMapItemSendRecv::
memSize() const
{
  return(sizeof(pMapItemSendRecv));
}


//_________________________________________________________________________________________________
// ASSOCIATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
size_t dynamicSizeOf(const pMapItemSendRecv & A)
{ return(A.memSize()); }
