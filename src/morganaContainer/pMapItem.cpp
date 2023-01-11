/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "pMapItem.h"

pMapItem::
pMapItem() : lid(0),
             gid(0),
             pid(0),
             bufLid(0)
{
}

pMapItem::
pMapItem(const UInt & Lid, const UInt & Gid) : lid(Lid),
                                               gid(Gid),
                                               pid(0),
                                               bufLid(0)
{
}

pMapItem::
pMapItem(const UInt & Lid, const UInt & Gid, const UInt & Pid) : lid(Lid),
                                                                 gid(Gid),
                                                                 pid(Pid),
                                                                 bufLid(0)
{
}

pMapItem::
pMapItem(const pMapItem & M) : lid(M.lid),
                               gid(M.gid),
                               pid(M.pid),
                               bufLid(M.bufLid)
{
}

pMapItem &
pMapItem::
operator=(const pMapItem & M)
{
  lid = M.lid;
  gid = M.gid;
  pid = M.pid;
  
  bufLid = M.bufLid;
  
  return(*this);
}

bool
pMapItem::
operator<(const pMapItem & E) const
{
  return(gid < E.gid);
}

bool
pMapItem::
operator!=(const pMapItem & E) const
{
  return(gid != E.gid);
}
    
void
pMapItem::
setLid(const UInt & Lid)
{
  lid = Lid;
}

void
pMapItem::
setGid(const UInt & Gid)
{
  gid = Gid;
}

void
pMapItem::
setPid(const UInt & Pid)
{
  pid = Pid;
}

void
pMapItem::
setBufLid(const UInt & BufLid)
{
  bufLid = BufLid;
}

UInt &
pMapItem::
getLid()
{
  return(lid);
}

UInt &
pMapItem::
getGid()
{
  return(gid);
}

UInt &
pMapItem::
getPid()
{
  return(pid);
}

UInt &
pMapItem::
getBufLid()
{
  return(bufLid);
}

const UInt &
pMapItem::
getLid() const
{
  return(lid);
}

const UInt &
pMapItem::
getGid() const
{
  return(gid);
}

const UInt &
pMapItem::
getPid() const
{
  return(pid);
}

const UInt &
pMapItem::
getBufLid() const
{
  return(bufLid);
}

void
pMapItem::
bufferLid()
{
  bufLid = lid;
}

void
pMapItem::
restoreLid()
{
  lid = bufLid;
}
    
ostream & operator<<(ostream & f, const pMapItem & M)
{
  f << "pid: " << M.getPid() << " lid: " << M.getLid() << " gid: " << M.getGid() << endl;
  return(f);
}

size_t
pMapItem::
memSize() const
{
  return(sizeof(pMapItem));
}


//_________________________________________________________________________________________________
// ASSOCIATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
size_t dynamicSizeOf(const pMapItem & A)
{ return(A.memSize()); }
