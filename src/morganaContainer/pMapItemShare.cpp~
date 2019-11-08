/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "pMapItemShare.h"

pMapItemShare::
pMapItemShare() : pMapItem(), shared(false), owned(false)
{
}

pMapItemShare::
pMapItemShare(const UInt & Lid, const UInt & Gid) : pMapItem(Lid,Gid), shared(false), owned(false)
{
}

pMapItemShare::
pMapItemShare(const UInt & Lid, const UInt & Gid, const bool & Shared, const bool & Owned) : pMapItem(Lid,Gid), shared(Shared), owned(Owned)
{
}

pMapItemShare::
pMapItemShare(const UInt & Lid, const UInt & Gid, const UInt & Pid, const bool & Shared, const bool & Owned) : pMapItem(Lid,Gid,Pid), shared(Shared), owned(Owned)
{
}

pMapItemShare::
pMapItemShare(const pMapItemShare & M) : pMapItem(M), shared(M.shared), owned(M.owned)
{
}

pMapItemShare &
pMapItemShare::
operator=(const pMapItemShare & M)
{
  lid    = M.lid;
  gid    = M.gid;
  pid    = M.pid;
  shared = M.shared;
  owned  = M.owned;
  
  bufLid = M.bufLid;
  
  return(*this);
}

void
pMapItemShare::
setShared(const bool & Shared)
{
  shared = Shared;
}

void
pMapItemShare::
setOwned(const bool & Owned)
{
  owned = Owned;
}

bool &
pMapItemShare::
getShared()
{
  return(shared);
}

bool &
pMapItemShare::
getOwned()
{
  return(owned);
}

const bool &
pMapItemShare::
getShared() const
{
  return(shared);
}

const bool &
pMapItemShare::
getOwned() const
{
  return(owned);
}

ostream & operator<<(ostream & f, const pMapItemShare & M)
{
  f << "pid: " << M.getPid() << " lid: " << M.getLid() << " gid: " << M.getGid() << " shared: " << M.shared << " owned " << M.owned << endl;
  return(f);
}
