/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "sOctTreeItem.h"

sOctTreeItem::
sOctTreeItem()
{
  father = 0;
  isLeaf = true;
  
  for(UInt i=0; i<8; ++i)
  { sons[i] = 0; }
  
  id = 0;
}

sOctTreeItem::
sOctTreeItem(const UInt & Id)
{
  father = 0;
  isLeaf = true;
  
  for(UInt i=0; i<8; ++i)
  { sons[i] = 0; }
  
  id = Id;
}

sOctTreeItem::
sOctTreeItem(const sOctTreeItem & O)
{
  father = O.father;
  isLeaf = O.isLeaf;
  
  for(UInt i=0; i<8; ++i)
  { sons[i] = O.sons[i]; }
  
  id = O.id;
}

sOctTreeItem &
sOctTreeItem::
operator=(const sOctTreeItem & O)
{
  father = O.father;
  isLeaf = O.isLeaf;
  
  for(UInt i=0; i<8; ++i)
  { sons[i] = O.sons[i]; }
  
  id = O.id;
  
  return(*this);
}

UInt
sOctTreeItem::
octMap(const bool & ix, const bool & iy, const bool & iz)
{
  return(UInt(ix) + UInt(iy) * 2 + UInt(iz) * 4);
}

bool
sOctTreeItem::
operator<(const sOctTreeItem & V) const
{
  return(id < V.id);
}

bool
sOctTreeItem::
operator!=(const sOctTreeItem & V) const
{
  return(id != V.id);
}

void
sOctTreeItem::
setFather(const UInt & Father)
{
  father = Father;
}

void
sOctTreeItem::
setSons(const UInt VVV, const UInt FVV, const UInt VFV, const UInt FFV,
        const UInt VVF, const UInt FVF, const UInt VFF, const UInt FFF)
{
  sons[0] = VVV;
  sons[1] = FVV;
  sons[2] = VFV;
  sons[3] = FFV;
  
  sons[4] = VVF;
  sons[5] = FVF;
  sons[6] = VFF;
  sons[7] = FFF;
}

void
sOctTreeItem::
setSons(const sVect<UInt> Sons)
{
  assert(Sons.size() == 8);
  
  for(UInt i=0; i<8; ++i)
  { sons[i] = Sons(i+1); }
}

void
sOctTreeItem::
setIsLeaf(const bool & IsLeaf)
{
  isLeaf = IsLeaf;
}

void
sOctTreeItem::
setId(const UInt & Id)
{
  id = Id;
}

const UInt &
sOctTreeItem::
getFather() const
{
  return(father);
}

const UInt &
sOctTreeItem::
getSons(const bool & ix, const bool & iy, const bool & iz) const
{
  return(sons[ octMap(ix,iy,iz) ]);
}

const UInt &
sOctTreeItem::
getSons(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= 8);
  
  return(sons[i-1]);
}

const bool &
sOctTreeItem::
getIsLeaf() const
{
  return(isLeaf);
}

const UInt &
sOctTreeItem::
getId() const
{
  return(id);
}
    
UInt &
sOctTreeItem::
getFather()
{
  return(father);
}

UInt &
sOctTreeItem::
getSons(const bool & ix, const bool & iy, const bool & iz)
{
  return(sons[ octMap(ix,iy,iz) ]);
}

UInt &
sOctTreeItem::
getSons(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 8);
  
  return(sons[i-1]);
}

bool &
sOctTreeItem::
getIsLeaf()
{
  return(isLeaf);
}

UInt &
sOctTreeItem::
getId()
{
  return(id);
}
    
ostream & operator<<(ostream & f, const sOctTreeItem & P)
{
  f << "father : " << P.father << endl;
  f << "Id     : " << P.id << endl;
  
  f << "VVV : " << P.sons[0] << endl;
  f << "FVV : " << P.sons[1] << endl;
  f << "VFV : " << P.sons[2] << endl;
  f << "FFV : " << P.sons[3] << endl;
  
  f << "VVF : " << P.sons[4] << endl;
  f << "FVF : " << P.sons[5] << endl;
  f << "VFF : " << P.sons[6] << endl;
  f << "FFF : " << P.sons[7] << endl;
  
  f << "isLeaf : " << P.isLeaf << endl;
  
  return(f);
}

size_t
sOctTreeItem::
memSize() const
{
  return(sizeof(sOctTreeItem));
}
