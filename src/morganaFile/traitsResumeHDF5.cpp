/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "traitsResumeHDF5.h"


//_________________________________________________________________________________________________
// PMAPITEM
//-------------------------------------------------------------------------------------------------
UInt
traitsResumeHDF5<pMapItem>::
size()
{
  return(4);
}

traitsResumeHDF5<pMapItem>::OUTTYPE
traitsResumeHDF5<pMapItem>::
getValue(const DATA & data, const UInt & i)
{
  assert(i >= 1);
  assert(i <= 4);
  
  switch(i)
  {
    case 1 : return(data.getLid());
    case 2 : return(data.getGid());
    case 3 : return(data.getPid());
    case 4 : return(data.getBufLid());
  }
  
  return(0);
}

void
traitsResumeHDF5<pMapItem>::
setValue(const OUTTYPE & item, const UInt & i, DATA & data)
{
  assert(i >= 1);
  assert(i <= 4);
  
  switch(i)
  {
    case 1 : data.setLid(item);    break;
    case 2 : data.setGid(item);    break;
    case 3 : data.setPid(item);    break;
    case 4 : data.setBufLid(item); break;
  }
}


//_________________________________________________________________________________________________
// PMAPITEMSHARE
//-------------------------------------------------------------------------------------------------
UInt
traitsResumeHDF5<pMapItemShare>::
size()
{
  return(6);
}

traitsResumeHDF5<pMapItemShare>::OUTTYPE
traitsResumeHDF5<pMapItemShare>::
getValue(const DATA & data, const UInt & i)
{
  assert(i >= 1);
  assert(i <= 6);
  
  switch(i)
  {
    case 1 : return(data.getLid());
    case 2 : return(data.getGid());
    case 3 : return(data.getPid());
    case 4 : return(UInt(data.getShared()));
    case 5 : return(UInt(data.getOwned()));
    case 6 : return(data.getBufLid());
  }
  
  return(0);
}

void
traitsResumeHDF5<pMapItemShare>::
setValue(const OUTTYPE & item, const UInt & i, DATA & data)
{
  assert(i >= 1);
  assert(i <= 6);
  
  switch(i)
  {
    case 1 : data.setLid(item);          break;
    case 2 : data.setGid(item);          break;
    case 3 : data.setPid(item);          break;
    case 4 : data.setShared(bool(item)); break;
    case 5 : data.setOwned(bool(item));  break;
    case 6 : data.setBufLid(item);       break;
  }
}


//_________________________________________________________________________________________________
// PMAPITEMSENDRECV
//-------------------------------------------------------------------------------------------------
UInt
traitsResumeHDF5<pMapItemSendRecv>::
size()
{
  return(6);
}

traitsResumeHDF5<pMapItemSendRecv>::OUTTYPE
traitsResumeHDF5<pMapItemSendRecv>::
getValue(const DATA & data, const UInt & i)
{
  assert(i >= 1);
  assert(i <= 6);
  
  switch(i)
  {
    case 1 : return(data.getLid());
    case 2 : return(data.getGid());
    case 3 : return(data.getPid());
    case 4 : return(data.getSid());
    case 5 : return(data.getRid());
    case 6 : return(data.getBufLid());
  }
  
  return(0);
}

void
traitsResumeHDF5<pMapItemSendRecv>::
setValue(const OUTTYPE & item, const UInt & i, DATA & data)
{
  assert(i >= 1);
  assert(i <= 6);
  
  switch(i)
  {
    case 1 : data.setLid(item);    break;
    case 2 : data.setGid(item);    break;
    case 3 : data.setPid(item);    break;
    case 4 : data.setSid(item);    break;
    case 5 : data.setRid(item);    break;
    case 6 : data.setBufLid(item); break;
  }
}


//_________________________________________________________________________________________________
// UInt
//-------------------------------------------------------------------------------------------------
UInt
traitsResumeHDF5<UInt>::
size()
{
  return(1);
}

traitsResumeHDF5<UInt>::OUTTYPE
traitsResumeHDF5<UInt>::
getValue(const DATA & data, const UInt & i)
{
  assert(i == 1); 
  return(data);
}

void
traitsResumeHDF5<UInt>::
setValue(const OUTTYPE & item, const UInt & i, DATA & data)
{
  assert(i == 1);
  data = item;
}



//_________________________________________________________________________________________________
// REAL
//-------------------------------------------------------------------------------------------------
UInt
traitsResumeHDF5<Real>::
size()
{
  return(1);
}

traitsResumeHDF5<Real>::OUTTYPE
traitsResumeHDF5<Real>::
getValue(const DATA & data, const UInt & i)
{
  assert(i == 1); 
  return(data);
}

void
traitsResumeHDF5<Real>::
setValue(const OUTTYPE & item, const UInt & i, DATA & data)
{
  assert(i == 1);
  data = item;
}


//_________________________________________________________________________________________________
// POINT2D
//-------------------------------------------------------------------------------------------------
UInt
traitsResumeHDF5<point2d>::
size()
{
  return(2);
}

traitsResumeHDF5<point2d>::OUTTYPE
traitsResumeHDF5<point2d>::
getValue(const DATA & data, const UInt & i)
{
  assert(i >= 1);
  assert(i <= 2);
  
  return(data.getI(i));
}

void
traitsResumeHDF5<point2d>::
setValue(const OUTTYPE & item, const UInt & i, DATA & data)
{
  assert(i >= 1);
  assert(i <= 2);
  
  data.setI(i,item);
}


//_________________________________________________________________________________________________
// POINT3D
//-------------------------------------------------------------------------------------------------
UInt
traitsResumeHDF5<point3d>::
size()
{
  return(3);
}

traitsResumeHDF5<point3d>::OUTTYPE
traitsResumeHDF5<point3d>::
getValue(const DATA & data, const UInt & i)
{
  assert(i >= 1);
  assert(i <= 3);
  
  return(data.getI(i));
}

void
traitsResumeHDF5<point3d>::
setValue(const OUTTYPE & item, const UInt & i, DATA & data)
{
  assert(i >= 1);
  assert(i <= 3);
  
  data.setI(i,item);
}


//_________________________________________________________________________________________________
// TENSOR2D
//-------------------------------------------------------------------------------------------------
UInt
traitsResumeHDF5<tensor2d>::
size()
{
  return(4);
}

traitsResumeHDF5<tensor2d>::OUTTYPE
traitsResumeHDF5<tensor2d>::
getValue(const DATA & data, const UInt & i)
{
  assert(i >= 1);
  assert(i <= 4);
  
  switch(i)
  {
    case 1 : return(data.getIJ(1,1));
    case 2 : return(data.getIJ(1,2));
    case 3 : return(data.getIJ(2,1));
    case 4 : return(data.getIJ(2,2));
  }
  
  return(0);
}

void
traitsResumeHDF5<tensor2d>::
setValue(const OUTTYPE & item, const UInt & i, DATA & data)
{
  assert(i >= 1);
  assert(i <= 4);
  
  switch(i)
  {
    case 1 : data.setIJ(1,1,item); break;
    case 2 : data.setIJ(1,2,item); break;
    case 3 : data.setIJ(2,1,item); break;
    case 4 : data.setIJ(2,2,item); break;
  }
}


//_________________________________________________________________________________________________
// TENSOR2D
//-------------------------------------------------------------------------------------------------
UInt
traitsResumeHDF5<tensor3d>::
size()
{
  return(9);
}

traitsResumeHDF5<tensor3d>::OUTTYPE
traitsResumeHDF5<tensor3d>::
getValue(const DATA & data, const UInt & i)
{
  assert(i >= 1);
  assert(i <= 9);
  
  switch(i)
  {
    case 1 : return(data.getIJ(1,1));
    case 2 : return(data.getIJ(1,2));
    case 3 : return(data.getIJ(1,3));
    
    case 4 : return(data.getIJ(2,1));
    case 5 : return(data.getIJ(2,2));
    case 6 : return(data.getIJ(2,3));
    
    case 7 : return(data.getIJ(3,1));
    case 8 : return(data.getIJ(3,2));
    case 9 : return(data.getIJ(3,3));
  }
  
  return(0);
}

void
traitsResumeHDF5<tensor3d>::
setValue(const OUTTYPE & item, const UInt & i, DATA & data)
{
  assert(i >= 1);
  assert(i <= 9);
  
  switch(i)
  {
    case 1 : data.setIJ(1,1,item); break;
    case 2 : data.setIJ(1,2,item); break;
    case 3 : data.setIJ(1,3,item); break;
    
    case 4 : data.setIJ(2,1,item); break;
    case 5 : data.setIJ(2,2,item); break;
    case 6 : data.setIJ(2,3,item); break;
    
    case 7 : data.setIJ(3,1,item); break;
    case 8 : data.setIJ(3,2,item); break;
    case 9 : data.setIJ(3,3,item); break;
  }
}

