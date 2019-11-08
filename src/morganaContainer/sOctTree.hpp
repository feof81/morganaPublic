/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef SOCTTREE_H
#define SOCTTREE_H

#include "sOctTreeItem.h"

template<typename DATA>
class sOctTree
{
    /*! @name Typedefs */ //@{
  public:
    typedef std::map<sOctTreeItem,DATA>      DATAMAP;
    typedef typename DATAMAP::iterator       ITER;
    typedef typename DATAMAP::const_iterator CONST_ITER;
    //@}
  
    /*! @name Internal data */ //@{
  public:
    UInt k;
    ITER iter;
    DATAMAP data;
    //@}
    
    /*! @name Constructors and functions */ //@{
  public:
    sOctTree();
    sOctTree(const sOctTree & O);
    sOctTree & operator=(const sOctTree & O);
    void clear();
    //@}
    
    /*! @name Tree handle */ //@{
  public:
    void restart();
    void goDown(const UInt & i);
    void goDown(const bool & ix, const bool & iy, const bool & iz);
    void goUp();
    void goTo(const UInt & Id);
    //@}
    
    /*! @name Set functions */ //@{
  public:
    void setData(const DATA & Data);
    void setData(const UInt & Id, const DATA & Data);
    sVect<UInt> addLeaves();
    sVect<UInt> addLeaves(const UInt & Id);
    sVect<UInt> addLeaves(const DATA & VVV, const DATA & FVV, const DATA & VFV, const DATA & FFV,
                          const DATA & VVF, const DATA & FVF, const DATA & VFF, const DATA & FFF);
    sVect<UInt> addLeaves(const UInt & Id,
                          const DATA & VVV, const DATA & FVV, const DATA & VFV, const DATA & FFV,
                          const DATA & VVF, const DATA & FVF, const DATA & VFF, const DATA & FFF);
    //@}
    
    /*! @name Const get functions */ //@{
  public:
    const sOctTreeItem & getMap() const;
    const sOctTreeItem & getMap(const UInt & Id) const;
    DATA       & getData();
    const DATA & getData() const;
    DATA       & getData(const UInt & Id);
    const DATA & getData(const UInt & Id) const;
    //@}
    
    /*! @name Outstream operator */ //@{
  public:
    template<typename D>
    friend ostream & operator<<(ostream & f, const sOctTree<D> & O);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename DATA>
sOctTree<DATA>::
sOctTree()
{
  k = 1;
}

template<typename DATA>
sOctTree<DATA>::
sOctTree(const sOctTree & O)
{
  k = O.iter;
  iter  = O.iter;
  data  = O.data;
}

template<typename DATA>
sOctTree<DATA> &
sOctTree<DATA>::
operator=(const sOctTree & O)
{
  k = O.k;
  iter  = O.iter;
  data  = O.data;
  
  return(*this);
}

template<typename DATA>
void
sOctTree<DATA>::
clear()
{
  k = 1;
  data.clear();
}


//_________________________________________________________________________________________________
// TREE HANDLE
//-------------------------------------------------------------------------------------------------
template<typename DATA>
void
sOctTree<DATA>::
restart()
{
  assert(data.size() != 0);
  iter = data.begin();
}

template<typename DATA>
void
sOctTree<DATA>::
goDown(const UInt & i)
{
  assert(data.size() != 0);
  
  UInt id = iter->first.getSons(i);
  iter = data.find(sOctTreeItem(id));
}

template<typename DATA>
void
sOctTree<DATA>::
goDown(const bool & ix, const bool & iy, const bool & iz)
{
  assert(data.size() != 0);
  
  UInt i  = sOctTreeItem::octMap(ix,iy,iz) + 1;
  UInt id = iter->first.getSons(i);
  iter = data.find(sOctTreeItem(id));
}

template<typename DATA>
void
sOctTree<DATA>::
goUp()
{
  assert(data.size() != 0);
  
  UInt id = iter->first.getFather();
  iter = data.find(sOctTreeItem(id));
}

template<typename DATA>
void
sOctTree<DATA>::
goTo(const UInt & Id)
{
  assert(data.size() != 0);
  
  iter = data.find(sOctTreeItem(Id));
}
    
    

//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DATA>
void
sOctTree<DATA>::
setData(const DATA & Data)
{
  iter->second = Data;
}

template<typename DATA>
void
sOctTree<DATA>::
setData(const UInt & Id, const DATA & Data)
{
  iter = data.find(sOctTreeItem(Id));
  iter->second = Data;
}

template<typename DATA>
sVect<UInt>
sOctTree<DATA>::
addLeaves()
{
  DATA empty;
  sOctTreeItem item;
  std::pair<sOctTreeItem,DATA> paio;
  paio.second = empty;
  sVect<UInt> out;
  
  if(data.size() == 0)
  {
    item.setFather(0);
    item.setId(k);
    item.setIsLeaf(true);
    
    paio.first = item;
    data.insert(paio);
    
    iter = data.begin();
    k++;
  }
  else
  {
    assert(iter->first.getIsLeaf());
    
    //Modify the father
    UInt id = iter->first.getId();
    
    iter->first.setSons(k,   k+1, k+2, k+3,
                        k+4, k+5, k+6, k+7);
    iter->first.setIsLeaf(false);
    
    out.resize(8);
    for(UInt i=1; i <=8; ++i)
    { out(i) = k +i -1 ; }
    
    //Add sons
    item.setIsLeaf(true);
    item.setFather(id);
    
    item.setId(k);
    paio.first = item;
    data.insert(paio);
    
    item.setId(k+1);
    paio.first = item;
    data.insert(paio);
    
    item.setId(k+2);
    paio.first = item;
    data.insert(paio);
    
    item.setId(k+3);
    paio.first = item;
    data.insert(paio);
    
    
    item.setId(k+4);
    paio.first = item;
    data.insert(paio);
    
    item.setId(k+5);
    paio.first = item;
    data.insert(paio);
    
    item.setId(k+6);
    paio.first = item;
    data.insert(paio);
    
    item.setId(k+7);
    paio.first = item;
    data.insert(paio);
    
    k += 8;
  }
  
  return(out);
}

template<typename DATA>
sVect<UInt>
sOctTree<DATA>::
addLeaves(const UInt & Id)
{
  sOctTreeItem item;
  sVect<UInt> out;
  
  if(data.size() == 0)
  {
    item.setFather(0);
    item.setId(k);
    item.setIsLeaf(true);
    
    data.insert(item);
    iter = data.begin();
    
    k++;
  }
  else
  {
    iter = data.find(sOctTreeItem(Id));
    
    //Modify the father
    UInt father = iter->first.getFather();
    UInt id     = iter->first.getId();
    
    iter->first.setSons(k,   k+1, k+2, k+3,
                        k+4, k+5, k+6, k+7);
    iter->first.setIsLeaf(false);
    
    out.resize(8);
    for(UInt i=1; i <=8; ++i)
    { out(i) = k +i -1 ; }
    
    //Add sons
    item.setIsLeaf(true);
    item.setFather(id);
    
    
    item.setId(k);
    data.insert(item);
    
    item.setId(k+1);
    data.insert(item);
    
    item.setId(k+2);
    data.insert(item);
    
    item.setId(k+3);
    data.insert(item);
    
    
    item.setId(k+4);
    data.insert(item);
    
    item.setId(k+5);
    data.insert(item);
    
    item.setId(k+6);
    data.insert(item);
    
    item.setId(k+7);
    data.insert(item);
    
    k += 8;
  }
  
  return(out);
}

template<typename DATA>
sVect<UInt>
sOctTree<DATA>::
addLeaves(const DATA & VVV, const DATA & FVV, const DATA & VFV, const DATA & FFV,
          const DATA & VVF, const DATA & FVF, const DATA & VFF, const DATA & FFF)
{
  DATA empty;
  sOctTreeItem item;
  std::pair<sOctTreeItem,DATA> paio;
  paio.second = empty;
  sVect<UInt> out;
  
  if(data.size() == 0)
  {
    item.setFather(0);
    item.setId(k);
    item.setIsLeaf(true);
    
    paio.first = item;
    data.insert(paio);
    
    iter = data.begin();
    k++;
  }
  else
  {
    assert(iter->first.getIsLeaf());
    
    //Modify the father
    UInt id = iter->first.getId();
    
    iter->first.setSons(k,   k+1, k+2, k+3,
                        k+4, k+5, k+6, k+7);
    iter->first.setIsLeaf(false);
    
    out.resize(8);
    for(UInt i=1; i <=8; ++i)
    { out(i) = k +i -1 ; }
    
    //Add sons
    item.setIsLeaf(true);
    item.setFather(id);
    
    item.setId(k);
    paio.first = item;
    paio.second = VVV;
    data.insert(paio);
    
    item.setId(k+1);
    paio.first = item;
    paio.second = FVV;
    data.insert(paio);
    
    item.setId(k+2);
    paio.first = item;
    paio.second = VFV;
    data.insert(paio);
    
    item.setId(k+3);
    paio.first = item;
    paio.second = FFV;
    data.insert(paio);
    
    
    item.setId(k+4);
    paio.first = item;
    paio.second = VVF;
    data.insert(paio);
    
    item.setId(k+5);
    paio.first = item;
    paio.second = FVF;
    data.insert(paio);
    
    item.setId(k+6);
    paio.first = item;
    paio.second = VFF;
    data.insert(paio);
    
    item.setId(k+7);
    paio.first = item;
    paio.second = FFF;
    data.insert(paio);
    
    k += 8;
  }
  
  return(out);
}


template<typename DATA>
sVect<UInt>
sOctTree<DATA>::
addLeaves(const UInt & Id,
          const DATA & VVV, const DATA & FVV, const DATA & VFV, const DATA & FFV,
          const DATA & VVF, const DATA & FVF, const DATA & VFF, const DATA & FFF)
{
  DATA empty;
  sOctTreeItem item;
  std::pair<sOctTreeItem,DATA> paio;
  paio.second = empty;
  sVect<UInt> out;
  
  if(data.size() == 0)
  {
    item.setFather(0);
    item.setId(k);
    item.setIsLeaf(true);
    
    paio.first = item;
    data.insert(paio);
    
    iter = data.begin();
    k++;
  }
  else
  {
    iter = data.find(sOctTreeItem(Id));
    
    //Modify the father
    UInt father = iter->first.getFather();
    UInt id     = iter->first.getId();
    
    iter->first.setSons(k,   k+1, k+2, k+3,
                        k+4, k+5, k+6, k+7);
    iter->first.setIsLeaf(false);
    
    out.resize(8);
    for(UInt i=1; i <=8; ++i)
    { out(i) = k +i -1 ; }
    
    //Add sons
    item.setIsLeaf(true);
    item.setFather(id);
    
    item.setId(k);
    paio.first = item;
    paio.second = VVV;
    data.insert(paio);
    
    item.setId(k+1);
    paio.first = item;
    paio.second = FVV;
    data.insert(paio);
    
    item.setId(k+2);
    paio.first = item;
    paio.second = VFV;
    data.insert(paio);
    
    item.setId(k+3);
    paio.first = item;
    paio.second = FFV;
    data.insert(paio);
    
    
    item.setId(k+4);
    paio.first = item;
    paio.second = VVF;
    data.insert(paio);
    
    item.setId(k+5);
    paio.first = item;
    paio.second = FVF;
    data.insert(paio);
    
    item.setId(k+6);
    paio.first = item;
    paio.second = VFF;
    data.insert(paio);
    
    item.setId(k+7);
    paio.first = item;
    paio.second = FFF;
    data.insert(paio);
    
    k += 8;
  }
  
  return(out);
}



//_________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DATA>
const sOctTreeItem &
sOctTree<DATA>::
getMap() const
{
  assert(data.size() != 0);
  return(iter->first);
}

template<typename DATA>
const sOctTreeItem &
sOctTree<DATA>::
getMap(const UInt & Id) const
{
  assert(data.size() != 0);
  assert(data.count(sOctTreeItem(Id)) >= 0);
  
  CONST_ITER tempIter = data.find(sOctTreeItem(Id));
  return(tempIter->first);
}

template<typename DATA>
DATA &
sOctTree<DATA>::
getData()
{
  assert(data.size() != 0);
  return(iter->second);
}

template<typename DATA>
const DATA &
sOctTree<DATA>::
getData() const
{
  assert(data.size() != 0);
  return(iter->second);
}

template<typename DATA>
DATA &
sOctTree<DATA>::
getData(const UInt & Id)
{
  assert(data.size() != 0);
  assert(data.count(sOctTreeItem(Id)) >= 0);
  
  ITER tempIter = data.find(sOctTreeItem(Id));
  return(tempIter->second);
}

template<typename DATA>
const DATA &
sOctTree<DATA>::
getData(const UInt & Id) const
{
  assert(data.size() != 0);
  assert(data.count(sOctTreeItem(Id)) >= 0);
  
  ITER tempIter = data.find(sOctTreeItem(Id));
  return(tempIter->second);
}
  

//_________________________________________________________________________________________________
// PRINTOUT
//-------------------------------------------------------------------------------------------------
template<typename D>
ostream & operator<<(ostream & f, const sOctTree<D> & O)
{
  typedef std::map<sOctTreeItem,D>      DATAMAP;
  typedef typename DATAMAP::const_iterator ITER;
  
  for(ITER tempIt = O.data.begin(); tempIt != O.data.end(); ++tempIt)
  {
    f << "Map - - - - - - -" << endl;
    f << tempIt->first << endl;
    f << "Data" << endl;
    f << tempIt->second << endl << endl;
  }
  
  return(f);
}


#endif
