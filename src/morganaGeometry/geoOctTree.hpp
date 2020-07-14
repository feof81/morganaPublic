/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef GEOOCTTREE_HPP
#define GEOOCTTREE_HPP

#include "Teuchos_RCPDecl.hpp"

#include "morganaTypes.hpp"
#include "searchData.hpp"
#include "sOctTree.hpp"
#include "searchBoundingBox.h"

template<typename MESH, typename CONNECT>
class geoOctTree
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename MESH::GRID_GEOSHAPE GEOSHAPE;
    typedef typename MESH::GRID_ELMAP    ELMAP;
    typedef searchData<ELMAP>            SEARCHDATA;
    typedef sOctTree<searchBoundingBox>  OCTTREE;
    //@}
    
    /*! @name Internal flags */ //@{
  public:
    bool startupLocal, meshLoaded;
    OCTTREE octTree;
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<const MESH>    grid;
    Teuchos::RCP<const CONNECT> connectGrid;
    //@}
    
    /*! @name Constructor and setting functions */ //@{
  public:
    geoOctTree();
    geoOctTree(const Teuchos::RCP<const MESH> & Grid, const Teuchos::RCP<const CONNECT> & ConnectGrid);
    void setMesh(const Teuchos::RCP<const MESH> & Grid, const Teuchos::RCP<const CONNECT> & ConnectGrid);
    void setMesh(const MESH & Grid, const CONNECT & ConnectGrid);
    //@}
    
    /*! @name Startup and search function */ //@{
  public:
    void localInit(const Real & OctToll = 2.0);
    const sVect<UInt> getMatchingElements(const point3d & P);
    const sVect<UInt> getMatchingElementsExt(const point3d & P);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
template<typename MESH, typename CONNECT>
geoOctTree<MESH,CONNECT>::
geoOctTree()
{
  typedef typename MESH::GRID_ELMAP    MESH_ELMAP;
  typedef typename MESH::GRID_NODEMAP  MESH_NODEMAP;
  typedef typename MESH::GRID_GEOSHAPE MESH_GEOSHAPE;
  
  typedef typename CONNECT::GRID_ELMAP    CONNECT_ELMAP;
  typedef typename CONNECT::GRID_NODEMAP  CONNECT_NODEMAP;
  typedef typename CONNECT::GRID_GEOSHAPE CONNECT_GEOSHAPE;
  
  assert(staticAssert<MESH_GEOSHAPE::geoName     == CONNECT_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<MESH_ELMAP::parallelType   == CONNECT_ELMAP::parallelType>::returnValue);
  assert(staticAssert<MESH_NODEMAP::parallelType == CONNECT_NODEMAP::parallelType>::returnValue);
  
  startupLocal = false;
  meshLoaded   = false;
}

template<typename MESH, typename CONNECT>
geoOctTree<MESH,CONNECT>::
geoOctTree(const Teuchos::RCP<const MESH> & Grid, const Teuchos::RCP<const CONNECT> & ConnectGrid)
{
  typedef typename MESH::GRID_ELMAP    MESH_ELMAP;
  typedef typename MESH::GRID_NODEMAP  MESH_NODEMAP;
  typedef typename MESH::GRID_GEOSHAPE MESH_GEOSHAPE;
  
  typedef typename CONNECT::GRID_ELMAP    CONNECT_ELMAP;
  typedef typename CONNECT::GRID_NODEMAP  CONNECT_NODEMAP;
  typedef typename CONNECT::GRID_GEOSHAPE CONNECT_GEOSHAPE;
  
  assert(staticAssert<MESH_GEOSHAPE::geoName     == CONNECT_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<MESH_ELMAP::parallelType   == CONNECT_ELMAP::parallelType>::returnValue);
  assert(staticAssert<MESH_NODEMAP::parallelType == CONNECT_NODEMAP::parallelType>::returnValue);
  
  grid        = Grid;
  connectGrid = ConnectGrid;
  
  startupLocal = false;
  meshLoaded   = true;
}

template<typename MESH, typename CONNECT>
void
geoOctTree<MESH,CONNECT>::
setMesh(const Teuchos::RCP<const MESH> & Grid, const Teuchos::RCP<const CONNECT> & ConnectGrid)
{
  grid        = Grid;
  connectGrid = ConnectGrid;
  
  meshLoaded   = true;
  startupLocal = false;
}

template<typename MESH, typename CONNECT>
void
geoOctTree<MESH,CONNECT>::
setMesh(const MESH & Grid, const CONNECT & ConnectGrid)
{
  grid        = Teuchos::rcp(new MESH(Grid));
  connectGrid = Teuchos::rcp(new CONNECT(ConnectGrid));
  
  meshLoaded   = true;
  startupLocal = false;
}



//_________________________________________________________________________________________________
// STARTUP FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename MESH, typename CONNECT>
void
geoOctTree<MESH,CONNECT>::
localInit(const Real & OctToll)
{
  //! Assert-------------------------------------------------------------------
  assert(meshLoaded);
  startupLocal = true;
  
  //! Building bounding boxes--------------------------------------------------
  sVect<point3d> upperBound(grid->getNumElements()), lowerBound(grid->getNumElements());
  sVect<point3d> elNodes;
  point3d Pmin, Pmax;
  point3d bbMax, bbMin;
  
  if(grid->getNumElements() != 0)
  {
    elNodes = grid->getElementNodesL(1);
    bbMax   = elNodes(1);
    bbMin   = elNodes(1);
  }
  
  for(UInt i=1; i <= grid->getNumElements(); ++i)
  {
    elNodes = grid->getElementNodesL(i);
    GEOSHAPE::boundingBox(elNodes,Pmin,Pmax);
    
    lowerBound(i) = Pmin;
    upperBound(i) = Pmax;
    
    bbMax.setX( std::max(Pmax.getX(), bbMax.getX()) );
    bbMax.setY( std::max(Pmax.getY(), bbMax.getY()) );
    bbMax.setZ( std::max(Pmax.getZ(), bbMax.getZ()) );
    
    bbMin.setX( std::min(Pmin.getX(), bbMin.getX()) );
    bbMin.setY( std::min(Pmin.getY(), bbMin.getY()) );
    bbMin.setZ( std::min(Pmin.getZ(), bbMin.getZ()) );
  }
  
  //! OctTree feeding----------------------------------------------------------
  //Typedefs
  typedef std::set<UInt> SET;
  typedef std::set<UInt>::iterator ITER;
  
  //Clearing
  octTree.clear();
  
  //Alloc
  UInt maxSize;
  UInt handle;
  std::set<UInt> boxIndices;
  sVect<UInt> elements(grid->getNumElements());
  sVect<sVect<UInt> > blocks;
  sVect<UInt> temp, newWip, wip(1); wip(1) = 1;
  sVect<searchBoundingBox> leaves(8);
  
  //Oct tree startup 
  for(UInt i=1; i <= grid->getNumElements(); ++i)
  {elements(i) = i;} 
  
  octTree.addLeaves();
  octTree.getData().setElements(elements);
  octTree.getData().setBoundingBox(bbMin,bbMax);
  
  //Oct tree segmentation
  while(wip.size() != 0)
  {    
    for(UInt i=1; i <= wip.size(); ++i)
    {
      //Startup arrays
      blocks.clear();
      blocks.resize(8);
      
      //Tree positioning
      handle = wip(i);
      octTree.goTo(handle);
      
      //Elements segmentation
      elements = octTree.getData().getElements();
      
      for(UInt k=1; k <= elements.size(); ++k)
      {
	boxIndices = octTree.getData().octPart(lowerBound(elements(k)), upperBound(elements(k)));
	
	for(ITER it = boxIndices.begin(); it != boxIndices.end(); ++it)
	{
	  assert(*it <= 8);
	  blocks(*it).push_back(elements(k));
	}
      }
      
      //Maximum size
      maxSize = 0;
      
      for(UInt k=1; k <= 8; k++)
      { maxSize += UInt(blocks(k).size()); }     
      
      //Create new leaves
      if(maxSize < (elements.size() * OctToll))
      {
	Pmin = octTree.getData().getPmin(false,false,false);
	Pmax = octTree.getData().getPmax(false,false,false);
	leaves(1).setBoundingBox(Pmin,Pmax);
	leaves(1).setElements(blocks(1));
	
	Pmin = octTree.getData().getPmin(true,false,false);
	Pmax = octTree.getData().getPmax(true,false,false);
	leaves(2).setBoundingBox(Pmin,Pmax);
	leaves(2).setElements(blocks(2));
	
	Pmin = octTree.getData().getPmin(false,true,false);
	Pmax = octTree.getData().getPmax(false,true,false);
	leaves(3).setBoundingBox(Pmin,Pmax);
	leaves(3).setElements(blocks(3));
	
	Pmin = octTree.getData().getPmin(true,true,false);
	Pmax = octTree.getData().getPmax(true,true,false);
	leaves(4).setBoundingBox(Pmin,Pmax);
	leaves(4).setElements(blocks(4));
	
	Pmin = octTree.getData().getPmin(false,false,true);
	Pmax = octTree.getData().getPmax(false,false,true);
	leaves(5).setBoundingBox(Pmin,Pmax);
	leaves(5).setElements(blocks(5));
	
	Pmin = octTree.getData().getPmin(true,false,true);
	Pmax = octTree.getData().getPmax(true,false,true);
	leaves(6).setBoundingBox(Pmin,Pmax);
	leaves(6).setElements(blocks(6));
	
	Pmin = octTree.getData().getPmin(false,true,true);
	Pmax = octTree.getData().getPmax(false,true,true);
	leaves(7).setBoundingBox(Pmin,Pmax);
	leaves(7).setElements(blocks(7));
	
	Pmin = octTree.getData().getPmin(true,true,true);
	Pmax = octTree.getData().getPmax(true,true,true);
	leaves(8).setBoundingBox(Pmin,Pmax);
	leaves(8).setElements(blocks(8));
	
	temp = octTree.addLeaves(leaves(1),
				 leaves(2),
				 leaves(3),
				 leaves(4),
				 leaves(5),
				 leaves(6),
				 leaves(7),
				 leaves(8));
	
	newWip.push_back(temp(1));
	newWip.push_back(temp(2));
	newWip.push_back(temp(3));
	newWip.push_back(temp(4));
	newWip.push_back(temp(5));
	newWip.push_back(temp(6));
	newWip.push_back(temp(7));
	newWip.push_back(temp(8));
      }
    }
    
    //Update the leaves to be worked
    wip = newWip;
    newWip.clear();
  }
}

template<typename MESH, typename CONNECT>
const sVect<UInt>
geoOctTree<MESH,CONNECT>::
getMatchingElements(const point3d & P)
{
  bool ix, iy, iz;
  octTree.restart();
  
  while(!octTree.getMap().getIsLeaf())
  {
    octTree.getData().octPart(ix,iy,iz,P);
    octTree.goDown(ix,iy,iz);
  }
  
  return(octTree.getData().getElements());
}

template<typename MESH, typename CONNECT>
const sVect<UInt>
geoOctTree<MESH,CONNECT>::
getMatchingElementsExt(const point3d & P)
{
  //Assert---------------------------------------------------------------------
  assert(startupLocal);
  
  //Typedef
  typedef std::pair<Real,UInt>          DIST_ID;
  typedef multimap<Real,DIST_ID>::iterator ITER;
  typedef std::set<UInt>::iterator    UINT_ITER;
  
  //Alloc----------------------------------------------------------------------
  UInt count;
  UInt idFather, idSon, numElements;
  Real maxDist, minDist, target;
  multimap<Real,DIST_ID> blockList;
  DIST_ID distId;
  ITER itFather;
  sVect<ITER> deadList;
  
  //Startup--------------------------------------------------------------------
  octTree.restart();
  idFather = octTree.getMap().getId();
  minDist  = octTree.getData().getMinDist(P);
  maxDist  = octTree.getData().getMaxDist(P);
  target   = maxDist;
  
  octTree.getData(idSon).getMaxDist(P);
  
  distId.first  = minDist;
  distId.second = idFather;
  blockList.insert(std::pair<Real,DIST_ID>(maxDist,distId));
  
  //Selection loop-------------------------------------------------------------  
  while(true)
  {
    //Identify the first non-leaf
    count = 0;
    
    for(ITER it=blockList.begin(); it != blockList.end(); ++it)
    {
      //Extrac the block
      itFather = it;
      idFather = it->second.second;
      
      if(!octTree.getMap(idFather).getIsLeaf())
      { break; }
      
      count++;
    }
    
    //If all leaves exit
    if(blockList.size() == count)
    { break; }
    
    //Insert the non-empty childs
    blockList.erase(itFather);
    
    for(UInt k=1; k <= 8; ++k)
    {
      idSon       = octTree.getMap(idFather).getSons(k);
      minDist     = octTree.getData(idSon).getMinDist(P);
      maxDist     = octTree.getData(idSon).getMaxDist(P);
      numElements = octTree.getData(idSon).getElements().size();
      
      if(numElements != 0)
      {
	distId.first  = minDist;
	distId.second = idSon;
	blockList.insert(std::pair<Real,DIST_ID>(maxDist,distId));
      }
    }
    
    //Filtering
    target = std::min(target, blockList.begin()->first);
    deadList.clear();
    
    for(ITER it=blockList.begin(); it != blockList.end(); ++it)
    {
      if(it->second.first > target)
      { deadList.push_back(it); }
    }
    
    for(UInt k=1; k <= deadList.size(); ++k) 
    { blockList.erase(deadList(k)); }
  }
  
  //Output---------------------------------------------------------------------
  set<UInt>   outElements;
  sVect<UInt> tempElements;
  
  for(ITER it=blockList.begin(); it != blockList.end(); ++it)
  {    
    idFather     = it->second.second;
    tempElements = octTree.getData(idFather).getElements();
    
    for(UInt k=1; k <= tempElements.size(); ++k)
    { outElements.insert(tempElements(k)); }
  }
  
  tempElements.clear();
  
  for(UINT_ITER it = outElements.begin(); it != outElements.end(); ++it)
  { tempElements.push_back(*it); }
  
  return(tempElements);
}



#endif