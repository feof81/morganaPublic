/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MEMIMAGE_H
#define MEMIMAGE_H

#include "typesInterface.hpp"

#include "pMapItem.h"
#include "pMapItemSendRecv.h"
#include "pMapItemShare.h"
#include "simpleFormats.hpp"
#include "pMap.hpp"
#include "pMapManip.hpp"
#include "pVect.hpp"
#include "sOctTreeItem.h"
#include "sOctTree.hpp"
#include "dataBaseLin1d.hpp"

#include "pGraphItem.h"
#include "pGraphItemOriented.h"
#include "pGraphItemSubLoc.h"
#include "pGraph.hpp"
#include "pGraphManip.hpp"

#include "searchBoundingBox.h"
#include "searchCard.h"
#include "geoElement.hpp"
#include "pointElement.hpp"
#include "mesh1d.hpp"
#include "mesh2d.hpp"
#include "mesh3d.hpp"
#include "connect1d.hpp"
#include "connect2d.hpp"
#include "connect3d.hpp"
#include "geoOctTree.hpp"
#include "searchData.hpp"
#include "localSearch1dA.hpp"
#include "localSearch2dA.hpp"
#include "localSearch3dA.hpp"
#include "search1dA.hpp"
#include "search2dA.hpp"
#include "search3dA.hpp"


/*! Class to save to file and load from file run-time data */
class memImage
{
  /*! @name Links and flags */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructor and functions */ //@{
  public:
    memImage();
    memImage(const Teuchos::RCP<communicator> & CommDev);
    memImage(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Print to file - Dofs */ //@{
  public:
    void printMem(fstream & f, const bool & A) const;
    void printMem(fstream & f, const Int  & A) const;
    void printMem(fstream & f, const UInt & A) const;
    void printMem(fstream & f, const Real & A) const;
    void printMem(fstream & f, const komplex  & A) const;
    void printMem(fstream & f, const point2d  & A) const;
    void printMem(fstream & f, const point3d  & A) const;
    void printMem(fstream & f, const tensor2d & A) const;
    void printMem(fstream & f, const tensor3d & A) const;
    void printMem(fstream & f, const stateVector & A) const;
    void printMem(fstream & f, const stateMatrix & A) const;
    
    template<size_t N> 
    void printMem(fstream & f, const staticVector<N> & A) const;
    //@}
    
    /*! @name Load from file - Dofs */ //@{
  public:
    void loadMem(fstream & f, bool & A) const;
    void loadMem(fstream & f, Int  & A) const;
    void loadMem(fstream & f, UInt & A) const;
    void loadMem(fstream & f, Real & A) const;
    void loadMem(fstream & f, komplex  & A) const;
    void loadMem(fstream & f, point2d  & A) const;
    void loadMem(fstream & f, point3d  & A) const;
    void loadMem(fstream & f, tensor2d & A) const;
    void loadMem(fstream & f, tensor3d & A) const;
    void loadMem(fstream & f, stateVector & A) const;
    void loadMem(fstream & f, stateMatrix & A) const;
    
    template<size_t N> 
    void loadMem(fstream & f, staticVector<N> & A) const;
    //@}
    
    /*! @name Print to file - Container */ //@{
  public:
    void printMem(fstream & f, const pMapItem         & A) const;
    void printMem(fstream & f, const pMapItemShare    & A) const;
    void printMem(fstream & f, const pMapItemSendRecv & A) const;
    
    template<typename T>
    void printMem(fstream & f, sVect<T> & A) const;
    
    template<typename T>
    void printMem(fstream & f, sArray<T> & A) const;
    
    template<typename ITEM>
    void printMem(fstream & f, const pMap<ITEM> & A) const;
    
    template<typename ITEM>
    void printMem(fstream & f, const pMapManip<ITEM> & A) const;
    
    template<typename DATA, typename MAP>
    void printMem(fstream & f, const pVect<DATA,MAP> & A) const;
    
    void printMem(fstream & f, const sOctTreeItem & A) const;
    
    template<typename DATA>
    void printMem(fstream & f, const sOctTree<DATA> & A) const;
    
    template<typename X, typename Y>
    void printMem(fstream & f, const dataBaseLin1d<X,Y> & A) const;
    //@}
    
    /*! @name Load from file - Container */ //@{
  public:
    void loadMem(fstream & f, pMapItem & A) const;
    void loadMem(fstream & f, pMapItemShare & A) const;
    void loadMem(fstream & f, pMapItemSendRecv & A) const;
    
    template<typename T>
    void loadMem(fstream & f, sVect<T> & A) const;
    
    template<typename T>
    void loadMem(fstream & f, sArray<T> & A) const;
    
    template<typename ITEM>
    void loadMem(fstream & f, pMap<ITEM> & A) const;
    
    template<typename ITEM>
    void loadMem(fstream & f, pMapManip<ITEM> & A) const;
    
    template<typename DATA, typename MAP>
    void loadMem(fstream & f, pVect<DATA,MAP> & A) const;
    
    void loadMem(fstream & f, sOctTreeItem & A) const;
    
    template<typename DATA>
    void loadMem(fstream & f, sOctTree<DATA> & A) const;
    
    template<typename X, typename Y>
    void loadMem(fstream & f, dataBaseLin1d<X,Y> & A) const;
    //@}
    
    /*! @name Print to file - pGraph */ //@{
  public:
    void printMem(fstream & f, const pGraphItem         & A) const;
    void printMem(fstream & f, const pGraphItemOriented & A) const;
    void printMem(fstream & f, const pGraphItemSubLoc   & A) const;
    
    template<typename ITEM, typename ROWMAP, typename COLMAP>
    void printMem(fstream & f, const pGraph<ITEM,ROWMAP,COLMAP> & A) const;
    
    template<typename ITEM, typename ROWMAP, typename COLMAP>
    void printMem(fstream & f, const pGraphManip<ITEM,ROWMAP,COLMAP> & A) const;
    //@}
    
    /*! @name Load from file - pGraph */ //@{
  public:
    void loadMem(fstream & f, pGraphItem         & A) const;
    void loadMem(fstream & f, pGraphItemOriented & A) const;
    void loadMem(fstream & f, pGraphItemSubLoc   & A) const;
    
    template<typename ITEM, typename ROWMAP, typename COLMAP>
    void loadMem(fstream & f, pGraph<ITEM,ROWMAP,COLMAP> & A) const;
    
    template<typename ITEM, typename ROWMAP, typename COLMAP>
    void loadMem(fstream & f, pGraphManip<ITEM,ROWMAP,COLMAP> & A) const;
    //@}
    
    /*! @name Print to file - Geometry */ //@{
  public:
    template<typename GEOSHAPE>
    void printMem(fstream & f, const geoElement<GEOSHAPE> & A) const;
    
    template<typename GEOSHAPE>
    void printMem(fstream & f, const pointElement<GEOSHAPE> & A) const;
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void printMem(fstream & f, const mesh1d<GEOSHAPE,ELMAP,NODEMAP> & A) const;
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void printMem(fstream & f, const mesh2d<GEOSHAPE,ELMAP,NODEMAP> & A) const;
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void printMem(fstream & f, const mesh3d<GEOSHAPE,ELMAP,NODEMAP> & A) const;
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void printMem(fstream & f, const connect1d<GEOSHAPE,ELMAP,NODEMAP> & A) const;
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void printMem(fstream & f, const connect2d<GEOSHAPE,ELMAP,NODEMAP> & A) const;
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void printMem(fstream & f, const connect3d<GEOSHAPE,ELMAP,NODEMAP> & A) const;
    //@}
    
    /*! @name Load from file - Geometry */ //@{
  public:
    template<typename GEOSHAPE>
    void loadMem(fstream & f, geoElement<GEOSHAPE> & A) const;
    
    template<typename GEOSHAPE>
    void loadMem(fstream & f, pointElement<GEOSHAPE> & A) const;
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void loadMem(fstream & f, mesh1d<GEOSHAPE,ELMAP,NODEMAP> & A) const;
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void loadMem(fstream & f, mesh2d<GEOSHAPE,ELMAP,NODEMAP> & A) const;
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void loadMem(fstream & f, mesh3d<GEOSHAPE,ELMAP,NODEMAP> & A) const;
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void loadMem(fstream & f, connect1d<GEOSHAPE,ELMAP,NODEMAP> & A) const;
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void loadMem(fstream & f, connect2d<GEOSHAPE,ELMAP,NODEMAP> & A) const;
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void loadMem(fstream & f, connect3d<GEOSHAPE,ELMAP,NODEMAP> & A) const;
    //@}
};


//_________________________________________________________________________________________________
// DOFS
//-------------------------------------------------------------------------------------------------
template<size_t N> 
void
memImage::
printMem(fstream & f, const staticVector<N> & A) const
{
  for(UInt i=1; i <= UInt(N); ++i)
  { f.write((char*)&A(i), sizeof(Real)); }
}

template<size_t N> 
void
memImage::
loadMem(fstream & f, staticVector<N> & A) const
{
  for(UInt i=1; i <= UInt(N); ++i)
  { f.read((char*)&A(i), sizeof(Real)); }
}


//_________________________________________________________________________________________________
// PRINT TO FILE - CONTAINER
//-------------------------------------------------------------------------------------------------
template<typename T>
void
memImage::
printMem(fstream & f, sVect<T> & A) const
{
  UInt N = A.size();
  f.write((char*)&N, sizeof(UInt));
  
  for(UInt i=1; i <= N; ++i)
  { printMem(f,A(i)); }
}

template<typename T>
void
memImage::
printMem(fstream & f, sArray<T> & A) const
{
  UInt rows = A.nrows();
  UInt cols = A.ncols();
  
  f.write((char*)&rows, sizeof(UInt));
  f.write((char*)&cols, sizeof(UInt));
  
  for(UInt i=1; i <= rows; ++i)
  {
    for(UInt j=1; j <= cols; ++j)
    { printMem(f,A(i,j)); }
  }
}

template<typename ITEM>
void
memImage::
printMem(fstream & f, const pMap<ITEM> & A) const
{
  sVect<ITEM> outVect = A.items;
  printMem(f,outVect);
}

template<typename ITEM>
void
memImage::
printMem(fstream & f, const pMapManip<ITEM> & A) const
{
  typedef typename set<ITEM>::iterator ITER;
  
  printMem(f,A.mapLoaded);
  printMem(f,A.finderOk);
  
  UInt N = A.container.size();
  printMem(f,N);
  
  for(ITER iter = A.container.begin(); iter != A.container.end(); ++iter)
  { printMem(f,*iter); }
}

template<typename DATA, typename MAP>
void
memImage::
printMem(fstream & f, const pVect<DATA,MAP> & A) const
{
  printMem(f, A.startupOk);
  printMem(f,*A.data);
  printMem(f,*A.map);
  printMem(f,*A.mapManip);
}

template<typename DATA>
void
memImage::
printMem(fstream & f, const sOctTree<DATA> & A) const
{
  typedef std::map<sOctTreeItem,DATA>  DATAMAP;
  typedef typename DATAMAP::iterator   ITER;
  
  printMem(f,A.k);
  
  UInt N = A.data.size();
  printMem(f,N);
  
  DATAMAP data = A.data;

  for(ITER iter = data.begin(); iter != data.end(); ++iter)
  {
    printMem(f,iter->first);
    printMem(f,iter->second);
  }
}

template<typename X, typename Y>
void
memImage::
printMem(fstream & f, const dataBaseLin1d<X,Y> & A) const
{
  typedef std::map<X,Y>              DATAMAP;
  typedef typename DATAMAP::iterator ITER;
  
  printMem(f,A.scaleX);
  printMem(f,A.scaleY);
  
  UInt N = A.data.size();
  printMem(f,N);
  
  DATAMAP data = A.data;

  for(ITER iter = data.begin(); iter != data.end(); ++iter)
  {
    printMem(f,iter->first);
    printMem(f,iter->second);
  }
}


//_________________________________________________________________________________________________
// LOAD FROM FILE - CONTAINER
//-------------------------------------------------------------------------------------------------
template<typename T>
void
memImage::
loadMem(fstream & f, sVect<T> & A) const
{
  UInt N;
  f.read((char*)&N, sizeof(UInt));
  A.resize(N);
  
  T data;
  for(UInt i=1; i <= N; ++i)
  {
    loadMem(f,data);
    A(i) = data;
  }
}

template<typename T>
void
memImage::
loadMem(fstream & f, sArray<T> & A) const
{
  UInt rows; f.read((char*)&rows, sizeof(UInt));
  UInt cols; f.read((char*)&cols, sizeof(UInt));
  
  A.reshape(rows,cols);
  
  for(UInt i=1; i <= rows; ++i)
  {
    for(UInt j=1; j <= cols; ++j)
    { loadMem(f,A(i,j)); }
  }
}

template<typename ITEM>
void
memImage::
loadMem(fstream & f, pMap<ITEM> & A) const
{
  A.clear();
  loadMem(f,A.items);
}

template<typename ITEM>
void
memImage::
loadMem(fstream & f, pMapManip<ITEM> & A) const
{
  typedef typename set<ITEM>::iterator ITER;
  
  bool boolVal;
  UInt uintVal;
  ITEM itemVal;
  
  loadMem(f,boolVal);
  A.mapLoaded = boolVal;
  
  loadMem(f,boolVal);
  A.finderOk = boolVal;
  
  loadMem(f,uintVal);
  A.container.clear();
  
  for(UInt i=1; i <= uintVal; ++i)
  {
    loadMem(f,itemVal);
    A.container.insert(itemVal);
  }
}

template<typename DATA, typename MAP>
void
memImage::
loadMem(fstream & f, pVect<DATA,MAP> & A) const
{
  loadMem(f, A.startupOk);
  loadMem(f,*A.data);
  loadMem(f,*A.map);
  loadMem(f,*A.mapManip);
}

template<typename DATA>
void
memImage::
loadMem(fstream & f, sOctTree<DATA> & A) const
{
  typedef typename sOctTree<DATA>::DATAMAP  DATAMAP;
  typedef typename DATAMAP::iterator        ITER;
  
  UInt N;
  pair<sOctTreeItem,DATA> dataPair;
  
  loadMem(f,A.k);
  loadMem(f,N);
  
  A.data.clear();
  
  for(UInt i=1; i <= N; ++i)
  {
    loadMem(f,dataPair.first);
    loadMem(f,dataPair.second);
    
    A.data.insert(dataPair);
  }
  
  A.restart();
}

template<typename X, typename Y>
void
memImage::
loadMem(fstream & f, dataBaseLin1d<X,Y> & A) const
{
  typedef std::map<X,Y>              DATAMAP;
  typedef typename DATAMAP::iterator ITER;
  
  UInt N;
  pair<X,Y> dataPair;
  
  loadMem(f,A.scaleX);
  loadMem(f,A.scaleY);
  loadMem(f,N);
  
  A.data.clear();
  
  for(UInt i=1; i <= N; ++i)
  {
    loadMem(f,dataPair.first);
    loadMem(f,dataPair.second);
    
    A.data.insert(dataPair);
  }
}


//_________________________________________________________________________________________________
// PRINT TO FILE - PGRAPH
//-------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
void
memImage::
printMem(fstream & f, const pGraph<ITEM,ROWMAP,COLMAP> & A) const
{
  printMem(f, A.startupOk);
  printMem(f,*A.data);
  printMem(f,*A.map);
  printMem(f,*A.mapManip);
  
  printMem(f,A.isLocal);
  printMem(f,A.colStartupOk);
  printMem(f,A.colMap);
  printMem(f,A.colManip);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
memImage::
printMem(fstream & f, const pGraphManip<ITEM,ROWMAP,COLMAP> & A) const
{
  typedef std::map<ITEM,ROWMAP>      DATAMAP;
  typedef typename DATAMAP::iterator ITER;
  
  DATAMAP container = A.container;
  
  printMem(f,A.graphLoaded);
  printMem(f,A.finderOk);
  
  UInt N = container.size();
  printMem(f,N);
  
  for(ITER iter = container.begin(); iter != container.end(); ++iter)
  {
    printMem(f,iter->first);
    printMem(f,iter->second);
  }
}


//_________________________________________________________________________________________________
// LOAD FROM FILE - PGRAPH
//-------------------------------------------------------------------------------------------------
template<typename ITEM, typename ROWMAP, typename COLMAP>
void
memImage::
loadMem(fstream & f, pGraph<ITEM,ROWMAP,COLMAP> & A) const
{
  loadMem(f, A.startupOk);
  loadMem(f,*A.data);
  loadMem(f,*A.map);
  loadMem(f,*A.mapManip);
  
  loadMem(f,A.isLocal);
  loadMem(f,A.colStartupOk);
  loadMem(f,A.colMap);
  loadMem(f,A.colManip);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
memImage::
loadMem(fstream & f, pGraphManip<ITEM,ROWMAP,COLMAP> & A) const
{
  typedef std::map<ITEM,ROWMAP>      DATAMAP;
  typedef typename DATAMAP::iterator ITER;
  
  UInt N;
  pair<ITEM,ROWMAP> dataPair;
  
  loadMem(f,A.graphLoaded);
  loadMem(f,A.finderOk);
  loadMem(f,N);
  
  A.container.clear();
  
  for(UInt i=1; i <= N; ++i)
  {
    loadMem(f,dataPair.first);
    loadMem(f,dataPair.second);
    
    A.container.insert(dataPair);
  }
}


//_________________________________________________________________________________________________
// PRINT TO FILE - GEOMETRY
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE>
void
memImage::
printMem(fstream & f, const geoElement<GEOSHAPE> & A) const
{
  sVect<UInt> connected    = A.connected;
  sVect<UInt> orderedNodes = A.orderedNodes;
  
  printMem(f, connected);
  printMem(f, orderedNodes);
  printMem(f, A.nodesOrdered);
  printMem(f, A.fixedLength);
  printMem(f, A.geoId);
}
 
template<typename GEOSHAPE>
void
memImage::
printMem(fstream & f, const pointElement<GEOSHAPE> & A) const
{
  UInt geoId = A.geoId;
  sVect<point3d> points = A.points;
  
  printMem(f, geoId);
  printMem(f, points);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
memImage::
printMem(fstream & f, const mesh1d<GEOSHAPE,ELMAP,NODEMAP> & A) const
{
  printMem(f, A.verticesComputed);
  printMem(f, A.mapTransferred);
  
  UInt standard = UInt(A.standard);
  printMem(f, standard);
  printMem(f, A.numVertices);
  
  printMem(f, A.nodes);
  printMem(f, A.isVertex);
  printMem(f, A.elements);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
memImage::
printMem(fstream & f, const mesh2d<GEOSHAPE,ELMAP,NODEMAP> & A) const
{
  printMem(f, A.verticesComputed);
  printMem(f, A.mapTransferred);
  
  UInt standard = UInt(A.standard);
  printMem(f, standard);
  printMem(f, A.numVertices);
  
  printMem(f, A.nodes);
  printMem(f, A.isVertex);
  printMem(f, A.elements);
  printMem(f, A.edges);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
memImage::
printMem(fstream & f, const mesh3d<GEOSHAPE,ELMAP,NODEMAP> & A) const
{
  printMem(f, A.verticesComputed);
  printMem(f, A.mapTransferred);
  
  UInt standard = UInt(A.standard);
  printMem(f, standard);
  printMem(f, A.numVertices);
  
  printMem(f, A.nodes);
  printMem(f, A.isVertex);
  printMem(f, A.elements);
  printMem(f, A.faces);
  printMem(f, A.edges);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
memImage::
printMem(fstream & f, const connect1d<GEOSHAPE,ELMAP,NODEMAP> & A) const
{
  printMem(f, A.connectCreated);
  printMem(f, A.commDevLoaded);
  printMem(f, A.grid1dLoaded);
  
  printMem(f, A.vertexToVertex);
  printMem(f, A.vertexToElement);
  printMem(f, A.elementToElement);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
memImage::
printMem(fstream & f, const connect2d<GEOSHAPE,ELMAP,NODEMAP> & A) const
{
  printMem(f, A.connectCreated);
  printMem(f, A.boundaryConnectCreated);
  printMem(f, A.commDevLoaded);
  printMem(f, A.grid2dLoaded);
  printMem(f, A.grid1dLoaded);
  
  printMem(f, A.vertexToVertex);
  printMem(f, A.vertexToElement);
  printMem(f, A.vertexToEdge);
  printMem(f, A.edgeToElement);
  printMem(f, A.elementToEdge);
  printMem(f, A.elementToElement);
  
  printMem(f, A.vertexIsBoundary);
  printMem(f, A.elementIsBoundary);
  printMem(f, A.edgeIsBoundary);
   
  printMem(f, A.vertexBVertex);
  printMem(f, A.edgeBEdge);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
memImage::
printMem(fstream & f, const connect3d<GEOSHAPE,ELMAP,NODEMAP> & A) const
{
  printMem(f, A.connectCreated);
  printMem(f, A.boundaryConnectCreated);
  printMem(f, A.commDevLoaded);
  printMem(f, A.grid3dLoaded);
  printMem(f, A.grid2dLoaded);
  
  printMem(f, A.vertexToVertex);
  printMem(f, A.vertexToElement);
  printMem(f, A.vertexToEdge);
  printMem(f, A.faceToElement);
  printMem(f, A.elementToEdge);
  printMem(f, A.elementToFace);
  printMem(f, A.elementToElement);
  
  printMem(f, A.vertexIsBoundary);
  printMem(f, A.elementIsBoundary);
  printMem(f, A.faceIsBoundary);
  printMem(f, A.edgeIsBoundary);
    
  printMem(f, A.vertexBVertex);
  printMem(f, A.faceBFace);
  printMem(f, A.edgeBEdge);
}


//_________________________________________________________________________________________________
// LOAD FROM FILE - GEOMETRY
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE>
void
memImage::
loadMem(fstream & f, geoElement<GEOSHAPE> & A) const
{
  sVect<UInt> connected;
  sVect<UInt> orderedNodes;
  bool nodesOrdered;
  bool fixedLength;
  UInt geoId;
  
  loadMem(f, connected);
  loadMem(f, orderedNodes);
  loadMem(f, nodesOrdered);
  loadMem(f, fixedLength);
  loadMem(f, geoId);
  
  A.connected    = connected;
  A.orderedNodes = orderedNodes;
  A.nodesOrdered = nodesOrdered;
  A.fixedLength  = fixedLength;
  A.geoId        = geoId;
}

template<typename GEOSHAPE>
void
memImage::
loadMem(fstream & f, pointElement<GEOSHAPE> & A) const
{
  loadMem(f, A.geoId);
  loadMem(f, A.points);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
memImage::
loadMem(fstream & f, mesh1d<GEOSHAPE,ELMAP,NODEMAP> & A) const
{
  loadMem(f, A.verticesComputed);
  loadMem(f, A.mapTransferred);
  
  UInt standard;
  loadMem(f, standard);
  A.standard = MeshStandards(standard);
  
  loadMem(f, A.numVertices);
  
  loadMem(f, A.nodes);
  loadMem(f, A.isVertex);
  loadMem(f, A.elements);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
memImage::
loadMem(fstream & f, mesh2d<GEOSHAPE,ELMAP,NODEMAP> & A) const
{
  loadMem(f, A.verticesComputed);
  loadMem(f, A.mapTransferred);
  
  UInt standard;
  loadMem(f, standard);
  A.standard = MeshStandards(standard);
  
  loadMem(f, A.numVertices);
  
  loadMem(f, A.nodes);
  loadMem(f, A.isVertex);
  loadMem(f, A.elements);
  loadMem(f, A.edges);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
memImage::
loadMem(fstream & f, mesh3d<GEOSHAPE,ELMAP,NODEMAP> & A) const
{
  loadMem(f, A.verticesComputed);
  loadMem(f, A.mapTransferred);
  
  UInt standard;
  loadMem(f, standard);
  A.standard = MeshStandards(standard);
  
  loadMem(f, A.numVertices);
  
  loadMem(f, A.nodes);
  loadMem(f, A.isVertex);
  loadMem(f, A.elements);
  loadMem(f, A.faces);
  loadMem(f, A.edges);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
memImage::
loadMem(fstream & f, connect1d<GEOSHAPE,ELMAP,NODEMAP> & A) const
{
  loadMem(f, A.connectCreated);
  loadMem(f, A.commDevLoaded);
  loadMem(f, A.grid1dLoaded);
  
  loadMem(f, A.vertexToVertex);
  loadMem(f, A.vertexToElement);
  loadMem(f, A.elementToElement);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
memImage::
loadMem(fstream & f, connect2d<GEOSHAPE,ELMAP,NODEMAP> & A) const
{
  loadMem(f, A.connectCreated);
  loadMem(f, A.boundaryConnectCreated);
  loadMem(f, A.commDevLoaded);
  loadMem(f, A.grid2dLoaded);
  loadMem(f, A.grid1dLoaded);
  
  loadMem(f, A.vertexToVertex);
  loadMem(f, A.vertexToElement);
  loadMem(f, A.vertexToEdge);
  loadMem(f, A.edgeToElement);
  loadMem(f, A.elementToEdge);
  loadMem(f, A.elementToElement);
  
  loadMem(f, A.vertexIsBoundary);
  loadMem(f, A.elementIsBoundary);
  loadMem(f, A.edgeIsBoundary);
   
  loadMem(f, A.vertexBVertex);
  loadMem(f, A.edgeBEdge);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
memImage::
loadMem(fstream & f, connect3d<GEOSHAPE,ELMAP,NODEMAP> & A) const
{
  loadMem(f, A.connectCreated);
  loadMem(f, A.boundaryConnectCreated);
  loadMem(f, A.commDevLoaded);
  loadMem(f, A.grid3dLoaded);
  loadMem(f, A.grid2dLoaded);
  
  loadMem(f, A.vertexToVertex);
  loadMem(f, A.vertexToElement);
  loadMem(f, A.vertexToEdge);
  loadMem(f, A.faceToElement);
  loadMem(f, A.elementToEdge);
  loadMem(f, A.elementToFace);
  loadMem(f, A.elementToElement);
  
  loadMem(f, A.vertexIsBoundary);
  loadMem(f, A.elementIsBoundary);
  loadMem(f, A.faceIsBoundary);
  loadMem(f, A.edgeIsBoundary);
    
  loadMem(f, A.vertexBVertex);
  loadMem(f, A.faceBFace);
  loadMem(f, A.edgeBEdge);
}

#endif
