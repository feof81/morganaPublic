/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef LGMSINGLELIST_HPP
#define LGMSINGLELIST_HPP

#include "pMap.hpp"
#include "pointElement.hpp"

#include "mesh3d.hpp"
#include "mesh2d.hpp"
#include "mesh1d.hpp"

#include "connect3d.hpp"
#include "connect2d.hpp"
#include "connect1d.hpp"


//_________________________________________________________________________________________________
// DUMMY EMPTY CLASS
//-------------------------------------------------------------------------------------------------
template<typename MESH, typename CONNECT, UInt SRCDIM, UInt TGTDIM>
class lgmSingleList
{
  /*! @name Typedefs */ //@{
  public:
    typedef typename MESH::GRID_GEOSHAPE    GEOSHAPE;
    typedef typename MESH::GRID_ELMAP       PMAPTYPE;
    typedef pMap<PMAPTYPE>                  PMAP;
    typedef pointElement<GEOSHAPE>          POINTELEMENT;
    typedef std::map<POINTELEMENT,PMAPTYPE> STDMAP;
    typedef GEOSHAPE                        OUTGEOSHAPE;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool gridOk;
    Teuchos::RCP<MESH>    grid;
    Teuchos::RCP<CONNECT> connectGrid;
    //@}
    
    /*! @name Functions */ //@{
  public:
    lgmSingleList();
    lgmSingleList(const Teuchos::RCP<MESH>    & Grid,
                  const Teuchos::RCP<CONNECT> & ConnectGrid);
    
    void setMesh(const Teuchos::RCP<MESH>    & Grid,
                 const Teuchos::RCP<CONNECT> & ConnectGrid);
    
    STDMAP getList();
    sVect<point3d> getNodes(const PMAPTYPE & pmItem);
    //@}
};

template<typename MESH, typename CONNECT, UInt SRCDIM, UInt TGTDIM>
lgmSingleList<MESH,CONNECT,SRCDIM,TGTDIM>::
lgmSingleList()
{
}

template<typename MESH, typename CONNECT, UInt SRCDIM, UInt TGTDIM>
lgmSingleList<MESH,CONNECT,SRCDIM,TGTDIM>::
lgmSingleList(const Teuchos::RCP<MESH>    & Grid,
              const Teuchos::RCP<CONNECT> & ConnectGrid)
{
}

template<typename MESH, typename CONNECT, UInt SRCDIM, UInt TGTDIM>
void
lgmSingleList<MESH,CONNECT,SRCDIM,TGTDIM>::
setMesh(const Teuchos::RCP<MESH>    & Grid,
        const Teuchos::RCP<CONNECT> & ConnectGrid)
{
}

template<typename MESH, typename CONNECT, UInt SRCDIM, UInt TGTDIM>
typename lgmSingleList<MESH,CONNECT,SRCDIM,TGTDIM>::STDMAP
lgmSingleList<MESH,CONNECT,SRCDIM,TGTDIM>::
getList()
{
  STDMAP outMap;
  return(outMap);
}

template<typename MESH, typename CONNECT, UInt SRCDIM, UInt TGTDIM>
sVect<point3d>
lgmSingleList<MESH,CONNECT,SRCDIM,TGTDIM>::
getNodes(const PMAPTYPE & pmItem)
{
  sVect<point3d> voidVect;
  return(voidVect);
}


//_________________________________________________________________________________________________
// 3d Mesh - 3d Elements
//-------------------------------------------------------------------------------------------------
template<typename MESH, typename CONNECT>
class lgmSingleList<MESH,CONNECT,3,3>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename MESH::GEOSHAPE3D       GEOSHAPE3D;
    typedef typename MESH::GRID_ELMAP       PMAPTYPE;
    typedef pMap<PMAPTYPE>                  PMAP;
    typedef pointElement<GEOSHAPE3D>        POINTELEMENT;
    typedef std::map<POINTELEMENT,PMAPTYPE> STDMAP;
    typedef GEOSHAPE3D                      OUTGEOSHAPE;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool gridOk;
    Teuchos::RCP<MESH>    grid;
    Teuchos::RCP<CONNECT> connectGrid;
    //@}
    
    /*! @name Functions */ //@{
  public:
    lgmSingleList();
    lgmSingleList(const Teuchos::RCP<MESH>    & Grid,
                  const Teuchos::RCP<CONNECT> & ConnectGrid);
    
    void setMesh(const Teuchos::RCP<MESH>    & Grid,
                 const Teuchos::RCP<CONNECT> & ConnectGrid);
    
    STDMAP getList();
    sVect<point3d> getNodes(const PMAPTYPE & pmItem);
    //@}
};

template<typename MESH, typename CONNECT>
lgmSingleList<MESH,CONNECT,3,3>::
lgmSingleList()
{
  typedef typename MESH::GRID_GEOSHAPE  MESH_GEOSHAPE;
  typedef typename MESH::GRID_ELMAP     MESH_ELMAP;
  typedef typename MESH::GRID_NODEMAP   MESH_NODEMAP;
  
  typedef typename CONNECT::GRID_GEOSHAPE CONNECT_GEOSHAPE;
  typedef typename CONNECT::GRID_ELMAP    CONNECT_ELMAP;
  typedef typename CONNECT::GRID_NODEMAP  CONNECT_NODEMAP;
  
  assert(staticAssert<MESH_GEOSHAPE::geoName     == CONNECT_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<MESH_ELMAP::parallelType   == CONNECT_ELMAP::parallelType>::returnValue);
  assert(staticAssert<MESH_NODEMAP::parallelType == CONNECT_NODEMAP::parallelType>::returnValue);
  assert(staticAssert<GEOSHAPE3D::nDim == 3>::returnValue);
  
  gridOk = false;
}

template<typename MESH, typename CONNECT>
lgmSingleList<MESH,CONNECT,3,3>::
lgmSingleList(const Teuchos::RCP<MESH>    & Grid,
              const Teuchos::RCP<CONNECT> & ConnectGrid)
{
  typedef typename MESH::GRID_GEOSHAPE  MESH_GEOSHAPE;
  typedef typename MESH::GRID_ELMAP     MESH_ELMAP;
  typedef typename MESH::GRID_NODEMAP   MESH_NODEMAP;
  
  typedef typename CONNECT::GRID_GEOSHAPE CONNECT_GEOSHAPE;
  typedef typename CONNECT::GRID_ELMAP    CONNECT_ELMAP;
  typedef typename CONNECT::GRID_NODEMAP  CONNECT_NODEMAP;
  
  assert(staticAssert<MESH_GEOSHAPE::geoName     == CONNECT_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<MESH_ELMAP::parallelType   == CONNECT_ELMAP::parallelType>::returnValue);
  assert(staticAssert<MESH_NODEMAP::parallelType == CONNECT_NODEMAP::parallelType>::returnValue);
  assert(staticAssert<GEOSHAPE3D::nDim == 3>::returnValue);
  
  gridOk = true;
  
  grid        = Grid;
  connectGrid = ConnectGrid;
}

template<typename MESH, typename CONNECT>
void
lgmSingleList<MESH,CONNECT,3,3>::
setMesh(const Teuchos::RCP<MESH>    & Grid,
        const Teuchos::RCP<CONNECT> & ConnectGrid)
{
  gridOk = true;
  
  grid        = Grid;
  connectGrid = ConnectGrid;
}

template<typename MESH, typename CONNECT>
typename lgmSingleList<MESH,CONNECT,3,3>::STDMAP
lgmSingleList<MESH,CONNECT,3,3>::
getList()
{
  //Assert
  assert(gridOk);
  
  //Typedef
  typedef std::pair<POINTELEMENT,PMAPTYPE> PAIO;
  
  //Alloc
  PAIO paio;
  pointElement<GEOSHAPE3D> pElement;
  STDMAP outMap;
  
  //Loop
  for(UInt el=1; el <= grid->getNumElements(); ++el)
  {
    pElement.setPoints(grid->getElementNodesL(el));
    pElement.reorder();
    
    paio.first  = pElement;
    paio.second = grid->getElements().getRowMapL(el);
    
    outMap.insert(paio);
  }
  
  return(outMap); 
}

template<typename MESH, typename CONNECT>
sVect<point3d>
lgmSingleList<MESH,CONNECT,3,3>::
getNodes(const PMAPTYPE & pmItem)
{
  assert(gridOk);
  return(grid->getElementNodesL(pmItem.getLid()));
}


//_________________________________________________________________________________________________
// 3d Mesh - 2d Faces
//-------------------------------------------------------------------------------------------------
template<typename MESH, typename CONNECT>
class lgmSingleList<MESH,CONNECT,3,2>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename MESH::GEOSHAPE3D  GEOSHAPE3D;
    typedef typename MESH::GEOSHAPE2D  GEOSHAPE2D;
    typedef typename MESH::GRID_ELMAP  PMAPTYPE;
    typedef pMap<PMAPTYPE>             PMAP;
    typedef pointElement<GEOSHAPE2D>   POINTELEMENT;
    typedef std::map<POINTELEMENT,PMAPTYPE> STDMAP;
    typedef GEOSHAPE2D                 OUTGEOSHAPE;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool gridOk;
    Teuchos::RCP<MESH>    grid;
    Teuchos::RCP<CONNECT> connectGrid;
    //@}
    
    /*! @name Functions */ //@{
  public:
    lgmSingleList();
    lgmSingleList(const Teuchos::RCP<MESH>    & Grid,
                  const Teuchos::RCP<CONNECT> & ConnectGrid);
    
    void setMesh(const Teuchos::RCP<MESH>    & Grid,
                 const Teuchos::RCP<CONNECT> & ConnectGrid);
    
    STDMAP getList();
    sVect<point3d> getNodes(const PMAPTYPE & pmItem);
    //@}
};

template<typename MESH, typename CONNECT>
lgmSingleList<MESH,CONNECT,3,2>::
lgmSingleList()
{
  typedef typename MESH::GRID_GEOSHAPE  MESH_GEOSHAPE;
  typedef typename MESH::GRID_ELMAP     MESH_ELMAP;
  typedef typename MESH::GRID_NODEMAP   MESH_NODEMAP;
  
  typedef typename CONNECT::GRID_GEOSHAPE CONNECT_GEOSHAPE;
  typedef typename CONNECT::GRID_ELMAP    CONNECT_ELMAP;
  typedef typename CONNECT::GRID_NODEMAP  CONNECT_NODEMAP;
  
  assert(staticAssert<MESH_GEOSHAPE::geoName     == CONNECT_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<MESH_ELMAP::parallelType   == CONNECT_ELMAP::parallelType>::returnValue);
  assert(staticAssert<MESH_NODEMAP::parallelType == CONNECT_NODEMAP::parallelType>::returnValue);
  assert(staticAssert<GEOSHAPE3D::nDim == 3>::returnValue);
  assert(staticAssert<GEOSHAPE2D::nDim == 2>::returnValue);
  
  gridOk = false;
}

template<typename MESH, typename CONNECT>
lgmSingleList<MESH,CONNECT,3,2>::
lgmSingleList(const Teuchos::RCP<MESH>    & Grid,
              const Teuchos::RCP<CONNECT> & ConnectGrid)
{
  typedef typename MESH::GRID_GEOSHAPE  MESH_GEOSHAPE;
  typedef typename MESH::GRID_ELMAP     MESH_ELMAP;
  typedef typename MESH::GRID_NODEMAP   MESH_NODEMAP;
  
  typedef typename CONNECT::GRID_GEOSHAPE CONNECT_GEOSHAPE;
  typedef typename CONNECT::GRID_ELMAP    CONNECT_ELMAP;
  typedef typename CONNECT::GRID_NODEMAP  CONNECT_NODEMAP;
  
  assert(staticAssert<MESH_GEOSHAPE::geoName     == CONNECT_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<MESH_ELMAP::parallelType   == CONNECT_ELMAP::parallelType>::returnValue);
  assert(staticAssert<MESH_NODEMAP::parallelType == CONNECT_NODEMAP::parallelType>::returnValue);
  assert(staticAssert<GEOSHAPE3D::nDim == 3>::returnValue);
  assert(staticAssert<GEOSHAPE2D::nDim == 2>::returnValue);
  
  gridOk = true;
  
  grid        = Grid;
  connectGrid = ConnectGrid;
}

template<typename MESH, typename CONNECT>
void
lgmSingleList<MESH,CONNECT,3,2>::
setMesh(const Teuchos::RCP<MESH>    & Grid,
        const Teuchos::RCP<CONNECT> & ConnectGrid)
{
  gridOk = true;
  
  grid        = Grid;
  connectGrid = ConnectGrid;
}

template<typename MESH, typename CONNECT>
typename lgmSingleList<MESH,CONNECT,3,2>::STDMAP
lgmSingleList<MESH,CONNECT,3,2>::
getList()
{
  //Assert
  assert(gridOk);
  
  //Typedef
  typedef std::pair<POINTELEMENT,PMAPTYPE> PAIO;
  
  //Alloc
  PAIO paio;
  pointElement<GEOSHAPE2D> pElement;
  STDMAP outMap;
  UInt el;
  
  //Loop
  for(UInt fc=1; fc <= grid->getNumFaces(); ++fc)
  {
    pElement.setPoints(grid->getFaceNodesL(fc));
    pElement.reorder();
    
    el = connectGrid->getFaceToElement(fc,1);
    
    paio.first  = pElement;
    paio.second = grid->getElements().getRowMapL(el); 
    
    outMap.insert(paio);
  }
  
  return(outMap); 
}

template<typename MESH, typename CONNECT>
sVect<point3d>
lgmSingleList<MESH,CONNECT,3,2>::
getNodes(const PMAPTYPE & pmItem)
{
  assert(gridOk);
  return(grid->getElementNodesL(pmItem.getLid()));
}


//_________________________________________________________________________________________________
// 3d Mesh - 1d Edges
//-------------------------------------------------------------------------------------------------
template<typename MESH, typename CONNECT>
class lgmSingleList<MESH,CONNECT,3,1>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename MESH::GEOSHAPE3D  GEOSHAPE3D;
    typedef typename MESH::GEOSHAPE2D  GEOSHAPE2D;
    typedef typename MESH::GEOSHAPE1D  GEOSHAPE1D;
    typedef typename MESH::GRID_ELMAP  PMAPTYPE;
    typedef pMap<PMAPTYPE>             PMAP;
    typedef pointElement<GEOSHAPE1D>   POINTELEMENT;
    typedef std::map<POINTELEMENT,PMAPTYPE> STDMAP;
    typedef GEOSHAPE1D                 OUTGEOSHAPE;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool gridOk;
    Teuchos::RCP<MESH>    grid;
    Teuchos::RCP<CONNECT> connectGrid;
    //@}
    
    /*! @name Functions */ //@{
  public:
    lgmSingleList();
    lgmSingleList(const Teuchos::RCP<MESH>    & Grid,
                  const Teuchos::RCP<CONNECT> & ConnectGrid);
    
    void setMesh(const Teuchos::RCP<MESH>    & Grid,
                 const Teuchos::RCP<CONNECT> & ConnectGrid);
    
    STDMAP getList();
    sVect<point3d> getNodes(const PMAPTYPE & pmItem);
    //@}
};

template<typename MESH, typename CONNECT>
lgmSingleList<MESH,CONNECT,3,1>::
lgmSingleList()
{
  typedef typename MESH::GRID_GEOSHAPE  MESH_GEOSHAPE;
  typedef typename MESH::GRID_ELMAP     MESH_ELMAP;
  typedef typename MESH::GRID_NODEMAP   MESH_NODEMAP;
  
  typedef typename CONNECT::GRID_GEOSHAPE CONNECT_GEOSHAPE;
  typedef typename CONNECT::GRID_ELMAP    CONNECT_ELMAP;
  typedef typename CONNECT::GRID_NODEMAP  CONNECT_NODEMAP;
  
  assert(staticAssert<MESH_GEOSHAPE::geoName     == CONNECT_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<MESH_ELMAP::parallelType   == CONNECT_ELMAP::parallelType>::returnValue);
  assert(staticAssert<MESH_NODEMAP::parallelType == CONNECT_NODEMAP::parallelType>::returnValue);
  assert(staticAssert<GEOSHAPE3D::nDim == 3>::returnValue);
  assert(staticAssert<GEOSHAPE2D::nDim == 2>::returnValue);
  assert(staticAssert<GEOSHAPE1D::nDim == 1>::returnValue);
  
  gridOk = false;
}

template<typename MESH, typename CONNECT>
lgmSingleList<MESH,CONNECT,3,1>::
lgmSingleList(const Teuchos::RCP<MESH>    & Grid,
              const Teuchos::RCP<CONNECT> & ConnectGrid)
{
  typedef typename MESH::GRID_GEOSHAPE  MESH_GEOSHAPE;
  typedef typename MESH::GRID_ELMAP     MESH_ELMAP;
  typedef typename MESH::GRID_NODEMAP   MESH_NODEMAP;
  
  typedef typename CONNECT::GRID_GEOSHAPE CONNECT_GEOSHAPE;
  typedef typename CONNECT::GRID_ELMAP    CONNECT_ELMAP;
  typedef typename CONNECT::GRID_NODEMAP  CONNECT_NODEMAP;
  
  assert(staticAssert<MESH_GEOSHAPE::geoName     == CONNECT_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<MESH_ELMAP::parallelType   == CONNECT_ELMAP::parallelType>::returnValue);
  assert(staticAssert<MESH_NODEMAP::parallelType == CONNECT_NODEMAP::parallelType>::returnValue);
  assert(staticAssert<GEOSHAPE3D::nDim == 3>::returnValue);
  assert(staticAssert<GEOSHAPE2D::nDim == 2>::returnValue);
  assert(staticAssert<GEOSHAPE1D::nDim == 1>::returnValue);
  
  gridOk = true;
  
  grid        = Grid;
  connectGrid = ConnectGrid;
}

template<typename MESH, typename CONNECT>
void
lgmSingleList<MESH,CONNECT,3,1>::
setMesh(const Teuchos::RCP<MESH>    & Grid,
        const Teuchos::RCP<CONNECT> & ConnectGrid)
{
  gridOk = true;
  
  grid        = Grid;
  connectGrid = ConnectGrid;
}

template<typename MESH, typename CONNECT>
typename lgmSingleList<MESH,CONNECT,3,1>::STDMAP
lgmSingleList<MESH,CONNECT,3,1>::
getList()
{
  //Assert
  assert(gridOk);
  
  //Typedef
  typedef std::pair<POINTELEMENT,PMAPTYPE> PAIO;
  
  //Alloc
  PAIO   paio;
  STDMAP outMap;
  sVect<point3d> elementNodes;
  pointElement<GEOSHAPE1D>    pElement;
  geoMapInterface<GEOSHAPE1D> refShape1d;
  geoMapInterface<GEOSHAPE3D> refShape3d;
  
  //Loop
  for(UInt el=1; el <= grid->getNumElements(); ++el)
  {
    elementNodes = grid->getElementNodesL(el);
    
    for(UInt j=1; j <= refShape3d.getNumEdges(); ++j)
    {
      for(UInt k=1; k <= refShape1d.getNumPoints(); ++k)
      { pElement.setPoint(k,elementNodes(refShape3d.edgeToPoint(j,k))); }
      
      pElement.reorder();
      
      paio.first  = pElement;
      paio.second = grid->getElements().getRowMapL(el);
      
      outMap.insert(paio);
    }
  }
  
  return(outMap); 
}

template<typename MESH, typename CONNECT>
sVect<point3d>
lgmSingleList<MESH,CONNECT,3,1>::
getNodes(const PMAPTYPE & pmItem)
{
  assert(gridOk);
  return(grid->getElementNodesL(pmItem.getLid()));
}


//_________________________________________________________________________________________________
// 2d Mesh - 2d Elements
//-------------------------------------------------------------------------------------------------
template<typename MESH, typename CONNECT>
class lgmSingleList<MESH,CONNECT,2,2>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename MESH::GEOSHAPE2D  GEOSHAPE2D;
    typedef typename MESH::GRID_ELMAP  PMAPTYPE;
    typedef pMap<PMAPTYPE>             PMAP;
    typedef pointElement<GEOSHAPE2D>   POINTELEMENT;
    typedef std::map<POINTELEMENT,PMAPTYPE> STDMAP;
    typedef GEOSHAPE2D                 OUTGEOSHAPE;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool gridOk;
    Teuchos::RCP<MESH>    grid;
    Teuchos::RCP<CONNECT> connectGrid;
    //@}
    
    /*! @name Functions */ //@{
  public:
    lgmSingleList();
    lgmSingleList(const Teuchos::RCP<MESH>    & Grid,
                  const Teuchos::RCP<CONNECT> & ConnectGrid);
    
    void setMesh(const Teuchos::RCP<MESH>    & Grid,
                 const Teuchos::RCP<CONNECT> & ConnectGrid);
    
    STDMAP getList();
    sVect<point3d> getNodes(const PMAPTYPE & pmItem);
    //@}
};

template<typename MESH, typename CONNECT>
lgmSingleList<MESH,CONNECT,2,2>::
lgmSingleList()
{
  typedef typename MESH::GRID_GEOSHAPE  MESH_GEOSHAPE;
  typedef typename MESH::GRID_ELMAP     MESH_ELMAP;
  typedef typename MESH::GRID_NODEMAP   MESH_NODEMAP;
  
  typedef typename CONNECT::GRID_GEOSHAPE CONNECT_GEOSHAPE;
  typedef typename CONNECT::GRID_ELMAP    CONNECT_ELMAP;
  typedef typename CONNECT::GRID_NODEMAP  CONNECT_NODEMAP;
  
  assert(staticAssert<MESH_GEOSHAPE::geoName     == CONNECT_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<MESH_ELMAP::parallelType   == CONNECT_ELMAP::parallelType>::returnValue);
  assert(staticAssert<MESH_NODEMAP::parallelType == CONNECT_NODEMAP::parallelType>::returnValue);
  assert(staticAssert<GEOSHAPE2D::nDim == 2>::returnValue);
  
  gridOk = false;
}

template<typename MESH, typename CONNECT>
lgmSingleList<MESH,CONNECT,2,2>::
lgmSingleList(const Teuchos::RCP<MESH>    & Grid,
              const Teuchos::RCP<CONNECT> & ConnectGrid)
{
  typedef typename MESH::GRID_GEOSHAPE  MESH_GEOSHAPE;
  typedef typename MESH::GRID_ELMAP     MESH_ELMAP;
  typedef typename MESH::GRID_NODEMAP   MESH_NODEMAP;
  
  typedef typename CONNECT::GRID_GEOSHAPE CONNECT_GEOSHAPE;
  typedef typename CONNECT::GRID_ELMAP    CONNECT_ELMAP;
  typedef typename CONNECT::GRID_NODEMAP  CONNECT_NODEMAP;
  
  assert(staticAssert<MESH_GEOSHAPE::geoName     == CONNECT_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<MESH_ELMAP::parallelType   == CONNECT_ELMAP::parallelType>::returnValue);
  assert(staticAssert<MESH_NODEMAP::parallelType == CONNECT_NODEMAP::parallelType>::returnValue);
  assert(staticAssert<GEOSHAPE2D::nDim == 2>::returnValue);
  
  gridOk = true;
  
  grid        = Grid;
  connectGrid = ConnectGrid;
}

template<typename MESH, typename CONNECT>
void
lgmSingleList<MESH,CONNECT,2,2>::
setMesh(const Teuchos::RCP<MESH>    & Grid,
        const Teuchos::RCP<CONNECT> & ConnectGrid)
{
  gridOk = true;
  
  grid        = Grid;
  connectGrid = ConnectGrid;
}

template<typename MESH, typename CONNECT>
typename lgmSingleList<MESH,CONNECT,2,2>::STDMAP
lgmSingleList<MESH,CONNECT,2,2>::
getList()
{
  //Assert
  assert(gridOk);
  
  //Typedef
  typedef std::pair<POINTELEMENT,PMAPTYPE> PAIO;
  
  //Alloc
  PAIO paio;
  pointElement<GEOSHAPE2D> pElement;
  STDMAP outMap;
  
  //Loop
  for(UInt el=1; el <= grid->getNumElements(); ++el)
  {
    pElement.setPoints(grid->getElementNodesL(el));
    pElement.reorder();
    
    paio.first  = pElement;
    paio.second = grid->getElements().getRowMapL(el);
    
    outMap.insert(paio);
  }
  
  return(outMap); 
}

template<typename MESH, typename CONNECT>
sVect<point3d>
lgmSingleList<MESH,CONNECT,2,2>::
getNodes(const PMAPTYPE & pmItem)
{
  assert(gridOk);
  return(grid->getElementNodesL(pmItem.getLid()));
}


//_________________________________________________________________________________________________
// 2d Mesh - 1d Edges
//-------------------------------------------------------------------------------------------------
template<typename MESH, typename CONNECT>
class lgmSingleList<MESH,CONNECT,2,1>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename MESH::GEOSHAPE2D  GEOSHAPE2D;
    typedef typename MESH::GEOSHAPE1D  GEOSHAPE1D;
    typedef typename MESH::GRID_ELMAP  PMAPTYPE;
    typedef pMap<PMAPTYPE>             PMAP;
    typedef pointElement<GEOSHAPE1D>   POINTELEMENT;
    typedef std::map<POINTELEMENT,PMAPTYPE> STDMAP;
    typedef GEOSHAPE1D                 OUTGEOSHAPE;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool gridOk;
    Teuchos::RCP<MESH>    grid;
    Teuchos::RCP<CONNECT> connectGrid;
    //@}
    
    /*! @name Functions */ //@{
  public:
    lgmSingleList();
    lgmSingleList(const Teuchos::RCP<MESH>    & Grid,
                  const Teuchos::RCP<CONNECT> & ConnectGrid);
    
    void setMesh(const Teuchos::RCP<MESH>    & Grid,
                 const Teuchos::RCP<CONNECT> & ConnectGrid);
    
    STDMAP getList();
    sVect<point3d> getNodes(const PMAPTYPE & pmItem);
    //@}
};

template<typename MESH, typename CONNECT>
lgmSingleList<MESH,CONNECT,2,1>::
lgmSingleList()
{
  typedef typename MESH::GRID_GEOSHAPE  MESH_GEOSHAPE;
  typedef typename MESH::GRID_ELMAP     MESH_ELMAP;
  typedef typename MESH::GRID_NODEMAP   MESH_NODEMAP;
  
  typedef typename CONNECT::GRID_GEOSHAPE CONNECT_GEOSHAPE;
  typedef typename CONNECT::GRID_ELMAP    CONNECT_ELMAP;
  typedef typename CONNECT::GRID_NODEMAP  CONNECT_NODEMAP;
  
  assert(staticAssert<MESH_GEOSHAPE::geoName     == CONNECT_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<MESH_ELMAP::parallelType   == CONNECT_ELMAP::parallelType>::returnValue);
  assert(staticAssert<MESH_NODEMAP::parallelType == CONNECT_NODEMAP::parallelType>::returnValue);
  assert(staticAssert<GEOSHAPE2D::nDim == 2>::returnValue);
  assert(staticAssert<GEOSHAPE1D::nDim == 1>::returnValue);  
  
  gridOk = false;
}

template<typename MESH, typename CONNECT>
lgmSingleList<MESH,CONNECT,2,1>::
lgmSingleList(const Teuchos::RCP<MESH>    & Grid,
              const Teuchos::RCP<CONNECT> & ConnectGrid)
{
  typedef typename MESH::GRID_GEOSHAPE  MESH_GEOSHAPE;
  typedef typename MESH::GRID_ELMAP     MESH_ELMAP;
  typedef typename MESH::GRID_NODEMAP   MESH_NODEMAP;
  
  typedef typename CONNECT::GRID_GEOSHAPE CONNECT_GEOSHAPE;
  typedef typename CONNECT::GRID_ELMAP    CONNECT_ELMAP;
  typedef typename CONNECT::GRID_NODEMAP  CONNECT_NODEMAP;
  
  assert(staticAssert<MESH_GEOSHAPE::geoName     == CONNECT_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<MESH_ELMAP::parallelType   == CONNECT_ELMAP::parallelType>::returnValue);
  assert(staticAssert<MESH_NODEMAP::parallelType == CONNECT_NODEMAP::parallelType>::returnValue);
  assert(staticAssert<GEOSHAPE2D::nDim == 2>::returnValue);
  assert(staticAssert<GEOSHAPE1D::nDim == 1>::returnValue);
  
  gridOk = true;
  
  grid        = Grid;
  connectGrid = ConnectGrid;
}

template<typename MESH, typename CONNECT>
void
lgmSingleList<MESH,CONNECT,2,1>::
setMesh(const Teuchos::RCP<MESH>    & Grid,
        const Teuchos::RCP<CONNECT> & ConnectGrid)
{
  gridOk = true;
  
  grid        = Grid;
  connectGrid = ConnectGrid;
}

template<typename MESH, typename CONNECT>
typename lgmSingleList<MESH,CONNECT,2,1>::STDMAP
lgmSingleList<MESH,CONNECT,2,1>::
getList()
{
  //Assert
  assert(gridOk);
  
  //Typedef
  typedef std::pair<POINTELEMENT,PMAPTYPE> PAIO;
  
  //Alloc
  PAIO paio;
  pointElement<GEOSHAPE1D> pElement;
  STDMAP outMap;
  UInt el;
  
  //Loop
  for(UInt ed=1; ed <= grid->getNumEdges(); ++ed)
  {
    pElement.setPoints(grid->getEdgeNodesL(ed));
    pElement.reorder();
    
    el = connectGrid->getEdgeToElement(ed,1);
    
    paio.first  = pElement;
    paio.second = grid->getElements().getRowMapL(el);
    
    outMap.insert(paio);
  }
  
  return(outMap); 
}

template<typename MESH, typename CONNECT>
sVect<point3d>
lgmSingleList<MESH,CONNECT,2,1>::
getNodes(const PMAPTYPE & pmItem)
{
  assert(gridOk);
  return(grid->getElementNodesL(pmItem.getLid()));
}


//_________________________________________________________________________________________________
// 1d Mesh - 1d Elements
//-------------------------------------------------------------------------------------------------
template<typename MESH, typename CONNECT>
class lgmSingleList<MESH,CONNECT,1,1>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename MESH::GEOSHAPE1D  GEOSHAPE1D;
    typedef typename MESH::GRID_ELMAP  PMAPTYPE;
    typedef pMap<PMAPTYPE>             PMAP;
    typedef pointElement<GEOSHAPE1D>   POINTELEMENT;
    typedef std::map<POINTELEMENT,PMAPTYPE> STDMAP;
    typedef GEOSHAPE1D                 OUTGEOSHAPE;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool gridOk;
    Teuchos::RCP<MESH>    grid;
    Teuchos::RCP<CONNECT> connectGrid;
    //@}
    
    /*! @name Functions */ //@{
  public:
    lgmSingleList();
    lgmSingleList(const Teuchos::RCP<MESH>    & Grid,
                  const Teuchos::RCP<CONNECT> & ConnectGrid);
    
    void setMesh(const Teuchos::RCP<MESH>    & Grid,
                 const Teuchos::RCP<CONNECT> & ConnectGrid);
    
    STDMAP getList();
    sVect<point3d> getNodes(const PMAPTYPE & pmItem);
    //@}
};

template<typename MESH, typename CONNECT>
lgmSingleList<MESH,CONNECT,1,1>::
lgmSingleList()
{
  typedef typename MESH::GRID_GEOSHAPE  MESH_GEOSHAPE;
  typedef typename MESH::GRID_ELMAP     MESH_ELMAP;
  typedef typename MESH::GRID_NODEMAP   MESH_NODEMAP;
  
  typedef typename CONNECT::GRID_GEOSHAPE CONNECT_GEOSHAPE;
  typedef typename CONNECT::GRID_ELMAP    CONNECT_ELMAP;
  typedef typename CONNECT::GRID_NODEMAP  CONNECT_NODEMAP;
  
  assert(staticAssert<MESH_GEOSHAPE::geoName     == CONNECT_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<MESH_ELMAP::parallelType   == CONNECT_ELMAP::parallelType>::returnValue);
  assert(staticAssert<MESH_NODEMAP::parallelType == CONNECT_NODEMAP::parallelType>::returnValue);
  assert(staticAssert<GEOSHAPE1D::nDim == 1>::returnValue);
  
  gridOk = false;
}

template<typename MESH, typename CONNECT>
lgmSingleList<MESH,CONNECT,1,1>::
lgmSingleList(const Teuchos::RCP<MESH>    & Grid,
              const Teuchos::RCP<CONNECT> & ConnectGrid)
{
  typedef typename MESH::GRID_GEOSHAPE  MESH_GEOSHAPE;
  typedef typename MESH::GRID_ELMAP     MESH_ELMAP;
  typedef typename MESH::GRID_NODEMAP   MESH_NODEMAP;
  
  typedef typename CONNECT::GRID_GEOSHAPE CONNECT_GEOSHAPE;
  typedef typename CONNECT::GRID_ELMAP    CONNECT_ELMAP;
  typedef typename CONNECT::GRID_NODEMAP  CONNECT_NODEMAP;
  
  assert(staticAssert<MESH_GEOSHAPE::geoName     == CONNECT_GEOSHAPE::geoName>::returnValue);
  assert(staticAssert<MESH_ELMAP::parallelType   == CONNECT_ELMAP::parallelType>::returnValue);
  assert(staticAssert<MESH_NODEMAP::parallelType == CONNECT_NODEMAP::parallelType>::returnValue);
  assert(staticAssert<GEOSHAPE1D::nDim == 1>::returnValue);
  
  gridOk = true;
  
  grid        = Grid;
  connectGrid = ConnectGrid;
}

template<typename MESH, typename CONNECT>
void
lgmSingleList<MESH,CONNECT,1,1>::
setMesh(const Teuchos::RCP<MESH>    & Grid,
        const Teuchos::RCP<CONNECT> & ConnectGrid)
{
  gridOk = true;
  
  grid        = Grid;
  connectGrid = ConnectGrid;
}

template<typename MESH, typename CONNECT>
typename lgmSingleList<MESH,CONNECT,1,1>::STDMAP
lgmSingleList<MESH,CONNECT,1,1>::
getList()
{
  //Assert
  assert(gridOk);
  
  //Typedef
  typedef std::pair<POINTELEMENT,PMAPTYPE> PAIO;
  
  //Alloc
  PAIO paio;
  pointElement<GEOSHAPE1D> pElement;
  STDMAP outMap;
  
  //Loop
  for(UInt el=1; el <= grid->getNumElements(); ++el)
  {
    pElement.setPoints(grid->getElementNodesL(el));
    pElement.reorder();
    
    paio.first  = pElement;
    paio.second = grid->getElements().getRowMapL(el);
    
    outMap.insert(paio);
  }
  
  return(outMap); 
}

template<typename MESH, typename CONNECT>
sVect<point3d>
lgmSingleList<MESH,CONNECT,1,1>::
getNodes(const PMAPTYPE & pmItem)
{
  assert(gridOk);
  return(grid->getElementNodesL(pmItem.getLid()));
}

#endif
