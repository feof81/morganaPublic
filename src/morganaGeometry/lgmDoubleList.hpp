/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef LGMDOUBLELIST_HPP
#define LGMDOUBLELIST_HPP

#include "lgmSingleList.hpp"
#include "lgmElementTraits.hpp"

template<typename TGT_MESH, typename SRC_MESH>
class lgmDoubleList
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename TGT_MESH::GRID_GEOSHAPE  TGT_GEOSHAPE;
    typedef typename SRC_MESH::GRID_GEOSHAPE  SRC_GEOSHAPE;
    typedef typename TGT_MESH::GRID_ELMAP     TGT_ELMAP;
    typedef typename SRC_MESH::GRID_ELMAP     SRC_ELMAP;
    typedef typename TGT_MESH::GRID_NODEMAP   TGT_NODEMAP;
    typedef typename SRC_MESH::GRID_NODEMAP   SRC_NODEMAP;
    
    typedef lgmElementTraits<TGT_GEOSHAPE, TGT_ELMAP, TGT_NODEMAP, TGT_GEOSHAPE::nDim> TGT_TRAITS;
    typedef lgmElementTraits<SRC_GEOSHAPE, SRC_ELMAP, SRC_NODEMAP, SRC_GEOSHAPE::nDim> SRC_TRAITS;
    
    typedef typename TGT_TRAITS::CONNECT TGT_CONNECT;
    typedef typename SRC_TRAITS::CONNECT SRC_CONNECT;
    
    typedef lgmSingleList<TGT_MESH, TGT_CONNECT, TGT_GEOSHAPE::nDim, TGT_GEOSHAPE::nDim> TGT_SINGLELIST;
    typedef lgmSingleList<SRC_MESH, SRC_CONNECT, SRC_GEOSHAPE::nDim, TGT_GEOSHAPE::nDim> SRC_SINGLELIST;
    
    typedef typename TGT_SINGLELIST::STDMAP      TGT_STDMAP;
    typedef typename SRC_SINGLELIST::STDMAP      SRC_STDMAP;
    typedef typename TGT_SINGLELIST::OUTGEOSHAPE TGT_OUTGEOSHAPE;
    typedef typename SRC_SINGLELIST::OUTGEOSHAPE SRC_OUTGEOSHAPE;
    //@}
    
    /*! @name Internal structures */ //@{
  public:
    TGT_SINGLELIST tgtListGen;
    SRC_SINGLELIST srcListGen;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool gridOk;
    Teuchos::RCP<TGT_MESH>    tgtGrid;
    Teuchos::RCP<SRC_MESH>    srcGrid;
    Teuchos::RCP<TGT_CONNECT> tgtConnect;
    Teuchos::RCP<SRC_CONNECT> srcConnect;
    //@}
    
    /*! @name Functions */ //@{
  public:
    lgmDoubleList();
    lgmDoubleList(const Teuchos::RCP<TGT_MESH>    & TgtGrid,
                  const Teuchos::RCP<TGT_CONNECT> & TgtConnect,
                  const Teuchos::RCP<SRC_MESH>    & SrcGrid,
                  const Teuchos::RCP<SRC_CONNECT> & SrcConnect);
    
    void setMesh(const Teuchos::RCP<TGT_MESH>    & TgtGrid,
                 const Teuchos::RCP<TGT_CONNECT> & TgtConnect,
                 const Teuchos::RCP<SRC_MESH>    & SrcGrid,
                 const Teuchos::RCP<SRC_CONNECT> & SrcConnect);
    
    TGT_STDMAP getTgtList();
    SRC_STDMAP getSrcList();
    
    sVect<point3d> getTgtNodes(const TGT_ELMAP & pmItem);
    sVect<point3d> getSrcNodes(const SRC_ELMAP & pmItem);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TGT_MESH, typename SRC_MESH>
lgmDoubleList<TGT_MESH,SRC_MESH>::
lgmDoubleList()
{
  assert(staticAssert<TGT_OUTGEOSHAPE::nDim    == SRC_OUTGEOSHAPE::nDim>::returnValue);
  assert(staticAssert<TGT_OUTGEOSHAPE::geoName == SRC_OUTGEOSHAPE::geoName>::returnValue);
  
  gridOk = false;
}

template<typename TGT_MESH, typename SRC_MESH>
lgmDoubleList<TGT_MESH,SRC_MESH>::
lgmDoubleList(const Teuchos::RCP<TGT_MESH>    & TgtGrid,
              const Teuchos::RCP<TGT_CONNECT> & TgtConnect,
              const Teuchos::RCP<SRC_MESH>    & SrcGrid,
              const Teuchos::RCP<SRC_CONNECT> & SrcConnect)
{
  assert(staticAssert<TGT_OUTGEOSHAPE::nDim    == SRC_OUTGEOSHAPE::nDim>::returnValue);
  assert(staticAssert<TGT_OUTGEOSHAPE::geoName == SRC_OUTGEOSHAPE::geoName>::returnValue);
  
  gridOk = true;
  
  tgtGrid    = TgtGrid;
  tgtConnect = TgtConnect;
  srcGrid    = SrcGrid;
  srcConnect = SrcConnect;
  
  tgtListGen.setMesh(TgtGrid,TgtConnect);
  srcListGen.setMesh(SrcGrid,SrcConnect);
}

template<typename TGT_MESH, typename SRC_MESH>
void
lgmDoubleList<TGT_MESH,SRC_MESH>::
setMesh(const Teuchos::RCP<TGT_MESH>    & TgtGrid,
        const Teuchos::RCP<TGT_CONNECT> & TgtConnect,
        const Teuchos::RCP<SRC_MESH>    & SrcGrid,
        const Teuchos::RCP<SRC_CONNECT> & SrcConnect)
{
  gridOk = true;
  
  tgtGrid    = TgtGrid;
  tgtConnect = TgtConnect;
  srcGrid    = SrcGrid;
  srcConnect = SrcConnect;
  
  tgtListGen.setMesh(TgtGrid,TgtConnect);
  srcListGen.setMesh(SrcGrid,SrcConnect);
}
    
template<typename TGT_MESH, typename SRC_MESH>
typename lgmDoubleList<TGT_MESH,SRC_MESH>::TGT_STDMAP
lgmDoubleList<TGT_MESH,SRC_MESH>::
getTgtList()
{
  assert(gridOk);
  
  return(tgtListGen.getList());
}
    
template<typename TGT_MESH, typename SRC_MESH>
typename lgmDoubleList<TGT_MESH,SRC_MESH>::SRC_STDMAP
lgmDoubleList<TGT_MESH,SRC_MESH>::
getSrcList()
{
  assert(gridOk);
  
  return(srcListGen.getList());
}

template<typename TGT_MESH, typename SRC_MESH>
sVect<point3d>
lgmDoubleList<TGT_MESH,SRC_MESH>::
getTgtNodes(const TGT_ELMAP & pmItem)
{
  return(tgtListGen.getNodes(pmItem));
}

template<typename TGT_MESH, typename SRC_MESH>
sVect<point3d>
lgmDoubleList<TGT_MESH,SRC_MESH>::
getSrcNodes(const SRC_ELMAP & pmItem)
{
  return(srcListGen.getNodes(pmItem));
}


#endif
