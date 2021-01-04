/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef LGMGRIDSEARCH_HPP
#define LGMGRIDSEARCH_HPP

#include "lgmDoubleList.hpp"
#include "lgmSearchData.hpp"


template<typename TGT_MESH, typename SRC_MESH>
class lgmGridSearch
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename TGT_MESH::GRID_GEOSHAPE  TGT_GEOSHAPE;
    typedef typename SRC_MESH::GRID_GEOSHAPE  SRC_GEOSHAPE;
    typedef typename TGT_MESH::GRID_ELMAP     TGT_ELMAP;
    typedef typename SRC_MESH::GRID_ELMAP     SRC_ELMAP;
    typedef typename TGT_MESH::GRID_ELMAP     TGT_NODEMAP;
    typedef typename SRC_MESH::GRID_ELMAP     SRC_NODEMAP;
    
    typedef pointElement<TGT_GEOSHAPE> INDATA;
    typedef lgmSearchData<SRC_ELMAP>   OUTDATA;
    
    typedef pVect<INDATA,TGT_ELMAP>    INVECT;
    typedef pVect<OUTDATA,TGT_ELMAP>   OUTVECT;
    
    typedef lgmDoubleList<TGT_MESH,SRC_MESH> DOUBLELIST;
    typedef lgmElementTraits<TGT_GEOSHAPE,TGT_ELMAP,TGT_NODEMAP,TGT_GEOSHAPE::nDim> TGT_TRAITS;
    typedef lgmElementTraits<SRC_GEOSHAPE,SRC_ELMAP,SRC_NODEMAP,SRC_GEOSHAPE::nDim> SRC_TRATTS;
    
    typedef typename TGT_TRAITS::GEOSUPPORT TGT_PROJECTOR;
    typedef typename SRC_TRATTS::GEOSUPPORT SRC_PROJECTOR;
    
    typedef typename DOUBLELIST::TGT_STDMAP  TGT_STDMAP;
    typedef typename DOUBLELIST::SRC_STDMAP  SRC_STDMAP;
    typedef typename DOUBLELIST::TGT_CONNECT TGT_CONNECT;
    typedef typename DOUBLELIST::SRC_CONNECT SRC_CONNECT;
    //@}
    
    /*! @name Internal structures */ //@{
  public:
    DOUBLELIST    doubleList;
    TGT_PROJECTOR tgtProjector;
    SRC_PROJECTOR srcProjector;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool gridOk;
    Teuchos::RCP<TGT_MESH> tgtGrid;
    Teuchos::RCP<SRC_MESH> srcGrid;
    //@}
    
    /*! @name Functions */ //@{
  public:
    lgmGridSearch();
    lgmGridSearch(const Teuchos::RCP<TGT_MESH>    & TgtGrid,
                  const Teuchos::RCP<TGT_CONNECT> & TgtConnect,
                  const Teuchos::RCP<SRC_MESH>    & SrcGrid,
                  const Teuchos::RCP<SRC_CONNECT> & SrcConnect);
    
    void setMesh(const Teuchos::RCP<TGT_MESH>    & TgtGrid,
                 const Teuchos::RCP<TGT_CONNECT> & TgtConnect,
                 const Teuchos::RCP<SRC_MESH>    & SrcGrid,
                 const Teuchos::RCP<SRC_CONNECT> & SrcConnect);
    
    OUTVECT localSearch(const sVect<point3d> & locCoords,
                        const INVECT         & inVect);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TGT_MESH, typename SRC_MESH>
lgmGridSearch<TGT_MESH,SRC_MESH>::
lgmGridSearch()
{
  assert(staticAssert<TGT_ELMAP::parallelType   == SRC_ELMAP::parallelType>::returnValue);
  assert(staticAssert<TGT_NODEMAP::parallelType == SRC_NODEMAP::parallelType>::returnValue);
  
  gridOk = false;
}

template<typename TGT_MESH, typename SRC_MESH>
lgmGridSearch<TGT_MESH,SRC_MESH>::
lgmGridSearch(const Teuchos::RCP<TGT_MESH>    & TgtGrid,
              const Teuchos::RCP<TGT_CONNECT> & TgtConnect,
              const Teuchos::RCP<SRC_MESH>    & SrcGrid,
              const Teuchos::RCP<SRC_CONNECT> & SrcConnect)
{
  assert(staticAssert<TGT_ELMAP::parallelType   == SRC_ELMAP::parallelType>::returnValue);
  assert(staticAssert<TGT_NODEMAP::parallelType == SRC_NODEMAP::parallelType>::returnValue);
  
  gridOk = true;
  
  tgtGrid = TgtGrid;
  srcGrid = SrcGrid;
  
  doubleList.setMesh(TgtGrid,
                     TgtConnect,
                     SrcGrid,
                     SrcConnect);
}

template<typename TGT_MESH, typename SRC_MESH>
void
lgmGridSearch<TGT_MESH,SRC_MESH>::
setMesh(const Teuchos::RCP<TGT_MESH>    & TgtGrid,
        const Teuchos::RCP<TGT_CONNECT> & TgtConnect,
        const Teuchos::RCP<SRC_MESH>    & SrcGrid,
        const Teuchos::RCP<SRC_CONNECT> & SrcConnect)

{
  gridOk = true;
  
  tgtGrid = TgtGrid;
  srcGrid = SrcGrid;
  
  doubleList.setMesh(TgtGrid,
                     TgtConnect,
                     SrcGrid,
                     SrcConnect);
}

template<typename TGT_MESH, typename SRC_MESH>
typename lgmGridSearch<TGT_MESH,SRC_MESH>::OUTVECT
lgmGridSearch<TGT_MESH,SRC_MESH>::
localSearch(const sVect<point3d> & inLocCoords,
            const INVECT         & inVect)
{
  //Assert
  assert(gridOk);
  
  //Typedefs
  typedef typename SRC_STDMAP::iterator ITER;
  
  //Alloc
  OUTVECT        outVect;
  OUTDATA        outData;
  SRC_ELMAP      outMap;
  sVect<point3d> outNodes;
  sVect<point3d> outLocCoords(inLocCoords.size());
  INDATA         inElement;
  
  ITER  iter;
  bool  found;
  point3d P, Y;
  
  //Main loop
  if(TGT_GEOSHAPE::nDim <= SRC_GEOSHAPE::nDim)
  {
    //Create local list
    SRC_STDMAP srcStdMap = doubleList.getSrcList();
    
    //Loop
    for(UInt i=1; i <= inVect.size(); ++i)
    {
      //Element to be searched
      inElement = inVect(i);
      inElement.reorder();
      
      //Search 
      found = (srcStdMap.count(inElement) >= 1);
      
      if(found)
      {
        iter   = srcStdMap.find(inElement);
        outMap = iter->second;

        outNodes = doubleList.getSrcNodes(outMap);
        srcProjector.setPoints(outNodes);

        inElement = inVect(i);

        for(UInt j=1; j <= inLocCoords.size(); ++j)
        {
          P     = tgtProjector.getPosition(inElement.getPoints(),inLocCoords(j));
          found = found && srcProjector.projection(P,Y);
          outLocCoords(j) = Y;
        }
        
        outData.setFound(found);
        outData.setElMap(outMap);
        outData.setLocCoords(outLocCoords);
      }
      else
      {
        outData.setFound(false);
        outData.setElMap(SRC_ELMAP());
        outData.clearLocCoords();
      }
      
      //Upload data 
      outVect.push_back(inVect.getMapL(i),outData);
    }
  }
  
  //Update finder
  outVect.updateFinder();
  
  //Return
  return(outVect);
}

#endif
