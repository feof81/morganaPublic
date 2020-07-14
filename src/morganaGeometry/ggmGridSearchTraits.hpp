/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef GGMGRIDSEARCHTRAITS_HPP
#define GGMGRIDSEARCHTRAITS_HPP

#include "lgmGridSearch.hpp"
#include "search3dA.hpp"
#include "search2dA.hpp"
#include "search1dA.hpp"


//_________________________________________________________________________________________________
// VOID IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TGT_MESH, typename SRC_MESH, UInt DIM_SRC>
class ggmGridSearchTraits
{
};


//_________________________________________________________________________________________________
// 3D - IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TGT_MESH, typename SRC_MESH>
class ggmGridSearchTraits<TGT_MESH,SRC_MESH,3> : public lgmGridSearch<TGT_MESH,SRC_MESH>
{
    /*! @name Typedefs */ //@{
  public:
    typedef lgmGridSearch<TGT_MESH,SRC_MESH>    LOCALSEARCH; 
    typedef typename LOCALSEARCH::SRC_GEOSHAPE  SRC_GEOSHAPE;
    typedef typename LOCALSEARCH::SRC_ELMAP     SRC_ELMAP;
    typedef typename LOCALSEARCH::SRC_ELMAP     SRC_NODEMAP;
    
    typedef typename SRC_GEOSHAPE::GEOBSHAPE                   SRC_SURF_GEOSHAPE;
    typedef mesh2d<SRC_SURF_GEOSHAPE,SRC_ELMAP,SRC_NODEMAP>    SRC_SURF_MESH;
    typedef connect2d<SRC_SURF_GEOSHAPE,SRC_ELMAP,SRC_NODEMAP> SRC_SURF_CONNECT;
    
    typedef typename LOCALSEARCH::TGT_CONNECT TGT_CONNECT;
    typedef typename LOCALSEARCH::SRC_CONNECT SRC_CONNECT;
    
    typedef search3dA<SRC_GEOSHAPE,SRC_ELMAP,SRC_NODEMAP> POINTSEARCH;
    //@}
    
    /*! @name Links and internal structures */ //@{
  public:
    POINTSEARCH pointSearch;
    
    bool commDevOk;
    Teuchos::RCP<const communicator> commDev;
    
    bool gridOk;
    Teuchos::RCP<TGT_MESH>    tgtGrid;
    Teuchos::RCP<SRC_MESH>    srcGrid;
    Teuchos::RCP<TGT_CONNECT> tgtConnect;
    Teuchos::RCP<SRC_CONNECT> srcConnect;
    //@}
    
    /*! @name Functions */ //@{
  public:
    ggmGridSearchTraits();
    ggmGridSearchTraits(const Teuchos::RCP<const communicator> & CommDev);
    ggmGridSearchTraits(const communicator & CommDev);
    
    void setCommDev(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    
    void setMesh(const Teuchos::RCP<TGT_MESH>         & TgtGrid,
                 const Teuchos::RCP<TGT_CONNECT>      & TgtConnect,
                 const Teuchos::RCP<SRC_MESH>         & SrcGrid,
                 const Teuchos::RCP<SRC_CONNECT>      & SrcConnect,
                 const Teuchos::RCP<SRC_SURF_MESH>    & SrcSurfGrid,
                 const Teuchos::RCP<SRC_SURF_CONNECT> & SrcSurfConnect);
    //@}
};

template<typename TGT_MESH, typename SRC_MESH>
ggmGridSearchTraits<TGT_MESH,SRC_MESH,3>::
ggmGridSearchTraits() : lgmGridSearch<TGT_MESH,SRC_MESH>::lgmGridSearch()
{
  commDevOk = false;
  gridOk    = false;
}

template<typename TGT_MESH, typename SRC_MESH>
ggmGridSearchTraits<TGT_MESH,SRC_MESH,3>::
ggmGridSearchTraits(const Teuchos::RCP<const communicator> & CommDev) : lgmGridSearch<TGT_MESH,SRC_MESH>::lgmGridSearch()
{
  commDevOk = true;
  gridOk    = false;
  commDev   = CommDev;
  
  pointSearch.setCommunicator(CommDev);
}

template<typename TGT_MESH, typename SRC_MESH>
ggmGridSearchTraits<TGT_MESH,SRC_MESH,3>::
ggmGridSearchTraits(const communicator & CommDev)
{
  commDevOk = true;
  gridOk    = false;
  commDev   = Teuchos::rcpFromRef(CommDev);
  
  pointSearch.setCommunicator(CommDev);
}

template<typename TGT_MESH, typename SRC_MESH>
void
ggmGridSearchTraits<TGT_MESH,SRC_MESH,3>::
setCommDev(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevOk = true;
  commDev   = CommDev;
  
  pointSearch.setCommunicator(CommDev);
}

template<typename TGT_MESH, typename SRC_MESH>
void
ggmGridSearchTraits<TGT_MESH,SRC_MESH,3>::
setCommunicator(const communicator & CommDev)
{
  commDevOk = true;
  commDev   = Teuchos::rcpFromRef(CommDev);
  
  pointSearch.setCommunicator(CommDev);
}

template<typename TGT_MESH, typename SRC_MESH>
void
ggmGridSearchTraits<TGT_MESH,SRC_MESH,3>::
setMesh(const Teuchos::RCP<TGT_MESH>         & TgtGrid,
        const Teuchos::RCP<TGT_CONNECT>      & TgtConnect,
        const Teuchos::RCP<SRC_MESH>         & SrcGrid,
        const Teuchos::RCP<SRC_CONNECT>      & SrcConnect,
        const Teuchos::RCP<SRC_SURF_MESH>    & SrcSurfGrid,
        const Teuchos::RCP<SRC_SURF_CONNECT> & SrcSurfConnect)
{
  assert(commDevOk);
  
  gridOk = true;
  tgtGrid    = TgtGrid;
  srcGrid    = SrcGrid;
  tgtConnect = TgtConnect;
  srcConnect = SrcConnect;
  
  LOCALSEARCH::setMesh(TgtGrid,
                       TgtConnect,
                       SrcGrid,
                       SrcConnect);
  
  pointSearch.setMesh(SrcSurfGrid,
                      SrcGrid,
                      SrcSurfConnect,
                      SrcConnect);
  
  pointSearch.localInit();
  pointSearch.globalInit();
}


//_________________________________________________________________________________________________
// 2D - IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TGT_MESH, typename SRC_MESH>
class ggmGridSearchTraits<TGT_MESH,SRC_MESH,2> : public lgmGridSearch<TGT_MESH,SRC_MESH>
{
    /*! @name Typedefs */ //@{
  public:
    typedef lgmGridSearch<TGT_MESH,SRC_MESH>    LOCALSEARCH; 
    typedef typename LOCALSEARCH::SRC_GEOSHAPE  SRC_GEOSHAPE;
    typedef typename LOCALSEARCH::SRC_ELMAP     SRC_ELMAP;
    typedef typename LOCALSEARCH::SRC_ELMAP     SRC_NODEMAP;
    
    typedef typename LOCALSEARCH::TGT_CONNECT TGT_CONNECT;
    typedef typename LOCALSEARCH::SRC_CONNECT SRC_CONNECT;
    
    typedef search2dA<SRC_GEOSHAPE,SRC_ELMAP,SRC_NODEMAP> POINTSEARCH;
    //@}
    
    /*! @name Links and internal structures */ //@{
  public:
    POINTSEARCH pointSearch;
    
    bool commDevOk;
    Teuchos::RCP<const communicator> commDev;
    
    bool gridOk;
    Teuchos::RCP<TGT_MESH>    tgtGrid;
    Teuchos::RCP<SRC_MESH>    srcGrid;
    Teuchos::RCP<TGT_CONNECT> tgtConnect;
    Teuchos::RCP<SRC_CONNECT> srcConnect;
    //@}
    
    /*! @name Functions */ //@{
  public:
    ggmGridSearchTraits();
    ggmGridSearchTraits(const Teuchos::RCP<const communicator> & CommDev);
    ggmGridSearchTraits(const communicator & CommDev);
    
    void setCommDev(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    
    void setMesh(const Teuchos::RCP<TGT_MESH>         & TgtGrid,
                 const Teuchos::RCP<TGT_CONNECT>      & TgtConnect,
                 const Teuchos::RCP<SRC_MESH>         & SrcGrid,
                 const Teuchos::RCP<SRC_CONNECT>      & SrcConnect);
    //@}
};

template<typename TGT_MESH, typename SRC_MESH>
ggmGridSearchTraits<TGT_MESH,SRC_MESH,2>::
ggmGridSearchTraits() : lgmGridSearch<TGT_MESH,SRC_MESH>::lgmGridSearch()
{
  commDevOk = false;
  gridOk    = false;
}

template<typename TGT_MESH, typename SRC_MESH>
ggmGridSearchTraits<TGT_MESH,SRC_MESH,2>::
ggmGridSearchTraits(const Teuchos::RCP<const communicator> & CommDev) : lgmGridSearch<TGT_MESH,SRC_MESH>::lgmGridSearch()
{
  commDevOk = true;
  gridOk    = false;
  commDev   = CommDev;
  
  pointSearch.setCommunicator(CommDev);
}

template<typename TGT_MESH, typename SRC_MESH>
ggmGridSearchTraits<TGT_MESH,SRC_MESH,2>::
ggmGridSearchTraits(const communicator & CommDev)
{
  commDevOk = true;
  gridOk    = false;
  commDev   = Teuchos::rcpFromRef(CommDev);
  
  pointSearch.setCommunicator(CommDev);
}

template<typename TGT_MESH, typename SRC_MESH>
void
ggmGridSearchTraits<TGT_MESH,SRC_MESH,2>::
setCommDev(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevOk = true;
  gridOk    = false;
  commDev   = CommDev;
  
  pointSearch.setCommunicator(CommDev);
}

template<typename TGT_MESH, typename SRC_MESH>
void
ggmGridSearchTraits<TGT_MESH,SRC_MESH,2>::
setCommunicator(const communicator & CommDev)
{
  commDevOk = true;
  commDev   = Teuchos::rcpFromRef(CommDev);
  
  pointSearch.setCommunicator(CommDev);
}

template<typename TGT_MESH, typename SRC_MESH>
void
ggmGridSearchTraits<TGT_MESH,SRC_MESH,2>::
setMesh(const Teuchos::RCP<TGT_MESH>         & TgtGrid,
        const Teuchos::RCP<TGT_CONNECT>      & TgtConnect,
        const Teuchos::RCP<SRC_MESH>         & SrcGrid,
        const Teuchos::RCP<SRC_CONNECT>      & SrcConnect)
{
  assert(commDevOk);
  
  gridOk = true;
  tgtGrid    = TgtGrid;
  srcGrid    = SrcGrid;
  tgtConnect = TgtConnect;
  srcConnect = SrcConnect;
  
  LOCALSEARCH::setMesh(TgtGrid,
                       TgtConnect,
                       SrcGrid,
                       SrcConnect);
  
  pointSearch.setMesh(SrcGrid,
                      SrcConnect);
  
  pointSearch.localInit();
  pointSearch.globalInit();
}


//_________________________________________________________________________________________________
// 1D - IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TGT_MESH, typename SRC_MESH>
class ggmGridSearchTraits<TGT_MESH,SRC_MESH,1> : public lgmGridSearch<TGT_MESH,SRC_MESH>
{
    /*! @name Typedefs */ //@{
  public:
    typedef lgmGridSearch<TGT_MESH,SRC_MESH>    LOCALSEARCH; 
    typedef typename LOCALSEARCH::SRC_GEOSHAPE  SRC_GEOSHAPE;
    typedef typename LOCALSEARCH::SRC_ELMAP     SRC_ELMAP;
    typedef typename LOCALSEARCH::SRC_ELMAP     SRC_NODEMAP;
    
    typedef typename LOCALSEARCH::TGT_CONNECT TGT_CONNECT;
    typedef typename LOCALSEARCH::SRC_CONNECT SRC_CONNECT;
    
    typedef search1dA<SRC_GEOSHAPE,SRC_ELMAP,SRC_NODEMAP> POINTSEARCH;
    //@}
    
    /*! @name Links and internal structures */ //@{
  public:
    POINTSEARCH pointSearch;
    
    bool commDevOk;
    Teuchos::RCP<const communicator> commDev;
    
    bool gridOk;
    Teuchos::RCP<TGT_MESH>    tgtGrid;
    Teuchos::RCP<SRC_MESH>    srcGrid;
    Teuchos::RCP<TGT_CONNECT> tgtConnect;
    Teuchos::RCP<SRC_CONNECT> srcConnect;
    //@}
    
    /*! @name Functions */ //@{
  public:
    ggmGridSearchTraits();
    ggmGridSearchTraits(const Teuchos::RCP<const communicator> & CommDev);
    ggmGridSearchTraits(const communicator & CommDev);
    
    void setCommDev(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    
    void setMesh(const Teuchos::RCP<TGT_MESH>         & TgtGrid,
                 const Teuchos::RCP<TGT_CONNECT>      & TgtConnect,
                 const Teuchos::RCP<SRC_MESH>         & SrcGrid,
                 const Teuchos::RCP<SRC_CONNECT>      & SrcConnect);
    //@}
};

template<typename TGT_MESH, typename SRC_MESH>
ggmGridSearchTraits<TGT_MESH,SRC_MESH,1>::
ggmGridSearchTraits() : lgmGridSearch<TGT_MESH,SRC_MESH>::lgmGridSearch()
{
  commDevOk = false;
  gridOk    = false;
}

template<typename TGT_MESH, typename SRC_MESH>
ggmGridSearchTraits<TGT_MESH,SRC_MESH,1>::
ggmGridSearchTraits(const Teuchos::RCP<const communicator> & CommDev) : lgmGridSearch<TGT_MESH,SRC_MESH>::lgmGridSearch()
{
  commDevOk = true;
  gridOk    = false;
  commDev   = CommDev;
  
  pointSearch.setCommunicator(CommDev);
}

template<typename TGT_MESH, typename SRC_MESH>
ggmGridSearchTraits<TGT_MESH,SRC_MESH,1>::
ggmGridSearchTraits(const communicator & CommDev)
{
  commDevOk = true;
  gridOk    = false;
  commDev   = Teuchos::rcpFromRef(CommDev);
  
  pointSearch.setCommunicator(CommDev);
}

template<typename TGT_MESH, typename SRC_MESH>
void
ggmGridSearchTraits<TGT_MESH,SRC_MESH,1>::
setCommDev(const Teuchos::RCP<const communicator> & CommDev)
{
  commDevOk = true;
  commDev   = CommDev;
  
  pointSearch.setCommunicator(CommDev);
}

template<typename TGT_MESH, typename SRC_MESH>
void
ggmGridSearchTraits<TGT_MESH,SRC_MESH,1>::
setCommunicator(const communicator & CommDev)
{
  commDevOk = true;
  commDev   = Teuchos::rcpFromRef(CommDev);
  
  pointSearch.setCommunicator(CommDev);
}

template<typename TGT_MESH, typename SRC_MESH>
void
ggmGridSearchTraits<TGT_MESH,SRC_MESH,1>::
setMesh(const Teuchos::RCP<TGT_MESH>         & TgtGrid,
        const Teuchos::RCP<TGT_CONNECT>      & TgtConnect,
        const Teuchos::RCP<SRC_MESH>         & SrcGrid,
        const Teuchos::RCP<SRC_CONNECT>      & SrcConnect)
{
  assert(commDevOk);
  
  gridOk = true;
  tgtGrid    = TgtGrid;
  srcGrid    = SrcGrid;
  tgtConnect = TgtConnect;
  srcConnect = SrcConnect;
  
  LOCALSEARCH::setMesh(TgtGrid,
                       TgtConnect,
                       SrcGrid,
                       SrcConnect);
  
  pointSearch.setMesh(SrcGrid,
                      SrcConnect);
  
  pointSearch.localInit();
  pointSearch.globalInit();
}

#endif
