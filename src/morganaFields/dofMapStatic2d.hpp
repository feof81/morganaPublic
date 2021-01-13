/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef DOFMAPSTATIC2D_HPP
#define DOFMAPSTATIC2D_HPP

enum dms2d_order {dms2d_vectMajor, dms2d_componentMajor};
enum dms2d_mode  {dms2d_allMode, dms2d_geoIdMode, dms2d_disjointMode};

#include "typesInterface.hpp"
#include "traitsBasic.h"

#include "pVectManip.hpp"

#include "morganaGeometry.hpp"
#include "geoMapInterface.hpp"
#include "connect2d.hpp"
#include "mesh2dGlobalManip.hpp"

#include "dofMapStatic2d_options.h"
#include "feStaticDofCard2d.h"


//! Forward declaration
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER = dms2d_vectMajor, dms2d_mode MODE = dms2d_allMode> class dofMapStatic2d;



//_________________________________________________________________________________________________
// ALL MODE
//-------------------------------------------------------------------------------------------------

/*! Support for the mapping of 2d static finite elements -> \c allMode */
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
class dofMapStatic2d<FETYPE,DOFTYPE, ORDER,dms2d_allMode>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE   PMAPTYPE;
    typedef typename FETYPE::GEOSHAPE          GEOSHAPE;
    typedef typename FETYPE::DOFCARD           DOFCARD;
    typedef typename FETYPE::FECARD            FECARD;
    
    typedef mesh2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH2D;
    typedef connect2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT2D;
    typedef geoMapInterface<GEOSHAPE>             GEOMAPINTERFACE;
    typedef pVect<FECARD,PMAPTYPE>                FECARDS;
    typedef dofMapStatic2d_options                OPTIONS;
    
    typedef pMap<PMAPTYPE> DOFMAP;
    typedef pMap<PMAPTYPE> LISTMAP;
    //@}
    
    /*! @name Internal structures */ //@{
  public:
    FETYPE          refFe;
    GEOMAPINTERFACE refShape;
    FECARDS         feCards;
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<communicator>  commDev;
    Teuchos::RCP<MESH2D>        grid2d;
    Teuchos::RCP<CONNECT2D>     connectGrid2d;
    Teuchos::RCP<OPTIONS>       options;
    //@}
    
    /*! @name Internal flags */ //@{
  public:
    bool commDevLoaded, geometryLoaded, optionsLoaded, startupOk;
    //@}
    
     /*! @name Static data */ //@{
  public:
    static const UInt block = traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
    //@}
    
    /*! @name Internal Numbers */ //@{
  public:
    UInt numVerticesL, numVerticesG;
    UInt numEdgesL,    numEdgesG;
    UInt numElementsL, numElementsG;
    
    UInt numVertexDofsL, numVertexDofsG;
    UInt numEdgesDofsL,  numEdgesDofsG;
    UInt numVolumeDofsL, numVolumeDofsG;
    
    UInt numDofsL, numDofsG;
    UInt sizeListL, sizeListG;
    //@}
    
    /*! @name Internal Maps */ //@{
  public:
    DOFMAP  dofMap;
    LISTMAP listMap;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    dofMapStatic2d();
    dofMapStatic2d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d);
    dofMapStatic2d(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d);
    dofMapStatic2d(const dofMapStatic2d & DofMap);
    dofMapStatic2d operator=(const dofMapStatic2d & DofMap);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d);
    void setGeometry(MESH2D & Grid2d, CONNECT2D & ConnedGrid2d);
    //@}
    
    /*! @name Options loading and startup */ //@{
  public:
    void setFeCardL(const UInt & lid, const FECARD & FeCards);
    void setFeCardG(const UInt & gid, const FECARD & FeCards);
    void setFeCards(const FECARDS & FECards);
    void setOptions(const Teuchos::RCP<OPTIONS> & Options);
    void setOptions(OPTIONS & Options);
    void startup();
    //@}
    
    /*! @name Get functions */ //@{
  public:
    bool isActive(const DOFCARD & dofCard) const;
    UInt dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const;
    UInt dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const;
    UInt mapDofL(const DOFCARD & dofCard) const;
    UInt mapDofG(const DOFCARD & dofCard) const;    
    UInt mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const;
    UInt mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const;
    const UInt & getNumDofsL() const;
    const UInt & getNumDofsG() const;
    const UInt & getSizeListL() const;
    const UInt & getSizeListG() const;
    //@}
    
     /*! @name Dump functions */ //@{
  public:
    const DOFMAP                     & getDofMap() const;
    const LISTMAP                    & getListMap() const;
    const FECARDS                    & getFeCards() const;
          FECARDS                    & getFeCards();
    const Teuchos::RCP<communicator> & getCommDev() const;
    const Teuchos::RCP<MESH2D>       & getGrid2d() const;
    const Teuchos::RCP<CONNECT2D>    & getConnectGrid2d() const;
    const Teuchos::RCP<OPTIONS>      & getOptions() const;
    const bool                       & getStartupOk() const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    template<typename F, typename D, dms2d_order O>
    friend ostream & operator<<(ostream & f, const dofMapStatic2d<F,D, O,dms2d_allMode> & V);
    //@}
};


template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
dofMapStatic2d()
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  
  commDevLoaded  = false;
  geometryLoaded = false;
  optionsLoaded  = false;
  startupOk      = false;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
dofMapStatic2d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = CommDev;
  grid2d        = Grid2d;
  connectGrid2d = ConnedGrid2d;
  
  feCards.resize(grid2d->getNumElements());
  feCards.setMap(grid2d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
dofMapStatic2d(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = Teuchos::rcpFromRef(CommDev);
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnedGrid2d));
  
  feCards.resize(grid2d->getNumElements());
  feCards.setMap(grid2d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
dofMapStatic2d(const dofMapStatic2d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid2d        = DofMap.grid2d;
  connectGrid2d = DofMap.connectGrid2d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;

  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numEdgesL    = DofMap.numEdgesL;
  numEdgesG    = DofMap.numEdgesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numEdgesDofsL  = DofMap.numEdgesDofsL;
  numEdgesDofsG  = DofMap.numEdgesDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  sizeListL = DofMap.sizeListL;
  sizeListG = DofMap.sizeListG;
    
  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
operator=(const dofMapStatic2d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid2d        = DofMap.grid2d;
  connectGrid2d = DofMap.connectGrid2d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;

  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numEdgesL    = DofMap.numEdgesL;
  numEdgesG    = DofMap.numEdgesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numEdgesDofsL  = DofMap.numEdgesDofsL;
  numEdgesDofsG  = DofMap.numEdgesDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  sizeListL = DofMap.sizeListL;
  sizeListG = DofMap.sizeListG;
    
  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
  
  return(*this);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded  = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded  = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d)
{
  geometryLoaded = true;
  grid2d         = Grid2d;
  connectGrid2d  = ConnedGrid2d;
  
  feCards.resize(grid2d->getNumElements());
  feCards.setMap(grid2d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
setGeometry(MESH2D & Grid2d, CONNECT2D & ConnedGrid2d)
{
  geometryLoaded = true; 
  grid2d         = Teuchos::rcp(new MESH2D(Grid2d));
  connectGrid2d  = Teuchos::rcp(new CONNECT2D(ConnedGrid2d));
  
  feCards.resize(grid2d->getNumElements());
  feCards.setMap(grid2d->getElements().getRowMap());
  feCards.updateFinder();
}



//_________________________________________________________________________________________________
// OPTIONS LOADING AND STARTUP
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
setFeCardL(const UInt & lid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(lid <= feCards.sizeL());
  feCards.getL(lid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
setFeCardG(const UInt & gid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(feCards.isG(gid));
  feCards.getG(gid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
setFeCards(const FECARDS & FECards)
{
  assert(geometryLoaded);
  assert(feCards.size() == FECards.size());
  feCards = FECards;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
setOptions(const Teuchos::RCP<OPTIONS> & Options)
{
  optionsLoaded = true;
  options = Options;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
setOptions(OPTIONS & Options)
{
  optionsLoaded = true;
  options = Teuchos::rcp(new OPTIONS(Options));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
startup()
{
  //Assert and flags_______________________________________________________________________________
  assert(commDevLoaded);
  assert(geometryLoaded);
  assert(optionsLoaded);
  
  assert(grid2d->getElements().colIsLocal());
  assert(grid2d->getEdges().colIsLocal());
  
  mesh2dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> manipulator(commDev);
  assert(manipulator.check(grid2d));
  
  startupOk = true;
  
  
  //Count the geotypes_____________________________________________________________________________ 
  numVerticesL = grid2d->getNumVertices();
  numEdgesL    = grid2d->getNumEdges();
  numElementsL = grid2d->getNumElements();
  
  numVerticesG = manipulator.getNumGlobalVertices(grid2d);
  numEdgesG    = manipulator.getNumGlobalEdges(grid2d);
  numElementsG = manipulator.getNumGlobalElements(grid2d);
  
  //Count the dofs
  numVertexDofsL = FETYPE::dofPerVertex * numVerticesL;
  numEdgesDofsL  = FETYPE::dofPerEdge   * numEdgesL;
  numVolumeDofsL = FETYPE::dofPerVolume * numElementsL;
  
  numVertexDofsG = FETYPE::dofPerVertex * numVerticesG;
  numEdgesDofsG  = FETYPE::dofPerEdge   * numEdgesG;
  numVolumeDofsG = FETYPE::dofPerVolume * numElementsG;
  
  
  //Total number of dofs and lists_________________________________________________________________
  numDofsL = numVertexDofsL + numEdgesDofsL + numVolumeDofsL;
  numDofsG = numVertexDofsG + numEdgesDofsG + numVolumeDofsG;
  
  sizeListL = numDofsL * traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
  sizeListG = numDofsG * traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
  
  
  //Build dof map__________________________________________________________________________________
  PMAPTYPE mapItem;
  UInt k=1;
  
  dofMap.resize(numDofsL);
  
  for(UInt lev=1; lev <= FETYPE::dofPerVertex; ++lev)  //Vertices
  {
    for(UInt i=1; i <= numVerticesL; ++i)
    {
      mapItem = grid2d->getNodes().getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numVerticesG) );
      
      dofMap(k) = mapItem;
      ++k;
    } 
  }
  
  for(UInt lev=1; lev <= FETYPE::dofPerEdge; ++lev)  //Edges
  {
    for(UInt i=1; i <= numEdgesL; ++i)
    {
      mapItem = grid2d->getEdges().getRowMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numEdgesG) + numVertexDofsG );
      
      dofMap(k) = mapItem;
      ++k;
    }
  }
  
  for(UInt lev=1; lev <= FETYPE::dofPerVolume; ++lev)  //Volume
  {
    for(UInt i=1; i <= numElementsL; ++i)
    {
      mapItem = grid2d->getElements().getRowMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numElementsG) + numVertexDofsG + numEdgesDofsG );
      
      dofMap(k) = mapItem;
      ++k;
    }
  }
  
  
  //Build list map_________________________________________________________________________________
  listMap.resize(sizeListL);
  
  UInt lid, gid;
  
  for(UInt i=1; i <= dofMap.size(); ++i)
  {   
    for(UInt K=1; K <= block; ++K)
    {
      mapItem = dofMap(i);
      
      lid = (block * (mapItem.getLid() - 1) + K) * (ORDER == dms2d_vectMajor)  + (mapItem.getLid() + (K-1) * numDofsL) * (ORDER == dms2d_componentMajor);
      gid = (block * (mapItem.getGid() - 1) + K) * (ORDER == dms2d_vectMajor)  + (mapItem.getGid() + (K-1) * numDofsG) * (ORDER == dms2d_componentMajor);
      
      mapItem.setLid(lid);
      mapItem.setGid(gid);
      
      listMap(lid) = mapItem;
    }
  }
  
  
  //Testing________________________________________________________________________________________
  pMapGlobalManip<PMAPTYPE> tester(commDev); 
  assert(tester.check(dofMap));
  assert(tester.check(listMap));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
bool
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
isActive(const DOFCARD & dofCard) const
{
  assert( (&dofCard) == (&dofCard));
  return(true);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
UInt
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofL - 1) + K)  * (ORDER == dms2d_vectMajor) +
          (dofL  + (K-1) * numDofsL) * (ORDER == dms2d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
UInt
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofG - 1) + K)   * (ORDER == dms2d_vectMajor) +
          (dofG  + (K-1) * numDofsG) * (ORDER == dms2d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
UInt
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
mapDofL(const DOFCARD & dofCard) const
{
  assert(startupOk);
  
  switch(dofCard.getGeoType())
  {
    case VERTEX :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVertex);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numVertices);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid2d->getNumElements());     
      return( grid2d->getElements().getCid_LL(dofCard.getLocalElId(), dofCard.getLocalId()) + (dofCard.getLevel() - 1) * numVerticesL);
    
    case EDGE   :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerEdge);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numEdges);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid2d->getNumElements());
      return( connectGrid2d->getElementToEdge(dofCard.getLocalElId(), dofCard.getLocalId()) + (dofCard.getLevel() - 1) * numEdgesL  + numVertexDofsL );
    
    case VOLUME :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVolume);
      assert(dofCard.getLocalId()   == 1);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid2d->getNumElements());
      return( dofCard.getLocalElId() + (dofCard.getLevel() - 1) * numElementsL + numVertexDofsL + numEdgesDofsL );
      
    default :
      assert(1==2);
  }
  
  return(0);             
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
UInt
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
mapDofG(const DOFCARD & dofCard) const
{
  assert(startupOk);
  return(dofMap(mapDofL(dofCard)).getGid());
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
UInt
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofL = mapDofL(dofCard);
  return(dofToListL(dofL,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
UInt
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofG = mapDofG(dofCard);
  return(dofToListG(dofG,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const typename dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::DOFMAP &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
getDofMap() const
{
  assert(startupOk);
  return(dofMap);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const typename dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::LISTMAP &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
getListMap() const
{
  assert(startupOk);
  return(listMap);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const UInt &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
getNumDofsL() const
{
  assert(startupOk);
  return(numDofsL);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const UInt &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
getNumDofsG() const
{
  assert(startupOk);
  return(numDofsG);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const UInt &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
getSizeListL() const
{
  assert(startupOk);
  return(sizeListL);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const UInt &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
getSizeListG() const
{
  assert(startupOk);
  return(sizeListG);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const typename dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::FECARDS &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
getFeCards() const
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
typename dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::FECARDS &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
getFeCards()
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const Teuchos::RCP<communicator> &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
getCommDev() const
{
  assert(commDevLoaded);
  return(commDev);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const Teuchos::RCP<typename dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::MESH2D> &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
getGrid2d() const
{
  assert(geometryLoaded);
  return(grid2d);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const Teuchos::RCP<typename dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::CONNECT2D> &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
getConnectGrid2d() const
{
  assert(geometryLoaded);
  return(connectGrid2d);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const Teuchos::RCP<typename dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::OPTIONS> &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
getOptions() const
{
  assert(optionsLoaded);
  return(options);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const bool &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_allMode>::
getStartupOk() const
{
  return(startupOk);
}

template<typename F, typename D, dms2d_order O>
ostream & operator<<(ostream & f, const dofMapStatic2d<F,D, O,dms2d_allMode> & V)
{
  f << "numVertices   : " << V.numVerticesL   << " " << V.numVerticesG << endl;
  f << "numEdges      : " << V.numEdgesL      << " " << V.numEdgesG << endl;
  f << "numElements   : " << V.numElementsL   << " " << V.numElementsG << endl << endl;
  
  f << "numVertexDofs : " << V.numVertexDofsL << " " << V.numVertexDofsG << endl;
  f << "numEdgesDofs  : " << V.numEdgesDofsL  << " " << V.numEdgesDofsG << endl;
  f << "numVolumeDofs : " << V.numVolumeDofsL << " " << V.numVolumeDofsG << endl << endl;
  
  f << "numTotalDofs  : " << V.numDofsL       << " " << V.numDofsG << endl;
  f << "sizeList      : " << V.sizeListL      << " " << V.sizeListG << endl;
  
  return(f);
}





//_________________________________________________________________________________________________
// GEO ID MODE
//-------------------------------------------------------------------------------------------------

/*! Support for the mapping of 2d static finite elements -> \c geoIdMode */
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
class dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>
{
  /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE   PMAPTYPE;
    typedef typename FETYPE::GEOSHAPE          GEOSHAPE;
    typedef typename FETYPE::DOFCARD           DOFCARD;
    typedef typename FETYPE::FECARD            FECARD;
    
    typedef mesh2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH2D;
    typedef connect2d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT2D;
    typedef geoMapInterface<GEOSHAPE>             GEOMAPINTERFACE;
    typedef pVect<FECARD,PMAPTYPE>                FECARDS;
    typedef dofMapStatic2d_options                OPTIONS;
    
    typedef pMap<PMAPTYPE> DOFMAP;
    typedef pMap<PMAPTYPE> LISTMAP;
    //@}
    
    /*! @name Internal structures */ //@{
  public:
    FETYPE          refFe;
    GEOMAPINTERFACE refShape;
    FECARDS         feCards;
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<communicator>  commDev;
    Teuchos::RCP<MESH2D>        grid2d;
    Teuchos::RCP<CONNECT2D>     connectGrid2d;
    Teuchos::RCP<OPTIONS>       options;
    //@}
    
    /*! @name Internal flags */ //@{
  public:
    bool commDevLoaded, geometryLoaded, optionsLoaded, startupOk;
    //@}
    
     /*! @name Static data */ //@{
  public:
    static const UInt block = traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    sVect<bool> vertexIsActive;
    sVect<bool> edgeIsActive;
    sVect<bool> elementIsActive;
    
    sVect<UInt> newVertexLid;
    sVect<UInt> newEdgeLid;
    sVect<UInt> newElementLid;
    
    pVect<UInt,PMAPTYPE>  globVertices;
    pVect<UInt,PMAPTYPE>  globEdges;
    pVect<UInt,PMAPTYPE>  globElements;
    
    UInt numVerticesL, numVerticesG;
    UInt numEdgesL,    numEdgesG;
    UInt numElementsL, numElementsG;
    
    UInt numVertexDofsL, numVertexDofsG;
    UInt numEdgesDofsL,  numEdgesDofsG;
    UInt numVolumeDofsL, numVolumeDofsG;
    
    UInt numDofsL, numDofsG;
    UInt sizeListL, sizeListG;
    //@}
    
    /*! @name Internal Maps */ //@{
  public:
    DOFMAP  dofMap;
    LISTMAP listMap;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    dofMapStatic2d();
    dofMapStatic2d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d);
    dofMapStatic2d(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d);
    dofMapStatic2d(const dofMapStatic2d & DofMap);
    dofMapStatic2d operator=(const dofMapStatic2d & DofMap);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d);
    void setGeometry(MESH2D & Grid2d, CONNECT2D & ConnedGrid2d);
    //@}
    
    /*! @name Options loading and startup */ //@{
  public:
    void setFeCardL(const UInt & lid, const FECARD & FeCards);
    void setFeCardG(const UInt & gid, const FECARD & FeCards);
    void setFeCards(const FECARDS & FECards);
    void setOptions(const Teuchos::RCP<OPTIONS> & Options);
    void setOptions(OPTIONS & Options);
    void startup();
    //@}
    
    /*! @name Get functions */ //@{
  public:
    bool isActive(const DOFCARD & dofCard) const;
    UInt dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const;
    UInt dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const;
    UInt mapDofL(const DOFCARD & dofCard) const;
    UInt mapDofG(const DOFCARD & dofCard) const;    
    UInt mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const;
    UInt mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const;
    const UInt & getNumDofsL() const;
    const UInt & getNumDofsG() const;
    const UInt & getSizeListL() const;
    const UInt & getSizeListG() const;
    //@}
    
    /*! @name Dump functions */ //@{
  public:
    const DOFMAP                     & getDofMap() const;
    const LISTMAP                    & getListMap() const;
    const FECARDS                    & getFeCards() const;
          FECARDS                    & getFeCards();
    const Teuchos::RCP<communicator> & getCommDev() const;
    const Teuchos::RCP<MESH2D>       & getGrid2d() const;
    const Teuchos::RCP<CONNECT2D>    & getConnectGrid2d() const;
    const Teuchos::RCP<OPTIONS>      & getOptions() const;
    const bool                       & getStartupOk() const;
    //@}
    
    /*! @name Extra dump functions */ //@{
  public:
    const sVect<bool> & getVertexIsActive() const;
    const sVect<bool> & getEdgeIsActive() const;
    const sVect<bool> & getElementIsActive() const;
    
    const sVect<UInt> & getNewVertexLid() const;
    const sVect<UInt> & getNewEdgeLid() const;
    const sVect<UInt> & getNewElementLid() const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    template<typename F, typename D, dms2d_order O>
    friend ostream & operator<<(ostream & f, const dofMapStatic2d<F,D, O,dms2d_allMode> & V);
    //@}
};



template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
dofMapStatic2d()
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  
  commDevLoaded  = false;
  geometryLoaded = false;
  optionsLoaded  = false;
  startupOk      = false;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
dofMapStatic2d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = CommDev;
  grid2d        = Grid2d;
  connectGrid2d = ConnedGrid2d;
  
  feCards.resize(grid2d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
dofMapStatic2d(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = Teuchos::rcpFromRef(CommDev);
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnedGrid2d));
  
  feCards.resize(grid2d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
dofMapStatic2d(const dofMapStatic2d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid2d        = DofMap.grid2d;
  connectGrid2d = DofMap.connectGrid2d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;
     
  vertexIsActive  = DofMap.vertexIsActive;
  edgeIsActive    = DofMap.edgeIsActive;
  elementIsActive = DofMap.elementIsActive;
   
  newVertexLid  = DofMap.newVertexLid;
  newEdgeLid    = DofMap.newEdgeLid; 
  newElementLid = DofMap.newElementLid;
    
  globVertices = DofMap.globVertices;
  globEdges    = DofMap.globEdges;
  globElements = DofMap.globElements;
    
  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numEdgesL    = DofMap.numEdgesL;
  numEdgesG    = DofMap.numEdgesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numEdgesDofsL  = DofMap.numEdgesDofsL;
  numEdgesDofsG  = DofMap.numEdgesDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  sizeListL = DofMap.sizeListL;
  sizeListG = DofMap.sizeListG;
    
  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
operator=(const dofMapStatic2d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid2d        = DofMap.grid2d;
  connectGrid2d = DofMap.connectGrid2d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;
     
  vertexIsActive  = DofMap.vertexIsActive;
  edgeIsActive    = DofMap.edgeIsActive;
  elementIsActive = DofMap.elementIsActive;
   
  newVertexLid  = DofMap.newVertexLid;
  newEdgeLid    = DofMap.newEdgeLid; 
  newElementLid = DofMap.newElementLid;
    
  globVertices = DofMap.globVertices;
  globEdges    = DofMap.globEdges;
  globElements = DofMap.globElements;
    
  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numEdgesL    = DofMap.numEdgesL;
  numEdgesG    = DofMap.numEdgesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numEdgesDofsL  = DofMap.numEdgesDofsL;
  numEdgesDofsG  = DofMap.numEdgesDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  sizeListL = DofMap.sizeListL;
  sizeListG = DofMap.sizeListG;
    
  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
  
  return(*this);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded  = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded  = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d)
{
  geometryLoaded = true;
  grid2d         = Grid2d;
  connectGrid2d  = ConnedGrid2d;
  
  feCards.resize(grid2d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
setGeometry(MESH2D & Grid2d, CONNECT2D & ConnedGrid2d)
{
  geometryLoaded = true; 
  grid2d         = Teuchos::rcp(new MESH2D(Grid2d));
  connectGrid2d  = Teuchos::rcp(new CONNECT2D(ConnedGrid2d));
  
  feCards.resize(grid2d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
setFeCardL(const UInt & lid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(lid <= feCards.sizeL());
  feCards.getL(lid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
setFeCardG(const UInt & gid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(feCards.isG(gid));
  feCards.getG(gid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
setFeCards(const FECARDS & FECards)
{
  assert(geometryLoaded);
  assert(feCards.size() == FECards.size());
  feCards = FECards;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
setOptions(const Teuchos::RCP<OPTIONS> & Options)
{
  optionsLoaded = true;
  options = Options;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
setOptions(OPTIONS & Options)
{
  optionsLoaded = true;
  options = Teuchos::rcp(new OPTIONS(Options));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
void
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
startup()
{
  assert(commDevLoaded);
  assert(optionsLoaded);
  assert(optionsLoaded);
  
  assert(grid2d->getElements().colIsLocal());
  assert(grid2d->getEdges().colIsLocal());
  
  mesh2dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> manipulator(commDev);
  assert(manipulator.check(grid2d));
  
  startupOk = true;
  
  //Alloc__________________________________________________________________________________________
  vertexIsActive.resize(  grid2d->getNumVertices() );
  edgeIsActive.resize(    grid2d->getNumEdges()    );
  elementIsActive.resize( grid2d->getNumElements() );
   
  newVertexLid.resize(  grid2d->getNumVertices() );
  newEdgeLid.resize(    grid2d->getNumEdges()    );
  newElementLid.resize( grid2d->getNumElements() );
  
  for(UInt i=1; i <= grid2d->getNumVertices(); ++i)
  { vertexIsActive(i) = false; }
  
  for(UInt i=1; i <= grid2d->getNumEdges(); ++i)
  { edgeIsActive(i) = false; }
  
  for(UInt i=1; i <= grid2d->getNumElements(); ++i)
  { elementIsActive(i) = false; }
  
  //Identify the active items______________________________________________________________________
  UInt geoId;
  
  for(UInt i=1; i <= grid2d->getNumElements(); ++i)
  {
    geoId = grid2d->getElementL(i).getGeoId();
    
    if(options->isGeoId(geoId))
    {
      for(UInt j=1; j <= GEOSHAPE::numVertices; ++j)
      {
	vertexIsActive(grid2d->getElementL(i).getCid(j)) = true;
      }
      
      for(UInt j=1; j <= GEOSHAPE::numEdges; ++j)
      {
	edgeIsActive(connectGrid2d->getElementToEdge(i,j)) = true;
      }
      
      elementIsActive(i) = true;
    }
  }
  
  //Assign new local ids___________________________________________________________________________
  numVerticesL = 0;
  numEdgesL    = 0;
  numElementsL = 0;
  
  globVertices.clear();
  globEdges.clear();
  globElements.clear();
  
  PMAPTYPE nodeMapItem;
  PMAPTYPE elMapItem;
  
  for(UInt i=1; i <= grid2d->getNumVertices(); ++i)  //Vertices
  {
    if(vertexIsActive(i))
    {
      ++numVerticesL;
      newVertexLid(i) = numVerticesL;
      
      nodeMapItem = grid2d->getNodes().getMapL(i);
      nodeMapItem.setLid(numVerticesL);
      
      globVertices.push_back(nodeMapItem, nodeMapItem.getGid());
    }
  }
  
  for(UInt i=1; i <= grid2d->getNumEdges(); ++i)     //Edges
  {
    if(edgeIsActive(i))
    {
      ++numEdgesL;
      newEdgeLid(i) = numEdgesL;
      
      elMapItem = grid2d->getEdges().getRowMapL(i);
      elMapItem.setLid(numEdgesL);
      
      globEdges.push_back(elMapItem, elMapItem.getGid());
    }
  }
  
  for(UInt i=1; i <= grid2d->getNumElements(); ++i)  //Elements
  {
    if(elementIsActive(i))
    {
      ++numElementsL;
      newElementLid(i) = numElementsL;
      
      elMapItem = grid2d->getElements().getRowMapL(i);
      elMapItem.setLid(numElementsL);
      
      globElements.push_back(elMapItem, elMapItem.getGid());
    }
  }
  
  //Assign new global ids__________________________________________________________________________
  pVectGlobalManip<UInt,PMAPTYPE> nodeManip(commDev);
  pVectGlobalManip<UInt,PMAPTYPE>   elManip(commDev);
   
  nodeManip.buildGlobalNumbering(globVertices);
  elManip.buildGlobalNumbering(globEdges);
  elManip.buildGlobalNumbering(globElements);
  
  numVerticesG = nodeManip.sizeG(globVertices);
  numEdgesG    = elManip.sizeG(globEdges);
  numElementsG = elManip.sizeG(globElements);
   
  //Count the dofs_________________________________________________________________________________
  numVertexDofsL = FETYPE::dofPerVertex * numVerticesL;
  numEdgesDofsL  = FETYPE::dofPerEdge   * numEdgesL;
  numVolumeDofsL = FETYPE::dofPerVolume * numElementsL;
  
  numVertexDofsG = FETYPE::dofPerVertex * numVerticesG;
  numEdgesDofsG  = FETYPE::dofPerEdge   * numEdgesG;
  numVolumeDofsG = FETYPE::dofPerVolume * numElementsG;
  
  //Total number of dofs and lists_________________________________________________________________
  numDofsL = numVertexDofsL + numEdgesDofsL + numVolumeDofsL;
  numDofsG = numVertexDofsG + numEdgesDofsG + numVolumeDofsG;
  
  sizeListL = numDofsL * traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
  sizeListG = numDofsG * traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
  
  
  //Build dof map__________________________________________________________________________________
  PMAPTYPE mapItem;
  UInt k=1;
  
  dofMap.resize(numDofsL);
  
  for(UInt lev=1; lev <= FETYPE::dofPerVertex; ++lev)  //Vertices
  {
    for(UInt i=1; i <= numVerticesL; ++i)
    {
      mapItem = globVertices.getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numVerticesG) );
      
      dofMap(k) = mapItem;
      ++k;
    } 
  }
  
  for(UInt lev=1; lev <= FETYPE::dofPerEdge; ++lev)  //Edges
  {
    for(UInt i=1; i <= numEdgesL; ++i)
    {
      mapItem = globEdges.getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numEdgesG) + numVertexDofsG );
      
      dofMap(k) = mapItem;
      ++k;
    }
  }
  
  for(UInt lev=1; lev <= FETYPE::dofPerVolume; ++lev)  //Volume
  {
    for(UInt i=1; i <= numElementsL; ++i)
    {
      mapItem = globElements.getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numElementsG) + numVertexDofsG + numEdgesDofsG );
      
      dofMap(k) = mapItem;
      ++k;
    }
  }
  
  
  //Build list map_________________________________________________________________________________
  listMap.resize(sizeListL);
  
  UInt lid, gid;
  
  for(UInt i=1; i <= dofMap.size(); ++i)
  {   
    for(UInt K=1; K <= block; ++K)
    {
      mapItem = dofMap(i);
      
      lid = (block * (mapItem.getLid() - 1) + K) * (ORDER == dms2d_vectMajor)  + (mapItem.getLid() + (K-1) * numDofsL) * (ORDER == dms2d_componentMajor);
      gid = (block * (mapItem.getGid() - 1) + K) * (ORDER == dms2d_vectMajor)  + (mapItem.getGid() + (K-1) * numDofsG) * (ORDER == dms2d_componentMajor);
      
      mapItem.setLid(lid);
      mapItem.setGid(gid);
      
      listMap(lid) = mapItem;
    }
  }
  
  
  //Testing________________________________________________________________________________________
  pMapGlobalManip<PMAPTYPE> tester(commDev); 
  assert(tester.check(dofMap));
  assert(tester.check(listMap));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
bool
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
isActive(const DOFCARD & dofCard) const
{
  assert(dofCard.getLocalElId() <= elementIsActive.size());
  return(elementIsActive(dofCard.getLocalElId()));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
UInt
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofL - 1) + K)  * (ORDER == dms2d_vectMajor) +
          (dofL  + (K-1) * numDofsL) * (ORDER == dms2d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
UInt
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofG - 1) + K)   * (ORDER == dms2d_vectMajor) +
          (dofG  + (K-1) * numDofsG) * (ORDER == dms2d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
UInt
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
mapDofL(const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(elementIsActive(dofCard.getLocalElId()));
  
  UInt oldId;
  
  switch(dofCard.getGeoType())
  {
    case VERTEX :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVertex);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numVertices);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid2d->getNumElements());
      
      oldId = grid2d->getElements().getCid_LL(dofCard.getLocalElId(), dofCard.getLocalId());
      return( newVertexLid(oldId) + (dofCard.getLevel() - 1) * numVerticesL);
    
    case EDGE   :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerEdge);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numEdges);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid2d->getNumElements());
      
      oldId = connectGrid2d->getElementToEdge(dofCard.getLocalElId(), dofCard.getLocalId());
      return( newEdgeLid(oldId) + (dofCard.getLevel() - 1) * numEdgesL + numVertexDofsL );
    
    case VOLUME :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVolume);
      assert(dofCard.getLocalId()   == 1);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid2d->getNumElements());
      
      oldId = dofCard.getLocalElId();
      return( newElementLid(oldId) + (dofCard.getLevel() - 1) * numElementsL + numVertexDofsL + numEdgesDofsL );
      
    default :
      assert(1==2);
   }
  
  return(0); 
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
UInt
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
mapDofG(const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(elementIsActive(dofCard.getLocalElId()));
  
  return(dofMap(mapDofL(dofCard)).getGid());
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
UInt
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofL = mapDofL(dofCard);
  return(dofToListL(dofL,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
UInt
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofG = mapDofG(dofCard);
  return(dofToListG(dofG,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const typename dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::DOFMAP &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getDofMap() const
{
  return(dofMap);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const typename dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::LISTMAP &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getListMap() const
{
  return(listMap);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const UInt &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getNumDofsL() const
{
  return(numDofsL);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const UInt &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getNumDofsG() const
{
  return(numDofsG);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const UInt &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getSizeListL() const
{
  return(sizeListL);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const UInt &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getSizeListG() const
{
  return(sizeListG);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const typename dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::FECARDS &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getFeCards() const
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
typename dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::FECARDS &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getFeCards()
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const Teuchos::RCP<communicator> &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getCommDev() const
{
  assert(commDevLoaded);
  return(commDev);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const Teuchos::RCP<typename dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::MESH2D> &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getGrid2d() const
{
  assert(geometryLoaded);
  return(grid2d);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const Teuchos::RCP<typename dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::CONNECT2D> &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getConnectGrid2d() const
{
  assert(geometryLoaded);
  return(connectGrid2d);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const Teuchos::RCP<typename dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::OPTIONS> &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getOptions() const
{
  assert(optionsLoaded);
  return(options);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const sVect<bool> &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getVertexIsActive() const
{
  assert(startupOk);
  return(vertexIsActive);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const sVect<bool> &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getEdgeIsActive() const
{
  assert(startupOk);
  return(edgeIsActive);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const sVect<bool> &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getElementIsActive() const
{
  assert(startupOk);
  return(elementIsActive);
}
    
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const sVect<UInt> &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getNewVertexLid() const
{
  assert(startupOk);
  return(newVertexLid);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const sVect<UInt> &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getNewEdgeLid() const
{
  assert(startupOk);
  return(newEdgeLid);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const sVect<UInt> &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getNewElementLid() const
{
  assert(startupOk);
  return(newElementLid);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER>
const bool &
dofMapStatic2d<FETYPE,DOFTYPE,ORDER,dms2d_geoIdMode>::
getStartupOk() const
{
  return(startupOk);
}

template<typename F, typename D, dms2d_order O>
ostream & operator<<(ostream & f, const dofMapStatic2d<F,D, O,dms2d_geoIdMode> & V)
{
  f << "numVertices   : " << V.numVerticesL   << " " << V.numVerticesG << endl;
  f << "numEdges      : " << V.numEdgesL      << " " << V.numEdgesG << endl;
  f << "numElements   : " << V.numElementsL   << " " << V.numElementsG << endl << endl;
  
  f << "numVertexDofs : " << V.numVertexDofsL << " " << V.numVertexDofsG << endl;
  f << "numEdgesDofs  : " << V.numEdgesDofsL  << " " << V.numEdgesDofsG << endl;
  f << "numVolumeDofs : " << V.numVolumeDofsL << " " << V.numVolumeDofsG << endl << endl;
  
  f << "numTotalDofs  : " << V.numDofsL       << " " << V.numDofsG << endl;
  f << "sizeList      : " << V.sizeListL      << " " << V.sizeListG << endl;
  
  return(f);
}




#endif
