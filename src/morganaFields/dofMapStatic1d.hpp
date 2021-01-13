/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef DOFMAPSTATIC1D_HPP
#define DOFMAPSTATIC1D_HPP

enum dms1d_order {dms1d_vectMajor, dms1d_componentMajor};
enum dms1d_mode  {dms1d_allMode, dms1d_geoIdMode, dms1d_disjointMode};

#include "typesInterface.hpp"
#include "traitsBasic.h"

#include "pVectManip.hpp"

#include "morganaGeometry.hpp"
#include "geoMapInterface.hpp"
#include "connect1d.hpp"
#include "mesh1dGlobalManip.hpp"

#include "dofMapStatic1d_options.h"
#include "feStaticDofCard1d.h"


//! Forward declaration
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER = dms1d_vectMajor, dms1d_mode MODE = dms1d_allMode> class dofMapStatic1d;



//_________________________________________________________________________________________________
// ALL MODE
//-------------------------------------------------------------------------------------------------

/*! Support for the mapping of 1d static finite elements -> \c allMode */
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
class dofMapStatic1d<FETYPE,DOFTYPE, ORDER,dms1d_allMode>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE   PMAPTYPE;
    typedef typename FETYPE::GEOSHAPE          GEOSHAPE;
    typedef typename FETYPE::DOFCARD           DOFCARD;
    typedef typename FETYPE::FECARD            FECARD;
    
    typedef mesh1d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH1D;
    typedef connect1d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT1D;
    typedef geoMapInterface<GEOSHAPE>             GEOMAPINTERFACE;
    typedef pVect<FECARD,PMAPTYPE>                FECARDS;
    typedef dofMapStatic1d_options                OPTIONS;
    
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
    Teuchos::RCP<MESH1D>        grid1d;
    Teuchos::RCP<CONNECT1D>     connectGrid1d;
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
    UInt numElementsL, numElementsG;
    
    UInt numVertexDofsL, numVertexDofsG;
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
    dofMapStatic1d();
    dofMapStatic1d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d);
    dofMapStatic1d(communicator & CommDev, MESH1D & Grid1d, CONNECT1D & ConnedGrid1d);
    dofMapStatic1d(const dofMapStatic1d & DofMap);
    dofMapStatic1d operator=(const dofMapStatic1d & DofMap);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d);
    void setGeometry(MESH1D & Grid1d, CONNECT1D & ConnedGrid1d);
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
    const Teuchos::RCP<MESH1D>       & getGrid1d() const;
    const Teuchos::RCP<CONNECT1D>    & getConnectGrid1d() const;
    const Teuchos::RCP<OPTIONS>      & getOptions() const;
    const bool                       & getStartupOk() const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    template<typename F, typename D, dms1d_order O>
    friend ostream & operator<<(ostream & f, const dofMapStatic1d<F,D, O,dms1d_allMode> & V);
    //@}
};


template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
dofMapStatic1d()
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  commDevLoaded  = false;
  geometryLoaded = false;
  optionsLoaded  = false;
  startupOk      = false;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
dofMapStatic1d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = CommDev;
  grid1d        = Grid1d;
  connectGrid1d = ConnedGrid1d;
  
  feCards.resize(grid1d->getNumElements());
  feCards.setMap(grid1d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
dofMapStatic1d(communicator & CommDev, MESH1D & Grid1d, CONNECT1D & ConnedGrid1d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = Teuchos::rcpFromRef(CommDev);
  grid1d        = Teuchos::rcp(new MESH1D(Grid1d));
  connectGrid1d = Teuchos::rcp(new CONNECT1D(ConnedGrid1d));
  
  feCards.resize(grid1d->getNumElements());
  feCards.setMap(grid1d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
dofMapStatic1d(const dofMapStatic1d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid1d        = DofMap.grid1d;
  connectGrid1d = DofMap.connectGrid1d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;

  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  sizeListL = DofMap.sizeListL;
  sizeListG = DofMap.sizeListG;
    
  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
operator=(const dofMapStatic1d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid1d        = DofMap.grid1d;
  connectGrid1d = DofMap.connectGrid1d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;

  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded  = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded  = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
setGeometry(const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d)
{
  geometryLoaded = true;
  grid1d         = Grid1d;
  connectGrid1d  = ConnedGrid1d;
  
  feCards.resize(grid1d->getNumElements());
  feCards.setMap(grid1d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
setGeometry(MESH1D & Grid1d, CONNECT1D & ConnedGrid1d)
{
  geometryLoaded = true; 
  grid1d        = Teuchos::rcp(new MESH1D(Grid1d));
  connectGrid1d = Teuchos::rcp(new CONNECT1D(ConnedGrid1d));
  
  feCards.resize(grid1d->getNumElements());
  feCards.setMap(grid1d->getElements().getRowMap());
  feCards.updateFinder();
}



//_________________________________________________________________________________________________
// OPTIONS LOADING AND STARTUP
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
setFeCardL(const UInt & lid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(lid <= feCards.sizeL());
  feCards.getL(lid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
setFeCardG(const UInt & gid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(feCards.isG(gid));
  feCards.getG(gid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
setFeCards(const FECARDS & FECards)
{
  assert(geometryLoaded);
  assert(feCards.size() == FECards.size());
  feCards = FECards;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
setOptions(const Teuchos::RCP<OPTIONS> & Options)
{
  optionsLoaded = true;
  options = Options;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
setOptions(OPTIONS & Options)
{
  optionsLoaded = true;
  options = Teuchos::rcp(new OPTIONS(Options));
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
startup()
{
  //Assert and flags_______________________________________________________________________________
  assert(commDevLoaded);
  assert(geometryLoaded);
  assert(optionsLoaded);
  
  assert(grid1d->getElements().colIsLocal());
  
  mesh1dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> manipulator(commDev);
  assert(manipulator.check(grid1d));
  
  startupOk = true;
  
  
  //Count the geotypes_____________________________________________________________________________ 
  numVerticesL = grid1d->getNumVertices();
  numElementsL = grid1d->getNumElements();
  
  numVerticesG = manipulator.getNumGlobalVertices(grid1d);
  numElementsG = manipulator.getNumGlobalElements(grid1d);
  
  //Count the dofs
  numVertexDofsL = FETYPE::dofPerVertex * numVerticesL;
  numVolumeDofsL = FETYPE::dofPerVolume * numElementsL;
  
  numVertexDofsG = FETYPE::dofPerVertex * numVerticesG;
  numVolumeDofsG = FETYPE::dofPerVolume * numElementsG;
  
  
  //Total number of dofs and lists_________________________________________________________________
  numDofsL = numVertexDofsL + numVolumeDofsL;
  numDofsG = numVertexDofsG + numVolumeDofsG;
  
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
      mapItem = grid1d->getNodes().getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numVerticesG) );
      
      dofMap(k) = mapItem;
      ++k;
    } 
  }
  
  for(UInt lev=1; lev <= FETYPE::dofPerVolume; ++lev)  //Volume
  {
    for(UInt i=1; i <= numElementsL; ++i)
    {
      mapItem = grid1d->getElements().getRowMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numElementsG) + numVertexDofsG);
      
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
      
      lid = (block * (mapItem.getLid() - 1) + K) * (ORDER == dms1d_vectMajor)  + (mapItem.getLid() + (K-1) * numDofsL) * (ORDER == dms1d_componentMajor);
      gid = (block * (mapItem.getGid() - 1) + K) * (ORDER == dms1d_vectMajor)  + (mapItem.getGid() + (K-1) * numDofsG) * (ORDER == dms1d_componentMajor);
      
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
bool
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
isActive(const DOFCARD & dofCard) const
{
  return(true);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
UInt
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofL - 1) + K)  * (ORDER == dms1d_vectMajor) +
          (dofL  + (K-1) * numDofsL) * (ORDER == dms1d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
UInt
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofG - 1) + K)   * (ORDER == dms1d_vectMajor) +
          (dofG  + (K-1) * numDofsG) * (ORDER == dms1d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
UInt
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
mapDofL(const DOFCARD & dofCard) const
{
  assert(startupOk);
  
  switch(dofCard.getGeoType())
  {
    case VERTEX :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVertex);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numVertices);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid1d->getNumElements());     
      return( grid1d->getElements().getCid_LL(dofCard.getLocalElId(), dofCard.getLocalId()) + (dofCard.getLevel() - 1) * numVerticesL);
    
    case VOLUME :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVolume);
      assert(dofCard.getLocalId()   == 1);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid1d->getNumElements());
      return( dofCard.getLocalElId() + (dofCard.getLevel() - 1) * numElementsL + numVertexDofsL );
      
    default :
      assert(1==2);
  }
  
  return(0);             
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
UInt
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
mapDofG(const DOFCARD & dofCard) const
{
  assert(startupOk);
  return(dofMap(mapDofL(dofCard)).getGid());
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
UInt
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofL = mapDofL(dofCard);
  return(dofToListL(dofL,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
UInt
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofG = mapDofG(dofCard);
  return(dofToListG(dofG,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const typename dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::DOFMAP &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
getDofMap() const
{
  assert(startupOk);
  return(dofMap);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const typename dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::LISTMAP &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
getListMap() const
{
  assert(startupOk);
  return(listMap);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const UInt &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
getNumDofsL() const
{
  assert(startupOk);
  return(numDofsL);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const UInt &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
getNumDofsG() const
{
  assert(startupOk);
  return(numDofsG);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const UInt &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
getSizeListL() const
{
  assert(startupOk);
  return(sizeListL);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const UInt &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
getSizeListG() const
{
  assert(startupOk);
  return(sizeListG);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const typename dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::FECARDS &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
getFeCards() const
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
typename dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::FECARDS &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
getFeCards()
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const Teuchos::RCP<communicator> &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
getCommDev() const
{
  assert(commDevLoaded);
  return(commDev);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const Teuchos::RCP<typename dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::MESH1D> &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
getGrid1d() const
{
  assert(geometryLoaded);
  return(grid1d);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const Teuchos::RCP<typename dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::CONNECT1D> &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
getConnectGrid1d() const
{
  assert(geometryLoaded);
  return(connectGrid1d);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const Teuchos::RCP<typename dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::OPTIONS> &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
getOptions() const
{
  assert(optionsLoaded);
  return(options);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const bool &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_allMode>::
getStartupOk() const
{
  return(startupOk);
}

template<typename F, typename D, dms1d_order O>
ostream & operator<<(ostream & f, const dofMapStatic1d<F,D, O,dms1d_allMode> & V)
{
  f << "numVertices   : " << V.numVerticesL   << " " << V.numVerticesG << endl;
  f << "numElements   : " << V.numElementsL   << " " << V.numElementsG << endl << endl;
  
  f << "numVertexDofs : " << V.numVertexDofsL << " " << V.numVertexDofsG << endl;
  f << "numVolumeDofs : " << V.numVolumeDofsL << " " << V.numVolumeDofsG << endl << endl;
  
  f << "numTotalDofs  : " << V.numDofsL       << " " << V.numDofsG << endl;
  f << "sizeList      : " << V.sizeListL      << " " << V.sizeListG << endl;
  
  return(f);
}





//_________________________________________________________________________________________________
// GEO ID MODE
//-------------------------------------------------------------------------------------------------

/*! Support for the mapping of 1d static finite elements -> \c geoIdMode */
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
class dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>
{
  /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE   PMAPTYPE;
    typedef typename FETYPE::GEOSHAPE          GEOSHAPE;
    typedef typename FETYPE::DOFCARD           DOFCARD;
    typedef typename FETYPE::FECARD            FECARD;
    
    typedef mesh1d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH1D;
    typedef connect1d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT1D;
    typedef geoMapInterface<GEOSHAPE>             GEOMAPINTERFACE;
    typedef pVect<FECARD,PMAPTYPE>                FECARDS;
    typedef dofMapStatic1d_options                OPTIONS;
    
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
    Teuchos::RCP<MESH1D>        grid1d;
    Teuchos::RCP<CONNECT1D>     connectGrid1d;
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
    sVect<bool> elementIsActive;
    
    sVect<UInt> newVertexLid;
    sVect<UInt> newElementLid;
    
    pVect<UInt,PMAPTYPE>  globVertices;
    pVect<UInt,PMAPTYPE>  globElements;
    
    UInt numVerticesL, numVerticesG;
    UInt numElementsL, numElementsG;
    
    UInt numVertexDofsL, numVertexDofsG;
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
    dofMapStatic1d();
    dofMapStatic1d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d);
    dofMapStatic1d(communicator & CommDev, MESH1D & Grid1d, CONNECT1D & ConnedGrid1d);
    dofMapStatic1d(const dofMapStatic1d & DofMap);
    dofMapStatic1d operator=(const dofMapStatic1d & DofMap);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d);
    void setGeometry(MESH1D & Grid1d, CONNECT1D & ConnedGrid1d);
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
    const Teuchos::RCP<MESH1D>       & getGrid1d() const;
    const Teuchos::RCP<CONNECT1D>    & getConnectGrid1d() const;
    const Teuchos::RCP<OPTIONS>      & getOptions() const;
    const bool                       & getStartupOk() const;
    //@}
    
    /*! @name Extra dump functions */ //@{
  public:
    const sVect<bool> & getVertexIsActive() const;
    const sVect<bool> & getElementIsActive() const;
    
    const sVect<UInt> & getNewVertexLid() const;
    const sVect<UInt> & getNewElementLid() const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    template<typename F, typename D, dms1d_order O>
    friend ostream & operator<<(ostream & f, const dofMapStatic1d<F,D, O,dms1d_allMode> & V);
    //@}
};



template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
dofMapStatic1d()
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  commDevLoaded  = false;
  geometryLoaded = false;
  optionsLoaded  = false;
  startupOk      = false;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
dofMapStatic1d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = CommDev;
  grid1d        = Grid1d;
  connectGrid1d = ConnedGrid1d;
  
  feCards.resize(grid1d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
dofMapStatic1d(communicator & CommDev, MESH1D & Grid1d, CONNECT1D & ConnedGrid1d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = Teuchos::rcpFromRef(CommDev);
  grid1d        = Teuchos::rcp(new MESH1D(Grid1d));
  connectGrid1d = Teuchos::rcp(new CONNECT1D(ConnedGrid1d));
  
  feCards.resize(grid1d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
dofMapStatic1d(const dofMapStatic1d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid1d        = DofMap.grid1d;
  connectGrid1d = DofMap.connectGrid1d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;
     
  vertexIsActive  = DofMap.vertexIsActive;
  elementIsActive = DofMap.elementIsActive;
   
  newVertexLid  = DofMap.newVertexLid;
  newElementLid = DofMap.newElementLid;
    
  globVertices = DofMap.globVertices;
  globElements = DofMap.globElements;
    
  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  sizeListL = DofMap.sizeListL;
  sizeListG = DofMap.sizeListG;
    
  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
operator=(const dofMapStatic1d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid1d        = DofMap.grid1d;
  connectGrid1d = DofMap.connectGrid1d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;
     
  vertexIsActive  = DofMap.vertexIsActive;
  elementIsActive = DofMap.elementIsActive;
   
  newVertexLid  = DofMap.newVertexLid;
  newElementLid = DofMap.newElementLid;
    
  globVertices = DofMap.globVertices;
  globElements = DofMap.globElements;
    
  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded  = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded  = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
setGeometry(const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d)
{
  geometryLoaded = true;
  grid1d         = Grid1d;
  connectGrid1d  = ConnedGrid1d;
  
  feCards.resize(grid1d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
setGeometry(MESH1D & Grid1d, CONNECT1D & ConnedGrid1d)
{
  geometryLoaded = true; 
  grid1d        = Teuchos::rcp(new MESH1D(Grid1d));
  connectGrid1d = Teuchos::rcp(new CONNECT1D(ConnedGrid1d));
  
  feCards.resize(grid1d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
setFeCardL(const UInt & lid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(lid <= feCards.sizeL());
  feCards.getL(lid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
setFeCardG(const UInt & gid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(feCards.isG(gid));
  feCards.getG(gid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
setFeCards(const FECARDS & FECards)
{
  assert(geometryLoaded);
  assert(feCards.size() == FECards.size());
  feCards = FECards;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
setOptions(const Teuchos::RCP<OPTIONS> & Options)
{
  optionsLoaded = true;
  options = Options;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
setOptions(OPTIONS & Options)
{
  optionsLoaded = true;
  options = Teuchos::rcp(new OPTIONS(Options));
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
void
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
startup()
{
  assert(commDevLoaded);
  assert(optionsLoaded);
  assert(optionsLoaded);
  
  assert(grid1d->getElements().colIsLocal());
  
  mesh1dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> manipulator(commDev);
  assert(manipulator.check(grid1d));
  
  startupOk = true;
  
  //Alloc__________________________________________________________________________________________
  vertexIsActive.resize(  grid1d->getNumVertices() );
  elementIsActive.resize( grid1d->getNumElements() );
   
  newVertexLid.resize(  grid1d->getNumVertices() );
  newElementLid.resize( grid1d->getNumElements() );
  
  for(UInt i=1; i <= grid1d->getNumVertices(); ++i)
  { vertexIsActive(i) = false; }
  
  for(UInt i=1; i <= grid1d->getNumElements(); ++i)
  { elementIsActive(i) = false; }
  
  //Identify the active items______________________________________________________________________
  UInt geoId;
  
  for(UInt i=1; i <= grid1d->getNumElements(); ++i)
  {
    geoId = grid1d->getElementL(i).getGeoId();
    
    if(options->isGeoId(geoId))
    {
      for(UInt j=1; j <= GEOSHAPE::numVertices; ++j)
      { vertexIsActive(grid1d->getElementL(i).getCid(j)) = true; }
      
      elementIsActive(i) = true;
    }
  }
  
  //Assign new local ids___________________________________________________________________________
  numVerticesL = 0;
  numElementsL = 0;
  
  globVertices.clear();
  globElements.clear();
  
  PMAPTYPE nodeMapItem;
  PMAPTYPE elMapItem;
  
  for(UInt i=1; i <= grid1d->getNumVertices(); ++i)  //Vertices
  {
    if(vertexIsActive(i))
    {
      ++numVerticesL;
      newVertexLid(i) = numVerticesL;
      
      nodeMapItem = grid1d->getNodes().getMapL(i);
      nodeMapItem.setLid(numVerticesL);
      
      globVertices.push_back(nodeMapItem, nodeMapItem.getGid());
    }
  }
  
  for(UInt i=1; i <= grid1d->getNumElements(); ++i)  //Elements
  {
    if(elementIsActive(i))
    {
      ++numElementsL;
      newElementLid(i) = numElementsL;
      
      elMapItem = grid1d->getElements().getRowMapL(i);
      elMapItem.setLid(numElementsL);
      
      globElements.push_back(elMapItem, elMapItem.getGid());
    }
  }
  
  //Assign new global ids__________________________________________________________________________
  pVectGlobalManip<UInt,PMAPTYPE> nodeManip(commDev);
  pVectGlobalManip<UInt,PMAPTYPE>   elManip(commDev);
   
  nodeManip.buildGlobalNumbering(globVertices);
  elManip.buildGlobalNumbering(globElements);
  
  numVerticesG = nodeManip.sizeG(globVertices);
  numElementsG = elManip.sizeG(globElements);
   
  //Count the dofs_________________________________________________________________________________
  numVertexDofsL = FETYPE::dofPerVertex * numVerticesL;
  numVolumeDofsL = FETYPE::dofPerVolume * numElementsL;
  
  numVertexDofsG = FETYPE::dofPerVertex * numVerticesG;
  numVolumeDofsG = FETYPE::dofPerVolume * numElementsG;
  
  //Total number of dofs and lists_________________________________________________________________
  numDofsL = numVertexDofsL + numVolumeDofsL;
  numDofsG = numVertexDofsG + numVolumeDofsG;
  
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
  
  for(UInt lev=1; lev <= FETYPE::dofPerVolume; ++lev)  //Volume
  {
    for(UInt i=1; i <= numElementsL; ++i)
    {
      mapItem = globElements.getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numElementsG) + numVertexDofsG );
      
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
      
      lid = (block * (mapItem.getLid() - 1) + K) * (ORDER == dms1d_vectMajor)  + (mapItem.getLid() + (K-1) * numDofsL) * (ORDER == dms1d_componentMajor);
      gid = (block * (mapItem.getGid() - 1) + K) * (ORDER == dms1d_vectMajor)  + (mapItem.getGid() + (K-1) * numDofsG) * (ORDER == dms1d_componentMajor);
      
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
bool
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
isActive(const DOFCARD & dofCard) const
{
  assert(dofCard.getLocalElId() <= elementIsActive.size());
  return(elementIsActive(dofCard.getLocalElId()));
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
UInt
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofL - 1) + K)  * (ORDER == dms1d_vectMajor) +
          (dofL  + (K-1) * numDofsL) * (ORDER == dms1d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
UInt
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofG - 1) + K)   * (ORDER == dms1d_vectMajor) +
          (dofG  + (K-1) * numDofsG) * (ORDER == dms1d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
UInt
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
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
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid1d->getNumElements());
      
      oldId = grid1d->getElements().getCid_LL(dofCard.getLocalElId(), dofCard.getLocalId());
      return( newVertexLid(oldId) + (dofCard.getLevel() - 1) * numVerticesL);
    
    case VOLUME :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVolume);
      assert(dofCard.getLocalId()   == 1);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid1d->getNumElements());
      
      oldId = dofCard.getLocalElId();
      return( newElementLid(oldId) + (dofCard.getLevel() - 1) * numElementsL + numVertexDofsL );
      
    default :
      assert(1==2);
   }
  
  return(0); 
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
UInt
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
mapDofG(const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(elementIsActive(dofCard.getLocalElId()));
  
  return(dofMap(mapDofL(dofCard)).getGid());
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
UInt
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofL = mapDofL(dofCard);
  return(dofToListL(dofL,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
UInt
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofG = mapDofG(dofCard);
  return(dofToListG(dofG,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const typename dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::DOFMAP &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getDofMap() const
{
  return(dofMap);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const typename dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::LISTMAP &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getListMap() const
{
  return(listMap);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const UInt &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getNumDofsL() const
{
  return(numDofsL);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const UInt &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getNumDofsG() const
{
  return(numDofsG);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const UInt &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getSizeListL() const
{
  return(sizeListL);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const UInt &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getSizeListG() const
{
  return(sizeListG);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const typename dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::FECARDS &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getFeCards() const
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
typename dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::FECARDS &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getFeCards()
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const Teuchos::RCP<communicator> &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getCommDev() const
{
  assert(commDevLoaded);
  return(commDev);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const Teuchos::RCP<typename dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::MESH1D> &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getGrid1d() const
{
  assert(geometryLoaded);
  return(grid1d);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const Teuchos::RCP<typename dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::CONNECT1D> &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getConnectGrid1d() const
{
  assert(geometryLoaded);
  return(connectGrid1d);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const Teuchos::RCP<typename dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::OPTIONS> &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getOptions() const
{
  assert(optionsLoaded);
  return(options);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const sVect<bool> &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getVertexIsActive() const
{
  assert(startupOk);
  return(vertexIsActive);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const sVect<bool> &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getElementIsActive() const
{
  assert(startupOk);
  return(elementIsActive);
}
    
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const sVect<UInt> &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getNewVertexLid() const
{
  assert(startupOk);
  return(newVertexLid);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const sVect<UInt> &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getNewElementLid() const
{
  assert(startupOk);
  return(newElementLid);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER>
const bool &
dofMapStatic1d<FETYPE,DOFTYPE,ORDER,dms1d_geoIdMode>::
getStartupOk() const
{
  return(startupOk);
}

template<typename F, typename D, dms1d_order O>
ostream & operator<<(ostream & f, const dofMapStatic1d<F,D, O,dms1d_geoIdMode> & V)
{
  f << "numVertices   : " << V.numVerticesL   << " " << V.numVerticesG << endl;
  f << "numElements   : " << V.numElementsL   << " " << V.numElementsG << endl << endl;
  
  f << "numVertexDofs : " << V.numVertexDofsL << " " << V.numVertexDofsG << endl;
  f << "numVolumeDofs : " << V.numVolumeDofsL << " " << V.numVolumeDofsG << endl << endl;
  
  f << "numTotalDofs  : " << V.numDofsL       << " " << V.numDofsG << endl;
  f << "sizeList      : " << V.sizeListL      << " " << V.sizeListG << endl;
  
  return(f);
}

#endif
