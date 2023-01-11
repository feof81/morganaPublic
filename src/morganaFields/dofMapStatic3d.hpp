/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef DOFMAPSTATIC3D_HPP
#define DOFMAPSTATIC3D_HPP

#include "typesInterface.hpp"
#include "traitsBasic.h"

#include "pVectManip.hpp"

#include "morganaGeometry.hpp"
#include "geoMapInterface.hpp"
#include "connect3d.hpp"
#include "mesh3dGlobalManip.hpp"

#include "dofMapStatic3d_options.h"
#include "feStaticDofCard3d.h"



enum dms3d_order {dms3d_vectMajor, dms3d_componentMajor};
enum dms3d_mode  {dms3d_allMode, dms3d_geoIdMode, dms3d_disjointMode};




//! Forward declaration
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER = dms3d_vectMajor, dms3d_mode MODE = dms3d_allMode> class dofMapStatic3d;


//_________________________________________________________________________________________________
// ALL MODE
//-------------------------------------------------------------------------------------------------

/*! Support for the mapping of 3d static finite elements -> \c allMode */
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
class dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE   PMAPTYPE;
    typedef typename FETYPE::GEOSHAPE          GEOSHAPE;
    typedef typename FETYPE::DOFCARD           DOFCARD;
    typedef typename FETYPE::FECARD            FECARD;
    
    typedef mesh3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH3D;
    typedef connect3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT3D;
    typedef geoMapInterface<GEOSHAPE>             GEOMAPINTERFACE;
    typedef pVect<FECARD,PMAPTYPE>                FECARDS;
    typedef dofMapStatic3d_options                OPTIONS;
    
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
    Teuchos::RCP<MESH3D>        grid3d;
    Teuchos::RCP<CONNECT3D>     connectGrid3d;
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
    UInt numFacesL,    numFacesG;
    UInt numElementsL, numElementsG;
    
    UInt numVertexDofsL, numVertexDofsG;
    UInt numEdgesDofsL,  numEdgesDofsG;
    UInt numFacesDofsL,  numFacesDofsG;
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
    dofMapStatic3d();
    dofMapStatic3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    dofMapStatic3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
    dofMapStatic3d(const dofMapStatic3d & DofMap);
    dofMapStatic3d operator=(const dofMapStatic3d & DofMap);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    void setGeometry(MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
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
    const DOFMAP                     & getDofMap()  const;
    const LISTMAP                    & getListMap() const;
    const FECARDS                    & getFeCards() const;
          FECARDS                    & getFeCards();
    const Teuchos::RCP<communicator> & getCommDev() const;
    const Teuchos::RCP<MESH3D>       & getGrid3d()  const;
    const Teuchos::RCP<CONNECT3D>    & getConnectGrid3d() const;
    const Teuchos::RCP<OPTIONS>      & getOptions() const;
    const bool                       & getStartupOk() const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    template<typename F, typename D, dms3d_order O>
    friend ostream & operator<<(ostream & f, const dofMapStatic3d<F,D, O,dms3d_allMode> & V);
    //@}
};


template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
dofMapStatic3d()
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = false;
  geometryLoaded = false;
  optionsLoaded  = false;
  startupOk      = false;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
dofMapStatic3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = CommDev;
  grid3d        = Grid3d;
  connectGrid3d = ConnedGrid3d;
  
  feCards.resize(grid3d->getNumElements());
  feCards.setMap(grid3d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
dofMapStatic3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = Teuchos::rcpFromRef(CommDev);
  grid3d        = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid3d = Teuchos::rcp(new CONNECT3D(ConnedGrid3d));
  
  feCards.resize(grid3d->getNumElements());
  feCards.setMap(grid3d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
dofMapStatic3d(const dofMapStatic3d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid3d        = DofMap.grid3d;
  connectGrid3d = DofMap.connectGrid3d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;

  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numEdgesL    = DofMap.numEdgesL;
  numEdgesG    = DofMap.numEdgesG;
  numFacesL    = DofMap.numFacesL;
  numFacesG    = DofMap.numFacesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numEdgesDofsL  = DofMap.numEdgesDofsL;
  numEdgesDofsG  = DofMap.numEdgesDofsG;
  numFacesDofsL  = DofMap.numFacesDofsL;
  numFacesDofsG  = DofMap.numFacesDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  sizeListL = DofMap.sizeListL;
  sizeListG = DofMap.sizeListG;
    
  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
operator=(const dofMapStatic3d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid3d        = DofMap.grid3d;
  connectGrid3d = DofMap.connectGrid3d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;

  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numEdgesL    = DofMap.numEdgesL;
  numEdgesG    = DofMap.numEdgesG;
  numFacesL    = DofMap.numFacesL;
  numFacesG    = DofMap.numFacesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numEdgesDofsL  = DofMap.numEdgesDofsL;
  numEdgesDofsG  = DofMap.numEdgesDofsG;
  numFacesDofsL  = DofMap.numFacesDofsL;
  numFacesDofsG  = DofMap.numFacesDofsG;
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

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded  = true;
  commDev        = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded  = true;
  commDev        = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d)
{
  geometryLoaded = true;
  grid3d         = Grid3d;
  connectGrid3d  = ConnedGrid3d;
  
  feCards.resize(grid3d->getNumElements());
  feCards.setMap(grid3d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setGeometry(MESH3D & Grid3d, CONNECT3D & ConnedGrid3d)
{
  geometryLoaded = true; 
  grid3d         = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnedGrid3d));
  
  feCards.resize(grid3d->getNumElements());
  feCards.setMap(grid3d->getElements().getRowMap());
  feCards.updateFinder();
}



//_________________________________________________________________________________________________
// OPTIONS LOADING AND STARTUP
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setFeCardL(const UInt & lid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(lid <= feCards.sizeL());
  feCards.getL(lid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setFeCardG(const UInt & gid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(feCards.isG(gid));
  feCards.getG(gid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setFeCards(const FECARDS & FECards)
{
  assert(geometryLoaded);
  assert(feCards.size() == FECards.size());
  feCards = FECards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setOptions(const Teuchos::RCP<OPTIONS> & Options)
{
  optionsLoaded = true;
  options = Options;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setOptions(OPTIONS & Options)
{
  optionsLoaded = true;
  options = Teuchos::rcp(new OPTIONS(Options));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
startup()
{
  //Assert and flags_______________________________________________________________________________
  assert(commDevLoaded);
  assert(geometryLoaded);
  assert(optionsLoaded);
  
  assert(grid3d->getElements().colIsLocal());
  assert(grid3d->getFaces().colIsLocal());
  assert(grid3d->getEdges().colIsLocal());
  
  mesh3dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> manipulator(commDev);
  assert(manipulator.check(grid3d));
  
  startupOk = true;
  
  
  //Count the geotypes_____________________________________________________________________________ 
  numVerticesL = grid3d->getNumVertices();
  numEdgesL    = grid3d->getNumEdges();
  numFacesL    = grid3d->getNumFaces();
  numElementsL = grid3d->getNumElements();
  
  numVerticesG = manipulator.getNumGlobalVertices(grid3d);
  numEdgesG    = manipulator.getNumGlobalEdges(grid3d);
  numFacesG    = manipulator.getNumGlobalFaces(grid3d);
  numElementsG = manipulator.getNumGlobalElements(grid3d);
  
  //Count the dofs
  numVertexDofsL = FETYPE::dofPerVertex * numVerticesL;
  numEdgesDofsL  = FETYPE::dofPerEdge   * numEdgesL;
  numFacesDofsL  = FETYPE::dofPerFace   * numFacesL;
  numVolumeDofsL = FETYPE::dofPerVolume * numElementsL;
  
  numVertexDofsG = FETYPE::dofPerVertex * numVerticesG;
  numEdgesDofsG  = FETYPE::dofPerEdge   * numEdgesG;
  numFacesDofsG  = FETYPE::dofPerFace   * numFacesG;
  numVolumeDofsG = FETYPE::dofPerVolume * numElementsG;
  
  
  //Total number of dofs and lists_________________________________________________________________
  numDofsL = numVertexDofsL + numEdgesDofsL + numFacesDofsL + numVolumeDofsL;
  numDofsG = numVertexDofsG + numEdgesDofsG + numFacesDofsG + numVolumeDofsG;
  
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
      mapItem = grid3d->getNodes().getMapL(i);
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
      mapItem = grid3d->getEdges().getRowMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numEdgesG) + numVertexDofsG );
      
      dofMap(k) = mapItem;
      ++k;
    }
  }
  
  for(UInt lev=1; lev <= FETYPE::dofPerFace; ++lev)  //Faces
  {
    for(UInt i=1; i <= numFacesL; ++i)
    {
      mapItem = grid3d->getFaces().getRowMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numFacesG) + numVertexDofsG + numEdgesDofsG );
      
      dofMap(k) = mapItem;
      ++k;
    }
  }
  
  for(UInt lev=1; lev <= FETYPE::dofPerVolume; ++lev)  //Volume
  {
    for(UInt i=1; i <= numElementsL; ++i)
    {
      mapItem = grid3d->getElements().getRowMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numElementsG) + numVertexDofsG + numEdgesDofsG + numFacesDofsG );
      
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
      
      lid = (block * (mapItem.getLid() - 1) + K) * (ORDER == dms3d_vectMajor)  + (mapItem.getLid() + (K-1) * numDofsL) * (ORDER == dms3d_componentMajor);
      gid = (block * (mapItem.getGid() - 1) + K) * (ORDER == dms3d_vectMajor)  + (mapItem.getGid() + (K-1) * numDofsG) * (ORDER == dms3d_componentMajor);
      
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

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
bool
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
isActive(const DOFCARD & dofCard) const
{
  assert( (&dofCard) == (&dofCard) );
  return(true);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofL - 1) + K)  * (ORDER == dms3d_vectMajor) +
          (dofL  + (K-1) * numDofsL) * (ORDER == dms3d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofG - 1) + K)   * (ORDER == dms3d_vectMajor) +
          (dofG  + (K-1) * numDofsG) * (ORDER == dms3d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
mapDofL(const DOFCARD & dofCard) const
{
  assert(startupOk);
  
  switch(dofCard.getGeoType())
  {
    case VERTEX :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVertex);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numVertices);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());     
      return( grid3d->getElements().getCid_LL(dofCard.getLocalElId(), dofCard.getLocalId()) + (dofCard.getLevel() - 1) * numVerticesL);
    
    case EDGE   :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerEdge);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numEdges);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      return( connectGrid3d->getElementToEdge(dofCard.getLocalElId(), dofCard.getLocalId()) + (dofCard.getLevel() - 1) * numEdgesL  + numVertexDofsL );
    
    case FACE   :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerFace);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numFaces);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      return( connectGrid3d->getElementToFace(dofCard.getLocalElId(), dofCard.getLocalId()) + (dofCard.getLevel() - 1) * numFacesL  + numVertexDofsL + numEdgesDofsL );
    
    case VOLUME :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVolume);
      assert(dofCard.getLocalId()   == 1);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      return(dofCard.getLocalElId() + (dofCard.getLevel() - 1) * numElementsL + numVertexDofsL + numEdgesDofsL + numFacesDofsL );
      
    default :
      assert(1==2);
  }
  
  return(0);             
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
mapDofG(const DOFCARD & dofCard) const
{
  assert(startupOk);
  return(dofMap(mapDofL(dofCard)).getGid());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofL = mapDofL(dofCard);
  return(dofToListL(dofL,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofG = mapDofG(dofCard);
  return(dofToListG(dofG,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::DOFMAP &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getDofMap() const
{
  assert(startupOk);
  return(dofMap);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::LISTMAP &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getListMap() const
{
  assert(startupOk);
  return(listMap);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getNumDofsL() const
{
  assert(startupOk);
  return(numDofsL);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getNumDofsG() const
{
  assert(startupOk);
  return(numDofsG);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getSizeListL() const
{
  assert(startupOk);
  return(sizeListL);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getSizeListG() const
{
  assert(startupOk);
  return(sizeListG);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::FECARDS &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getFeCards() const
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::FECARDS &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getFeCards()
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<communicator> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getCommDev() const
{
  assert(commDevLoaded);
  return(commDev);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::MESH3D> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getGrid3d() const
{
  assert(geometryLoaded);
  return(grid3d);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::CONNECT3D> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getConnectGrid3d() const
{
  assert(geometryLoaded);
  return(connectGrid3d);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::OPTIONS> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getOptions() const
{
  assert(optionsLoaded);
  return(options);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const bool &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getStartupOk() const
{
  return(startupOk);
}

template<typename F, typename D, dms3d_order O>
ostream & operator<<(ostream & f, const dofMapStatic3d<F,D, O,dms3d_allMode> & V)
{
  f << "numVertices   : " << V.numVerticesL   << " " << V.numVerticesG << endl;
  f << "numEdges      : " << V.numEdgesL      << " " << V.numEdgesG << endl;
  f << "numFaces      : " << V.numFacesL      << " " << V.numFacesG << endl;
  f << "numElements   : " << V.numElementsL   << " " << V.numElementsG << endl << endl;
  
  f << "numVertexDofs : " << V.numVertexDofsL << " " << V.numVertexDofsG << endl;
  f << "numEdgesDofs  : " << V.numEdgesDofsL  << " " << V.numEdgesDofsG << endl;
  f << "numFacesDofs  : " << V.numFacesDofsL  << " " << V.numFacesDofsG << endl;
  f << "numVolumeDofs : " << V.numVolumeDofsL << " " << V.numVolumeDofsG << endl << endl;
  
  f << "numTotalDofs  : " << V.numDofsL       << " " << V.numDofsG << endl;
  f << "sizeList      : " << V.sizeListL      << " " << V.sizeListG << endl;
  
  return(f);
}



//_________________________________________________________________________________________________
// GEO ID MODE
//-------------------------------------------------------------------------------------------------

/*! Support for the mapping of 3d static finite elements -> \c geoIdMode */
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
class dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE   PMAPTYPE;
    typedef typename FETYPE::GEOSHAPE          GEOSHAPE;
    typedef typename FETYPE::DOFCARD           DOFCARD;
    typedef typename FETYPE::FECARD            FECARD;
    
    typedef mesh3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH3D;
    typedef connect3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT3D;
    typedef geoMapInterface<GEOSHAPE>             GEOMAPINTERFACE;
    typedef pVect<FECARD,PMAPTYPE>                FECARDS;
    typedef dofMapStatic3d_options                OPTIONS;
    
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
    Teuchos::RCP<MESH3D>        grid3d;
    Teuchos::RCP<CONNECT3D>     connectGrid3d;
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
    sVect<bool> faceIsActive;
    sVect<bool> elementIsActive;
    
    sVect<UInt> newVertexLid;
    sVect<UInt> newEdgeLid;
    sVect<UInt> newFaceLid;
    sVect<UInt> newElementLid;
    
    pVect<UInt,PMAPTYPE>  globVertices;
    pVect<UInt,PMAPTYPE>  globEdges;
    pVect<UInt,PMAPTYPE>  globFaces;
    pVect<UInt,PMAPTYPE>  globElements;
    
    UInt numVerticesL, numVerticesG;
    UInt numEdgesL,    numEdgesG;
    UInt numFacesL,    numFacesG;
    UInt numElementsL, numElementsG;
    
    UInt numVertexDofsL, numVertexDofsG;
    UInt numEdgesDofsL,  numEdgesDofsG;
    UInt numFacesDofsL,  numFacesDofsG;
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
    dofMapStatic3d();
    dofMapStatic3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    dofMapStatic3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
    dofMapStatic3d(const dofMapStatic3d & DofMap);
    dofMapStatic3d operator=(const dofMapStatic3d & DofMap);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    void setGeometry(MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
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
    const Teuchos::RCP<MESH3D>       & getGrid3d() const;
    const Teuchos::RCP<CONNECT3D>    & getConnectGrid3d() const;
    const Teuchos::RCP<OPTIONS>      & getOptions() const;
    const bool                       & getStartupOk() const;
    //@}
    
    /*! @name Extra dump functions */ //@{
  public:
    const sVect<bool> & getVertexIsActive() const;
    const sVect<bool> & getEdgeIsActive() const;
    const sVect<bool> & getFaceIsActive() const;
    const sVect<bool> & getElementIsActive() const;
    
    const sVect<UInt> & getNewVertexLid() const;
    const sVect<UInt> & getNewEdgeLid() const;
    const sVect<UInt> & getNewFaceLid() const;
    const sVect<UInt> & getNewElementLid() const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    template<typename F, typename D, dms3d_order O>
    friend ostream & operator<<(ostream & f, const dofMapStatic3d<F,D, O,dms3d_allMode> & V);
    //@}
};



template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
dofMapStatic3d()
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = false;
  geometryLoaded = false;
  optionsLoaded  = false;
  startupOk      = false;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
dofMapStatic3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = CommDev;
  grid3d        = Grid3d;
  connectGrid3d = ConnedGrid3d;
  
  feCards.resize(grid3d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
dofMapStatic3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = Teuchos::rcpFromRef(CommDev);
  grid3d        = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid3d = Teuchos::rcp(new CONNECT3D(ConnedGrid3d));
  
  feCards.resize(grid3d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
dofMapStatic3d(const dofMapStatic3d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid3d        = DofMap.grid3d;
  connectGrid3d = DofMap.connectGrid3d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;
     
  vertexIsActive  = DofMap.vertexIsActive;
  edgeIsActive    = DofMap.edgeIsActive;
  faceIsActive    = DofMap.faceIsActive;
  elementIsActive = DofMap.elementIsActive;
   
  newVertexLid  = DofMap.newVertexLid;
  newEdgeLid    = DofMap.newEdgeLid; 
  newFaceLid    = DofMap.newFaceLid;
  newElementLid = DofMap.newElementLid;
    
  globVertices = DofMap.globVertices;
  globEdges    = DofMap.globEdges;
  globFaces    = DofMap.globFaces;
  globElements = DofMap.globElements;
    
  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numEdgesL    = DofMap.numEdgesL;
  numEdgesG    = DofMap.numEdgesG;
  numFacesL    = DofMap.numFacesL;
  numFacesG    = DofMap.numFacesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numEdgesDofsL  = DofMap.numEdgesDofsL;
  numEdgesDofsG  = DofMap.numEdgesDofsG;
  numFacesDofsL  = DofMap.numFacesDofsL;
  numFacesDofsG  = DofMap.numFacesDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  sizeListL = DofMap.sizeListL;
  sizeListG = DofMap.sizeListG;
    
  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
operator=(const dofMapStatic3d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid3d        = DofMap.grid3d;
  connectGrid3d = DofMap.connectGrid3d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;
     
  vertexIsActive  = DofMap.vertexIsActive;
  edgeIsActive    = DofMap.edgeIsActive;
  faceIsActive    = DofMap.faceIsActive;
  elementIsActive = DofMap.elementIsActive;
   
  newVertexLid  = DofMap.newVertexLid;
  newEdgeLid    = DofMap.newEdgeLid; 
  newFaceLid    = DofMap.newFaceLid;
  newElementLid = DofMap.newElementLid;
    
  globVertices = DofMap.globVertices;
  globEdges    = DofMap.globEdges;
  globFaces    = DofMap.globFaces;
  globElements = DofMap.globElements;
    
  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numEdgesL    = DofMap.numEdgesL;
  numEdgesG    = DofMap.numEdgesG;
  numFacesL    = DofMap.numFacesL;
  numFacesG    = DofMap.numFacesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numEdgesDofsL  = DofMap.numEdgesDofsL;
  numEdgesDofsG  = DofMap.numEdgesDofsG;
  numFacesDofsL  = DofMap.numFacesDofsL;
  numFacesDofsG  = DofMap.numFacesDofsG;
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

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded  = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded  = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d)
{
  geometryLoaded = true;
  grid3d         = Grid3d;
  connectGrid3d  = ConnedGrid3d;
  
  feCards.resize(grid3d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setGeometry(MESH3D & Grid3d, CONNECT3D & ConnedGrid3d)
{
  geometryLoaded = true; 
  grid3d         = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnedGrid3d));
  
  feCards.resize(grid3d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setFeCardL(const UInt & lid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(lid <= feCards.sizeL());
  feCards.getL(lid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setFeCardG(const UInt & gid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(feCards.isG(gid));
  feCards.getG(gid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setFeCards(const FECARDS & FECards)
{
  assert(geometryLoaded);
  assert(feCards.size() == FECards.size());
  feCards = FECards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setOptions(const Teuchos::RCP<OPTIONS> & Options)
{
  optionsLoaded = true;
  options = Options;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setOptions(OPTIONS & Options)
{
  optionsLoaded = true;
  options = Teuchos::rcp(new OPTIONS(Options));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
startup()
{
  assert(commDevLoaded);
  assert(optionsLoaded);
  
  assert(grid3d->getElements().colIsLocal());
  assert(grid3d->getFaces().colIsLocal());
  assert(grid3d->getEdges().colIsLocal());
  
  mesh3dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> manipulator(commDev);
  assert(manipulator.check(grid3d));
  
  startupOk = true;
  
  //Alloc__________________________________________________________________________________________
  vertexIsActive.resize(  grid3d->getNumVertices() );
  edgeIsActive.resize(    grid3d->getNumEdges()    );
  faceIsActive.resize(    grid3d->getNumFaces()    );
  elementIsActive.resize( grid3d->getNumElements() );
   
  newVertexLid.resize(  grid3d->getNumVertices() );
  newEdgeLid.resize(    grid3d->getNumEdges()    );
  newFaceLid.resize(    grid3d->getNumFaces()    );
  newElementLid.resize( grid3d->getNumElements() );
  
  for(UInt i=1; i <= grid3d->getNumVertices(); ++i)
  { vertexIsActive(i) = false; }
  
  for(UInt i=1; i <= grid3d->getNumEdges(); ++i)
  { edgeIsActive(i) = false; }
  
  for(UInt i=1; i <= grid3d->getNumFaces(); ++i)
  { faceIsActive(i) = false;}
  
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  { elementIsActive(i) = false; }

  
  //Identify the active items______________________________________________________________________
  UInt geoId;
  
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    geoId = grid3d->getElementL(i).getGeoId();
    
    if(options->isGeoId(geoId))
    {
      for(UInt j=1; j <= GEOSHAPE::numVertices; ++j)
      { vertexIsActive(grid3d->getElementL(i).getCid(j)) = true; }
      
      for(UInt j=1; j <= GEOSHAPE::numEdges; ++j)
      { edgeIsActive(connectGrid3d->getElementToEdge(i,j)) = true; }
      
      for(UInt j=1; j <= GEOSHAPE::numFaces; ++j)
      { faceIsActive(connectGrid3d->getElementToFace(i,j)) = true; }
      
      elementIsActive(i) = true;
    }
  }
  
  //Assign new local ids___________________________________________________________________________
  numVerticesL = 0;
  numEdgesL    = 0;
  numFacesL    = 0;
  numElementsL = 0;
  
  globVertices.clear();
  globEdges.clear();
  globFaces.clear();
  globElements.clear();
  
  PMAPTYPE nodeMapItem;
  PMAPTYPE elMapItem;
  
  for(UInt i=1; i <= grid3d->getNumVertices(); ++i)  //Vertices
  {
    if(vertexIsActive(i))
    {
      ++numVerticesL;
      newVertexLid(i) = numVerticesL;
      
      nodeMapItem = grid3d->getNodes().getMapL(i);
      nodeMapItem.setLid(numVerticesL);
      
      globVertices.push_back(nodeMapItem, nodeMapItem.getGid());
    }
  }
  
  for(UInt i=1; i <= grid3d->getNumEdges(); ++i)     //Edges
  {
    if(edgeIsActive(i))
    {
      ++numEdgesL;
      newEdgeLid(i) = numEdgesL;
      
      elMapItem = grid3d->getEdges().getRowMapL(i);
      elMapItem.setLid(numEdgesL);
      
      globEdges.push_back(elMapItem, elMapItem.getGid());
    }
  }
  
  for(UInt i=1; i <= grid3d->getNumFaces(); ++i)     //Faces
  {
    if(faceIsActive(i))
    {
      ++numFacesL;
      newFaceLid(i) = numFacesL;
      
      elMapItem = grid3d->getFaces().getRowMapL(i);
      elMapItem.setLid(numFacesL);
      
      globFaces.push_back(elMapItem, elMapItem.getGid());
    }
  }
  
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)  //Elements
  {
    if(elementIsActive(i))
    {
      ++numElementsL;
      newElementLid(i) = numElementsL;
      
      elMapItem = grid3d->getElements().getRowMapL(i);
      elMapItem.setLid(numElementsL);
      
      globElements.push_back(elMapItem, elMapItem.getGid());
    }
  }
  
  
  //Assign new global ids__________________________________________________________________________
  pVectGlobalManip<UInt,PMAPTYPE> nodeManip(commDev);
  pVectGlobalManip<UInt,PMAPTYPE>   elManip(commDev);
   
  nodeManip.buildGlobalNumbering(globVertices);
  elManip.buildGlobalNumbering(globEdges);
  elManip.buildGlobalNumbering(globFaces);
  elManip.buildGlobalNumbering(globElements);
  
  numVerticesG = nodeManip.sizeG(globVertices);
  numEdgesG    = elManip.sizeG(globEdges);
  numFacesG    = elManip.sizeG(globFaces);
  numElementsG = elManip.sizeG(globElements);
   
  //Count the dofs_________________________________________________________________________________
  numVertexDofsL = FETYPE::dofPerVertex * numVerticesL;
  numEdgesDofsL  = FETYPE::dofPerEdge   * numEdgesL;
  numFacesDofsL  = FETYPE::dofPerFace   * numFacesL;
  numVolumeDofsL = FETYPE::dofPerVolume * numElementsL;
  
  numVertexDofsG = FETYPE::dofPerVertex * numVerticesG;
  numEdgesDofsG  = FETYPE::dofPerEdge   * numEdgesG;
  numFacesDofsG  = FETYPE::dofPerFace   * numFacesG;
  numVolumeDofsG = FETYPE::dofPerVolume * numElementsG;
  
  //Total number of dofs and lists_________________________________________________________________
  numDofsL = numVertexDofsL + numEdgesDofsL + numFacesDofsL + numVolumeDofsL;
  numDofsG = numVertexDofsG + numEdgesDofsG + numFacesDofsG + numVolumeDofsG;
  
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
  
  for(UInt lev=1; lev <= FETYPE::dofPerFace; ++lev)  //Faces
  {
    for(UInt i=1; i <= numFacesL; ++i)
    {
      mapItem = globFaces.getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numFacesG) + numVertexDofsG + numEdgesDofsG );
      
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
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numElementsG) + numVertexDofsG + numEdgesDofsG + numFacesDofsG );
      
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
      
      lid = (block * (mapItem.getLid() - 1) + K) * (ORDER == dms3d_vectMajor)  + (mapItem.getLid() + (K-1) * numDofsL) * (ORDER == dms3d_componentMajor);
      gid = (block * (mapItem.getGid() - 1) + K) * (ORDER == dms3d_vectMajor)  + (mapItem.getGid() + (K-1) * numDofsG) * (ORDER == dms3d_componentMajor);
      
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

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
bool
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
isActive(const DOFCARD & dofCard) const
{
  assert(dofCard.getLocalElId() <= elementIsActive.size());
  return(elementIsActive(dofCard.getLocalElId()));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofL - 1) + K)  * (ORDER == dms3d_vectMajor) +
          (dofL  + (K-1) * numDofsL) * (ORDER == dms3d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofG - 1) + K)   * (ORDER == dms3d_vectMajor) +
          (dofG  + (K-1) * numDofsG) * (ORDER == dms3d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
mapDofL(const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(elementIsActive(dofCard.getLocalElId()));
  
  UInt oldId = 0;
  
  switch(dofCard.getGeoType())
  {
    case VERTEX :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVertex);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numVertices);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = grid3d->getElements().getCid_LL(dofCard.getLocalElId(), dofCard.getLocalId());
      return( newVertexLid(oldId) + (dofCard.getLevel() - 1) * numVerticesL);
    
    case EDGE   :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerEdge);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numEdges);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = connectGrid3d->getElementToEdge(dofCard.getLocalElId(), dofCard.getLocalId());
      return( newEdgeLid(oldId) + (dofCard.getLevel() - 1) * numEdgesL + numVertexDofsL );
    
    case FACE   :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerFace);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numFaces);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = connectGrid3d->getElementToFace(dofCard.getLocalElId(), dofCard.getLocalId());
      return( newFaceLid(oldId) + (dofCard.getLevel() - 1) * numFacesL + numVertexDofsL + numEdgesDofsL );
    
    case VOLUME :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVolume);
      assert(dofCard.getLocalId()   == 1);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = dofCard.getLocalElId();
      return( newElementLid(oldId) + (dofCard.getLevel() - 1) * numElementsL + numVertexDofsL + numEdgesDofsL + numFacesDofsL );
      
    default :
      assert(1==2);
   }
  
  return(0); 
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
mapDofG(const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(elementIsActive(dofCard.getLocalElId()));
  
  return(dofMap(mapDofL(dofCard)).getGid());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofL = mapDofL(dofCard);
  return(dofToListL(dofL,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofG = mapDofG(dofCard);
  return(dofToListG(dofG,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::DOFMAP &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getDofMap() const
{
  return(dofMap);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::LISTMAP &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getListMap() const
{
  return(listMap);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getNumDofsL() const
{
  return(numDofsL);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getNumDofsG() const
{
  return(numDofsG);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getSizeListL() const
{
  return(sizeListL);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getSizeListG() const
{
  return(sizeListG);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::FECARDS &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getFeCards() const
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::FECARDS &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getFeCards()
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<communicator> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getCommDev() const
{
  assert(commDevLoaded);
  return(commDev);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::MESH3D> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getGrid3d() const
{
  assert(geometryLoaded);
  return(grid3d);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::CONNECT3D> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getConnectGrid3d() const
{
  assert(geometryLoaded);
  return(connectGrid3d);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::OPTIONS> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getOptions() const
{
  assert(optionsLoaded);
  return(options);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<bool> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getVertexIsActive() const
{
  assert(startupOk);
  return(vertexIsActive);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<bool> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getEdgeIsActive() const
{
  assert(startupOk);
  return(edgeIsActive);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<bool> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getFaceIsActive() const
{
  assert(startupOk);
  return(faceIsActive);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<bool> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getElementIsActive() const
{
  assert(startupOk);
  return(elementIsActive);
}
    
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<UInt> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getNewVertexLid() const
{
  assert(startupOk);
  return(newVertexLid);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<UInt> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getNewEdgeLid() const
{
  assert(startupOk);
  return(newEdgeLid);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<UInt> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getNewFaceLid() const
{
  assert(startupOk);
  return(newFaceLid);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<UInt> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getNewElementLid() const
{
  assert(startupOk);
  return(newElementLid);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const bool &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getStartupOk() const
{
  return(startupOk);
}

template<typename F, typename D, dms3d_order O>
ostream & operator<<(ostream & f, const dofMapStatic3d<F,D, O,dms3d_geoIdMode> & V)
{
  f << "numVertices   : " << V.numVerticesL   << " " << V.numVerticesG << endl;
  f << "numEdges      : " << V.numEdgesL      << " " << V.numEdgesG << endl;
  f << "numFaces      : " << V.numFacesL      << " " << V.numFacesG << endl;
  f << "numElements   : " << V.numElementsL   << " " << V.numElementsG << endl << endl;
  
  f << "numVertexDofs : " << V.numVertexDofsL << " " << V.numVertexDofsG << endl;
  f << "numEdgesDofs  : " << V.numEdgesDofsL  << " " << V.numEdgesDofsG << endl;
  f << "numFacesDofs  : " << V.numFacesDofsL  << " " << V.numFacesDofsG << endl;
  f << "numVolumeDofs : " << V.numVolumeDofsL << " " << V.numVolumeDofsG << endl << endl;
  
  f << "numTotalDofs  : " << V.numDofsL       << " " << V.numDofsG << endl;
  f << "sizeList      : " << V.sizeListL      << " " << V.sizeListG << endl;
  
  return(f);
}


//_________________________________________________________________________________________________
// DISJOINT MODE
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
class dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE   PMAPTYPE;
    typedef typename FETYPE::GEOSHAPE          GEOSHAPE;
    typedef typename FETYPE::DOFCARD           DOFCARD;
    typedef typename FETYPE::FECARD            FECARD;
    
    typedef mesh3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH3D;
    typedef connect3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT3D;
    typedef geoMapInterface<GEOSHAPE>             GEOMAPINTERFACE;
    typedef pVect<FECARD,PMAPTYPE>                FECARDS;
    typedef dofMapStatic3d_options                OPTIONS;
    
    typedef pMap<PMAPTYPE> DOFMAP;
    typedef pMap<PMAPTYPE> LISTMAP;
    
    typedef sVect<bool>          SVECT_BOOL;
    typedef sVect<UInt>          SVECT_UINT;
    typedef sVect<SVECT_BOOL>    ARRAY_BOOL;
    typedef sVect<SVECT_UINT>    ARRAY_UINT;
    typedef pVect<UInt,PMAPTYPE> PVECT_UINT;
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
    Teuchos::RCP<MESH3D>        grid3d;
    Teuchos::RCP<CONNECT3D>     connectGrid3d;
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
    sVect<UInt> elementsToBlock;
    
    ARRAY_BOOL vertexIsActive;
    ARRAY_BOOL edgeIsActive;
    ARRAY_BOOL faceIsActive;
    ARRAY_BOOL elementIsActive;
      
    ARRAY_UINT newVertexLid;
    ARRAY_UINT newEdgeLid;
    ARRAY_UINT newFaceLid;
    ARRAY_UINT newElementLid;
    
    sVect<PVECT_UINT> globVertices;
    sVect<PVECT_UINT> globEdges;
    sVect<PVECT_UINT> globFaces;
    sVect<PVECT_UINT> globElements;
    
    SVECT_UINT numVerticesL, numVerticesG;
    SVECT_UINT numEdgesL,    numEdgesG;
    SVECT_UINT numFacesL,    numFacesG;
    SVECT_UINT numElementsL, numElementsG;
    
    SVECT_UINT numVertexDofsL, numVertexDofsG;
    SVECT_UINT numEdgesDofsL,  numEdgesDofsG;
    SVECT_UINT numFacesDofsL,  numFacesDofsG;
    SVECT_UINT numVolumeDofsL, numVolumeDofsG;
    
    SVECT_UINT numDofsL, numDofsG;
    SVECT_UINT offsetL,  offsetG;
    
    UInt numDofsTotL, numDofsTotG;
    UInt sizeListL, sizeListG;
    //@}
    
    /*! @name Internal Maps */ //@{
  public:
    DOFMAP  dofMap;
    LISTMAP listMap;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    dofMapStatic3d();
    dofMapStatic3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    dofMapStatic3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
    dofMapStatic3d(const dofMapStatic3d & DofMap);
    dofMapStatic3d operator=(const dofMapStatic3d & DofMap);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    void setGeometry(MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
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
    const Teuchos::RCP<MESH3D>       & getGrid3d() const;
    const Teuchos::RCP<CONNECT3D>    & getConnectGrid3d() const;
    const Teuchos::RCP<OPTIONS>      & getOptions() const;
    const bool                       & getStartupOk() const;
    //@}
    
    /*! @name Extra dump functions */ //@{
  public:
    const UInt & getNumBlocks() const;
    
    const SVECT_UINT & getNewVertexLid(const UInt & k) const;
    const SVECT_UINT & getNewEdgeLid(const UInt & k) const;
    const SVECT_UINT & getNewFaceLid(const UInt & k) const;
    const SVECT_UINT & getNewElementLid(const UInt & k) const;
    
    const PVECT_UINT & getGlobVertices(const UInt & k) const;
    const PVECT_UINT & getGlobEdges(const UInt & k) const;
    const PVECT_UINT & getGlobFaces(const UInt & k) const;
    const PVECT_UINT & getGlobElements(const UInt & k) const;
        
    const UInt & getNumVerticesL(const UInt & k) const;
    const UInt & getNumEdgesL(const UInt & k) const;
    const UInt & getNumFacesL(const UInt & k) const;
    const UInt & getNumElementsL(const UInt & k) const;
    
    const UInt & getNumVerticesG(const UInt & k) const;
    const UInt & getNumEdgesG(const UInt & k) const;
    const UInt & getNumFacesG(const UInt & k) const;
    const UInt & getNumElementsG(const UInt & k) const;
    
    const UInt & getNumVertexDofsL(const UInt & k) const;
    const UInt & getNumEdgesDofsL(const UInt & k) const;
    const UInt & getNumFacesDofsL(const UInt & k) const;
    const UInt & getNumVolumeDofsL(const UInt & k) const;
    
    const UInt & getNumVertexDofsG(const UInt & k) const;
    const UInt & getNumEdgesDofsG(const UInt & k) const;
    const UInt & getNumFacesDofsG(const UInt & k) const;
    const UInt & getNumVolumeDofsG(const UInt & k) const;
    
    const UInt & getOffsetL(const UInt & k) const;
    const UInt & getOffsetG(const UInt & k) const;
    //@}
    
    /*! @name Printout */ //@{
  /*public:
    template<typename F, typename D, dms3d_order O>
    friend ostream & operator<<(ostream & f, const dofMapStatic3d<F,D, O,dms3d_allMode> & V);*/
    //@}
};


template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
dofMapStatic3d()
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = false;
  geometryLoaded = false;
  optionsLoaded  = false;
  startupOk      = false;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
dofMapStatic3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = CommDev;
  grid3d        = Grid3d;
  connectGrid3d = ConnedGrid3d;
  
  feCards.resize(grid3d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
dofMapStatic3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = Teuchos::rcpFromRef(CommDev);
  grid3d        = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid3d = Teuchos::rcp(new CONNECT3D(ConnedGrid3d));
  
  feCards.resize(grid3d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
dofMapStatic3d(const dofMapStatic3d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid3d        = DofMap.grid3d;
  connectGrid3d = DofMap.connectGrid3d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;
  
  elementsToBlock = DofMap.elementsToBlock;
  
  vertexIsActive  = DofMap.vertexIsActive;
  edgeIsActive    = DofMap.edgeIsActive;
  faceIsActive    = DofMap.faceIsActive;
  elementIsActive = DofMap.elementIsActive;
      
  newVertexLid  = DofMap.newVertexLid;
  newEdgeLid    = DofMap.newEdgeLid;
  newFaceLid    = DofMap.newFaceLid;
  newElementLid = DofMap.newElementLid;
    
  globVertices = DofMap.globVertices;
  globEdges    = DofMap.globEdges;
  globFaces    = DofMap.globFaces;
  globElements = DofMap.globElements;
  
  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numEdgesL    = DofMap.numEdgesL;
  numEdgesG    = DofMap.numEdgesG;
  numFacesL    = DofMap.numFacesL;
  numFacesG    = DofMap.numFacesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numEdgesDofsL  = DofMap.numEdgesDofsL;
  numEdgesDofsG  = DofMap.numEdgesDofsG;
  numFacesDofsL  = DofMap.numFacesDofsL;
  numFacesDofsG  = DofMap.numFacesDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  
  offsetL = DofMap.offsetL;
  offsetG = DofMap.offsetG;
  
  numDofsTotL = DofMap.numDofsTotL;
  numDofsTotG = DofMap.numDofsTotG;
  sizeListL   = DofMap.sizeListL;
  sizeListG   = DofMap.sizeListG;
    
  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
operator=(const dofMapStatic3d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid3d        = DofMap.grid3d;
  connectGrid3d = DofMap.connectGrid3d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;
  
  elementsToBlock = DofMap.elementsToBlock;
  
  vertexIsActive  = DofMap.vertexIsActive;
  edgeIsActive    = DofMap.edgeIsActive;
  faceIsActive    = DofMap.faceIsActive;
  elementIsActive = DofMap.elementIsActive;
      
  newVertexLid  = DofMap.newVertexLid;
  newEdgeLid    = DofMap.newEdgeLid;
  newFaceLid    = DofMap.newFaceLid;
  newElementLid = DofMap.newElementLid;
    
  globVertices = DofMap.globVertices;
  globEdges    = DofMap.globEdges;
  globFaces    = DofMap.globFaces;
  globElements = DofMap.globElements;
  
  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numEdgesL    = DofMap.numEdgesL;
  numEdgesG    = DofMap.numEdgesG;
  numFacesL    = DofMap.numFacesL;
  numFacesG    = DofMap.numFacesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numEdgesDofsL  = DofMap.numEdgesDofsL;
  numEdgesDofsG  = DofMap.numEdgesDofsG;
  numFacesDofsL  = DofMap.numFacesDofsL;
  numFacesDofsG  = DofMap.numFacesDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  
  offsetL = DofMap.offsetL;
  offsetG = DofMap.offsetG;
  
  numDofsTotL = DofMap.numDofsTotL;
  numDofsTotG = DofMap.numDofsTotG;
  sizeListL   = DofMap.sizeListL;
  sizeListG   = DofMap.sizeListG;
    
  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
  
  return(*this);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded  = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded  = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d)
{
  geometryLoaded = true;
  grid3d         = Grid3d;
  connectGrid3d  = ConnedGrid3d;
  
  feCards.resize(grid3d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
setGeometry(MESH3D & Grid3d, CONNECT3D & ConnedGrid3d)
{
  geometryLoaded = true; 
  grid3d         = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnedGrid3d));
  
  feCards.resize(grid3d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
setFeCardL(const UInt & lid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(lid <= feCards.sizeL());
  feCards.getL(lid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
setFeCardG(const UInt & gid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(feCards.isG(gid));
  feCards.getG(gid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
setFeCards(const FECARDS & FECards)
{
  assert(geometryLoaded);
  assert(feCards.size() == FECards.size());
  feCards = FECards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
setOptions(const Teuchos::RCP<OPTIONS> & Options)
{
  optionsLoaded = true;
  options = Options;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
setOptions(OPTIONS & Options)
{
  optionsLoaded = true;
  options = Teuchos::rcp(new OPTIONS(Options));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
startup()
{
  assert(commDevLoaded);
  assert(optionsLoaded);
  
  assert(grid3d->getElements().colIsLocal());
  assert(grid3d->getFaces().colIsLocal());
  assert(grid3d->getEdges().colIsLocal());
  
  mesh3dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> manipulator(commDev);
  assert(manipulator.check(grid3d));
  
  startupOk = true;
  
  //Alloc------------------------------------------------------------------------------------------
  UInt geoId;
  UInt numBlocks = options->getNumBlocks();
  
  elementsToBlock.resize(grid3d->getNumElements());
  
  vertexIsActive.resize(numBlocks);
  edgeIsActive.resize(numBlocks);
  faceIsActive.resize(numBlocks);
  elementIsActive.resize(numBlocks);
  
  newVertexLid.resize(numBlocks);
  newEdgeLid.resize(numBlocks);
  newFaceLid.resize(numBlocks);
  newElementLid.resize(numBlocks);
  
  globVertices.resize(numBlocks);
  globEdges.resize(numBlocks);
  globFaces.resize(numBlocks);
  globElements.resize(numBlocks);
  
  numVerticesL.resize(numBlocks);
  numEdgesL.resize(numBlocks);
  numFacesL.resize(numBlocks);
  numElementsL.resize(numBlocks);
  
  numVerticesG.resize(numBlocks);
  numEdgesG.resize(numBlocks);
  numFacesG.resize(numBlocks);
  numElementsG.resize(numBlocks);
  
  numVertexDofsL.resize(numBlocks);
  numEdgesDofsL.resize(numBlocks);
  numFacesDofsL.resize(numBlocks);
  numVolumeDofsL.resize(numBlocks);
  
  numVertexDofsG.resize(numBlocks);
  numEdgesDofsG.resize(numBlocks);
  numFacesDofsG.resize(numBlocks);
  numVolumeDofsG.resize(numBlocks);
  
  offsetL.resize(numBlocks);
  offsetG.resize(numBlocks);
  
  numDofsL.resize(numBlocks);
  numDofsG.resize(numBlocks);
  
  for(UInt k=1; k <= numBlocks; ++k)
  {
    vertexIsActive(k).resize(grid3d->getNumVertices());
    edgeIsActive(k).resize(grid3d->getNumEdges());
    faceIsActive(k).resize(grid3d->getNumFaces());
    elementIsActive(k).resize(grid3d->getNumElements());
      
    newVertexLid(k).resize(grid3d->getNumVertices());
    newEdgeLid(k).resize(grid3d->getNumEdges());
    newFaceLid(k).resize(grid3d->getNumFaces());
    newElementLid(k).resize(grid3d->getNumElements());
    
    for(UInt i=1; i <= grid3d->getNumVertices(); ++i)
    { newVertexLid(k)(i) = 0; }
    
    for(UInt i=1; i <= grid3d->getNumEdges(); ++i)
    { newEdgeLid(k)(i) = 0; }
    
    for(UInt i=1; i <= grid3d->getNumFaces(); ++i)
    { newFaceLid(k)(i) = 0; }
    
    for(UInt i=1; i <= grid3d->getNumElements(); ++i)
    { newElementLid(k)(i) = 0; }
  }
  
  //Identify the active items----------------------------------------------------------------------
  for(UInt k=1; k <= numBlocks; ++k)
  {
    for(UInt i=1; i <= grid3d->getNumElements(); ++i)
    {
      geoId = grid3d->getElementL(i).getGeoId();
    
      if(options->isBlockGeoId(k,geoId))
      {
        elementsToBlock(i) = k;
      
        for(UInt j=1; j <= GEOSHAPE::numVertices; ++j)
        { vertexIsActive(k)(grid3d->getElementL(i).getCid(j)) = true; }
      
        for(UInt j=1; j <= GEOSHAPE::numEdges; ++j)
        { edgeIsActive(k)(connectGrid3d->getElementToEdge(i,j)) = true; }
      
        for(UInt j=1; j <= GEOSHAPE::numFaces; ++j)
        { faceIsActive(k)(connectGrid3d->getElementToFace(i,j)) = true; }
      
        elementIsActive(k)(i) = true;
      }
    }
  }
  
  //Assign new local ids---------------------------------------------------------------------------
  PMAPTYPE nodeMapItem;
  PMAPTYPE elMapItem;
  
  for(UInt k=1; k <= numBlocks; ++k)
  {
    numVerticesL(k) = 0;
    numEdgesL(k)    = 0;
    numFacesL(k)    = 0;
    numElementsL(k) = 0;
    
    globVertices(k).clear();
    globEdges(k).clear();
    globFaces(k).clear();
    globElements(k).clear();
    
    for(UInt i=1; i <= grid3d->getNumVertices(); ++i)  //Vertices
    {
      if(vertexIsActive(k)(i))
      {
        ++numVerticesL(k);
        newVertexLid(k)(i) = numVerticesL(k);
      
        nodeMapItem = grid3d->getNodes().getMapL(i);
        nodeMapItem.setLid(numVerticesL(k));
      
        globVertices(k).push_back(nodeMapItem, nodeMapItem.getGid());
      }
    }
    
    for(UInt i=1; i <= grid3d->getNumEdges(); ++i)     //Edges
    {
      if(edgeIsActive(k)(i))
      {
        ++numEdgesL(k);
        newEdgeLid(k)(i) = numEdgesL(k);
      
        elMapItem = grid3d->getEdges().getRowMapL(i);
        elMapItem.setLid(numEdgesL(k));
      
        globEdges(k).push_back(elMapItem, elMapItem.getGid());
      }
    }
    
    for(UInt i=1; i <= grid3d->getNumFaces(); ++i)     //Faces
    {
      if(faceIsActive(k)(i))
      {
        ++numFacesL(k);
        newFaceLid(k)(i) = numFacesL(k);
      
        elMapItem = grid3d->getFaces().getRowMapL(i);
        elMapItem.setLid(numFacesL(k));
      
        globFaces(k).push_back(elMapItem, elMapItem.getGid());
      }
    }
    
    for(UInt i=1; i <= grid3d->getNumElements(); ++i)  //Elements
    {
      if(elementIsActive(k)(i))
      {
        ++numElementsL(k);
        newElementLid(k)(i) = numElementsL(k);
      
        elMapItem = grid3d->getElements().getRowMapL(i);
        elMapItem.setLid(numElementsL(k));
      
        globElements(k).push_back(elMapItem, elMapItem.getGid());
      }
    }
  }
  
  //Assign new global ids--------------------------------------------------------------------------
  pVectGlobalManip<UInt,PMAPTYPE> nodeManip(commDev);
  pVectGlobalManip<UInt,PMAPTYPE>   elManip(commDev);
  
  for(UInt k=1; k <= numBlocks; ++k)
  {
    nodeManip.buildGlobalNumbering(globVertices(k));
    elManip.buildGlobalNumbering(globEdges(k));
    elManip.buildGlobalNumbering(globFaces(k));
    elManip.buildGlobalNumbering(globElements(k));
    
    numVerticesG(k) = nodeManip.sizeG(globVertices(k));
    numEdgesG(k)    = elManip.sizeG(globEdges(k));
    numFacesG(k)    = elManip.sizeG(globFaces(k));
    numElementsG(k) = elManip.sizeG(globElements(k));
  }
    
  //Count the dofs---------------------------------------------------------------------------------
  for(UInt k=1; k <= numBlocks; ++k)
  {
    numVertexDofsL(k) = FETYPE::dofPerVertex * numVerticesL(k);
    numEdgesDofsL(k)  = FETYPE::dofPerEdge   * numEdgesL(k);
    numFacesDofsL(k)  = FETYPE::dofPerFace   * numFacesL(k);
    numVolumeDofsL(k) = FETYPE::dofPerVolume * numElementsL(k);
    
    numVertexDofsG(k) = FETYPE::dofPerVertex * numVerticesG(k);
    numEdgesDofsG(k)  = FETYPE::dofPerEdge   * numEdgesG(k);
    numFacesDofsG(k)  = FETYPE::dofPerFace   * numFacesG(k);
    numVolumeDofsG(k) = FETYPE::dofPerVolume * numElementsG(k);
  }
  
  //Total number of dofs and lists-----------------------------------------------------------------
  numDofsTotL = 0.0;
  numDofsTotG = 0.0;
  
  for(UInt k=1; k <= numBlocks; ++k)
  {
    numDofsL(k) = numVertexDofsL(k) 
                + numEdgesDofsL(k)
                + numFacesDofsL(k)
                + numVolumeDofsL(k);
    
    numDofsG(k) = numVertexDofsG(k)
                + numEdgesDofsG(k)
                + numFacesDofsG(k)
                + numVolumeDofsG(k);
      
    numDofsTotL += numDofsL(k);
    numDofsTotG += numDofsG(k);
  }
  
  sizeListL = numDofsTotL * traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
  sizeListG = numDofsTotG * traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
  
  //Build dof map----------------------------------------------------------------------------------
  PMAPTYPE mapItem;
  UInt index      = 1;
  UInt offsetTotG = 0;
  
  dofMap.resize(numDofsTotL);
  
  for(UInt k=1; k <= numBlocks; ++k)
  {
    offsetL(k) = index - 1; //Update offset
      
    for(UInt lev=1; lev <= FETYPE::dofPerVertex; ++lev)  //Vertices
    {
      for(UInt i=1; i <= numVerticesL(k); ++i)
      {
        mapItem = globVertices(k).getMapL(i);
        mapItem.setLid(index);
        mapItem.setGid(mapItem.getGid()
                       + ((lev-1) * numVerticesG(k))
                       + offsetTotG);
      
        dofMap(index) = mapItem;
        ++index;
      } 
    }
    
    for(UInt lev=1; lev <= FETYPE::dofPerEdge; ++lev)  //Edges
    {
      for(UInt i=1; i <= numEdgesL(k); ++i)
      {
        mapItem = globEdges(k).getMapL(i);
        mapItem.setLid(index);
        mapItem.setGid(mapItem.getGid()
                       + ((lev-1) * numEdgesG(k))
                       + numVertexDofsG(k)
                       + offsetTotG);
      
        dofMap(index) = mapItem;
        ++index;
      }
    }
    
    for(UInt lev=1; lev <= FETYPE::dofPerFace; ++lev)  //Faces
    {
      for(UInt i=1; i <= numFacesL(k); ++i)
      {
        mapItem = globFaces(k).getMapL(i);
        mapItem.setLid(index);
        mapItem.setGid(mapItem.getGid()
                       + ((lev-1) * numFacesG(k))
                       + numVertexDofsG(k) + numEdgesDofsG(k)
                       + offsetTotG);
      
        dofMap(k) = mapItem;
        ++index;
      }
    }
    
    for(UInt lev=1; lev <= FETYPE::dofPerVolume; ++lev)  //Volume
    {
      for(UInt i=1; i <= numElementsL(k); ++i)
      {
        mapItem = globElements(k).getMapL(i);
        mapItem.setLid(index);
        mapItem.setGid(mapItem.getGid()
                       + ((lev-1) * numElementsG(k))
                       + numVertexDofsG(k) + numEdgesDofsG(k) + numFacesDofsG(k)
                       + offsetTotG);
      
        dofMap(k) = mapItem;
        ++index;
      }
    }
    
    offsetG(k)  = offsetTotG; //Update offset
    offsetTotG += numDofsG(k); 
  }
  
  //Build list map---------------------------------------------------------------------------------
  listMap.resize(sizeListL);
  
  UInt lid, gid;
  
  for(UInt i=1; i <= dofMap.size(); ++i)
  {   
    for(UInt K=1; K <= block; ++K)
    {
      mapItem = dofMap(i);
      
      lid = (block * (mapItem.getLid() - 1) + K) * (ORDER == dms3d_vectMajor)  + (mapItem.getLid() + (K-1) * numDofsTotL) * (ORDER == dms3d_componentMajor);
      gid = (block * (mapItem.getGid() - 1) + K) * (ORDER == dms3d_vectMajor)  + (mapItem.getGid() + (K-1) * numDofsTotG) * (ORDER == dms3d_componentMajor);
      
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

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
bool
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
isActive(const DOFCARD & dofCard) const
{
  assert(dofCard.getLocalElId() <= elementsToBlock.size());
  return(elementsToBlock(dofCard.getLocalElId()) != 0);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofL - 1) + K)      * (ORDER == dms3d_vectMajor) +
          (dofL  + (K-1) * numDofsTotL) * (ORDER == dms3d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofG - 1) + K)      * (ORDER == dms3d_vectMajor) +
          (dofG  + (K-1) * numDofsTotG) * (ORDER == dms3d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
mapDofL(const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(isActive(dofCard));
  assert(dofCard.getLocalElId() >= 1);
  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
  
  UInt oldId   = 0;
  UInt blockId = elementsToBlock(dofCard.getLocalElId());
  
  switch(dofCard.getGeoType())
  {
    case VERTEX :
      assert(dofCard.getLevel()   >= 1);  assert(dofCard.getLevel()   <= FETYPE::dofPerVertex);
      assert(dofCard.getLocalId() >= 1);  assert(dofCard.getLocalId() <= GEOSHAPE::numVertices);
      
      oldId = grid3d->getElements().getCid_LL(dofCard.getLocalElId(), dofCard.getLocalId());
      return( newVertexLid(blockId)(oldId) 
           + (dofCard.getLevel() - 1) * numVerticesL(blockId)
           +  offsetL(blockId));
    
    case EDGE   :
      assert(dofCard.getLevel()   >= 1);  assert(dofCard.getLevel()   <= FETYPE::dofPerEdge);
      assert(dofCard.getLocalId() >= 1);  assert(dofCard.getLocalId() <= GEOSHAPE::numEdges);
      
      oldId = connectGrid3d->getElementToEdge(dofCard.getLocalElId(), dofCard.getLocalId());
      return( newEdgeLid(blockId)(oldId) 
           + (dofCard.getLevel() - 1) * numEdgesL(blockId)
           +  numVertexDofsL(blockId)
           +  offsetL(blockId));
    
    case FACE   :
      assert(dofCard.getLevel()   >= 1);  assert(dofCard.getLevel()   <= FETYPE::dofPerFace);
      assert(dofCard.getLocalId() >= 1);  assert(dofCard.getLocalId() <= GEOSHAPE::numFaces);
      
      oldId = connectGrid3d->getElementToFace(dofCard.getLocalElId(), dofCard.getLocalId());
      return( newFaceLid(blockId)(oldId) 
           + (dofCard.getLevel() - 1) * numFacesL(blockId)
           +  numVertexDofsL(blockId) + numEdgesDofsL(blockId)
           +  offsetL(blockId));
    
    case VOLUME :
      assert(dofCard.getLevel()   >= 1);  assert(dofCard.getLevel() <= FETYPE::dofPerVolume);
      assert(dofCard.getLocalId() == 1);
      
      oldId = dofCard.getLocalElId();
      return( newElementLid(blockId)(oldId)
           + (dofCard.getLevel() - 1) * numElementsL(blockId)
           +  numVertexDofsL(blockId) + numEdgesDofsL(blockId) + numFacesDofsL(blockId)
           +  offsetL(blockId));
      
    default :
      assert(1==2);
   }
  
  return(0); 
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
mapDofG(const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(isActive(dofCard));
  
  return(dofMap(mapDofL(dofCard)).getGid());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofL = mapDofL(dofCard);
  return(dofToListL(dofL,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofG = mapDofG(dofCard);
  return(dofToListG(dofG,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::DOFMAP &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getDofMap() const
{
  return(dofMap);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::LISTMAP &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getListMap() const
{
  return(listMap);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumDofsL() const
{
  return(numDofsTotL);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumDofsG() const
{
  return(numDofsTotG);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getSizeListL() const
{
  return(sizeListL);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getSizeListG() const
{
  return(sizeListG);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::FECARDS &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getFeCards() const
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::FECARDS &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getFeCards()
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<communicator> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getCommDev() const
{
  assert(commDevLoaded);
  return(commDev);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::MESH3D> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getGrid3d() const
{
  assert(geometryLoaded);
  return(grid3d);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::CONNECT3D> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getConnectGrid3d() const
{
  assert(geometryLoaded);
  return(connectGrid3d);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::OPTIONS> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getOptions() const
{
  assert(optionsLoaded);
  return(options);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumBlocks() const
{
  return(options->getNumBlocks());
}
    
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::SVECT_UINT &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNewVertexLid(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(newVertexLid(k));
}
 
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::SVECT_UINT &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNewEdgeLid(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(newEdgeLid(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::SVECT_UINT &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNewFaceLid(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(newFaceLid(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::SVECT_UINT &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNewElementLid(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(newElementLid(k));
}
    
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::PVECT_UINT &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getGlobVertices(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(globVertices(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::PVECT_UINT &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getGlobEdges(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(globEdges(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::PVECT_UINT &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getGlobFaces(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(globFaces(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::PVECT_UINT &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getGlobElements(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(globElements(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumVerticesL(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numVerticesL(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumEdgesL(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numEdgesL(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumFacesL(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numFacesL(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumElementsL(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numElementsL(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumVerticesG(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numVertexDofsG(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumEdgesG(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numEdgesG(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumFacesG(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numFacesG(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumElementsG(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numElementsG(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumVertexDofsL(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numVertexDofsL(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumEdgesDofsL(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numEdgesDofsL(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumFacesDofsL(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numFacesDofsL(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumVolumeDofsL(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numVolumeDofsL(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumVertexDofsG(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numVertexDofsG(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumEdgesDofsG(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numEdgesDofsG(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumFacesDofsG(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numFacesDofsG(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getNumVolumeDofsG(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(numVolumeDofsG(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getOffsetL(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(offsetL(k));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_disjointMode>::
getOffsetG(const UInt & k) const
{
  assert(k <= options->getNumBlocks());
  return(offsetG(k));
}

template<typename F, typename D, dms3d_order O>
ostream & operator<<(ostream & f, const dofMapStatic3d<F,D, O,dms3d_disjointMode> & V)
{
  f << "numVertices   : " << V.numVerticesL   << " " << V.numVerticesG << endl;
  f << "numEdges      : " << V.numEdgesL      << " " << V.numEdgesG << endl;
  f << "numFaces      : " << V.numFacesL      << " " << V.numFacesG << endl;
  f << "numElements   : " << V.numElementsL   << " " << V.numElementsG << endl << endl;
  
  f << "numVertexDofs : " << V.numVertexDofsL << " " << V.numVertexDofsG << endl;
  f << "numEdgesDofs  : " << V.numEdgesDofsL  << " " << V.numEdgesDofsG << endl;
  f << "numFacesDofs  : " << V.numFacesDofsL  << " " << V.numFacesDofsG << endl;
  f << "numVolumeDofs : " << V.numVolumeDofsL << " " << V.numVolumeDofsG << endl << endl;
  
  f << "numTotalDofs  : " << V.numDofsL       << " " << V.numDofsG << endl;
  f << "sizeList      : " << V.sizeListL      << " " << V.sizeListG << endl;
  
  return(f);
}

#endif
